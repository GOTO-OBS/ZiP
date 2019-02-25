import subprocess
import inspect
from random import randint
import math
import glob
import os
import shutil
import multiprocessing
import multiprocessing.pool
from itertools import product
import sys
import time

import sep
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D as cut
from numpy import ndarray
from matplotlib.font_manager import FontProperties
from multiprocessing import Pool
from astropy.convolution import convolve
from astropy.io import ascii
from scipy import stats
import numpy as np
import scipy.interpolate as interp
from scipy import ndimage
try:
    import pyfftw.interfaces.numpy_fft as fft
except:
    print('For a faster performance install pyfftw')
    print(' ')
    print(' ')
    print(' ')
    import numpy.fft as fft  #Use if you don't have pyfftw

###from zogyp.zip import run_ZOGY

#Get correct call for sextractor
try:
    sex = 'sextractor'
    subprocess.call(sex, stderr=subprocess.DEVNULL)
except:
    sex = 'sex'
    subprocess.call(sex, stderr=subprocess.DEVNULL)

root_config = os.path.dirname(os.path.realpath(__file__))+'/configfls'
#####################################################################################################################
class NoDaemonProcess(multiprocessing.Process):
    """ 
    These two classes allow the daemons to have children... how nice
    Basically, allows a function in a pool to use its own pool
    Same as in zip
    """
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class NoDaemonProcessPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess
#####################################################################################################################

#####################################################################################################################
def get_fratio(psfcat_sci, psfcat_ref, sexcat_sci, sexcat_ref):
    """
    Function that takes in output catalogs of stars used in the PSFex runs on the new and the ref image, and 
    returns the arrays with pixel coordinates (!) x, y (in the new frame) and fratios (flux ratio) for 
    the matching stars. 

    In addition, it provides the difference in stars' RAs and DECs in arcseconds between the two catalogs.
    (This step is broken currently)!
    """
    def readcat (psfcat):
        table = ascii.read(psfcat, format='sextractor')
        number = table['SOURCE_NUMBER']
        x = table['X_IMAGE']
        y = table['Y_IMAGE']
        norm = table['NORM_PSF']
        return(number, x, y, norm)

    # read in psfcat_sci
    number_sci, x_sci, y_sci, norm_sci = readcat(psfcat_sci)
    # read in psfcat_ref
    number_ref, x_ref, y_ref, norm_ref = readcat(psfcat_ref)

    def xy2radec (number, sexcat):
        # read the Source Extractor fits table
        with fits.open(sexcat) as hdulist:
            data = hdulist[2].data
            ra_sex = data['ALPHAWIN_J2000']
            dec_sex = data['DELTAWIN_J2000']
            fwhm_sex = data['FWHM_IMAGE']
            Elon_sex = data['ELONGATION']
            X_sex = data['X_IMAGE']
            Y_sex = data['Y_IMAGE']
            # record in ra, dec for each source
            ra = []
            dec = []
            Y_pos = []
            X_pos = []
            FWHM = []
            ELON = []
            for n in number:
                ra.append(ra_sex[n-1])
                dec.append(dec_sex[n-1])
                X_pos.append(X_sex[n-1])
                Y_pos.append(Y_sex[n-1])
                FWHM.append(fwhm_sex[n-1])
                ELON.append(Elon_sex[n-1])
        return(np.array(ra), np.array(dec), np.array(X_pos), np.array(Y_pos), 
               np.array(FWHM), np.array(ELON))

    # get ra, dec in pixel coords
    #PARALLEL#
    ra_sci, dec_sci, X_pos_sci, Y_pos_sci, FWHM_sci, ELON_sci = xy2radec(number_sci, sexcat_sci)
    ra_ref, dec_ref, X_pos_ref, Y_pos_ref, FWHM_ref, ELON_ref = xy2radec(number_ref, sexcat_ref)

    # now find matching sources, this step needs improving
    x_sci_match = []
    y_ref_match = []
    dx = []
    dy = []
    dra_match = []
    ddec_match = []
    fratio = []
    nmatch = 0
    for i_sci in range(len(x_sci)):
        # calculate distance to ref objects
        dist = np.sqrt((X_pos_sci[i_sci] - X_pos_ref)**2 + (Y_pos_sci[i_sci] - Y_pos_ref)**2)
        # minimum distance and its index
        dist_min, i_ref = np.amin(dist), np.argmin(dist)

        if dist_min <5.: #This min distance is dependant on your registrtion. 
                          #The less confident you are in your registration the bigger it needs to be.
            nmatch += 1
            select = max(FWHM_sci[i_sci], FWHM_ref[i_ref])
            x_sci_match.append(x_sci[i_sci])
            y_ref_match.append(y_sci[i_sci])
            dx.append(select)
            dy.append(select)
            # append ratio of normalized counts to fratios
            fratio.append(norm_sci[i_sci] / norm_ref[i_ref])

    return(np.array(x_sci_match), np.array(y_ref_match), np.array(fratio), 
           np.array(dx), np.array(dy))

#####################################################################################################################

#####################################################################################################################
def rdfits(image):
    """ 
    Compressed photo images are not compatible with SExtractor and SWARP
    This function remakes the fits files so they are
    """
    try:
        hdulist =  fits.open(image)
        header= hdulist[1]._header
        data = hdulist[1].data
     
    except:
        hdulist =  fits.open(image)
        header= hdulist[0].header
        data = hdulist[0].data
  
    image_RD = image.replace('.fits', '_RD.fits')

    fits.writeto(image_RD, data, header, overwrite=True)
    return(image_RD)

#####################################################################################################################


#####################################################################################################################
def clean_norm_psf (psf_ar, clean_fact = 0.25):
    """
    Normalises the psf required for ZOGY, will remove any values smaller than the clean factor 
    to avoid clutter by precision
    """
    ysize, xsize = psf_ar.shape
    assert ysize == xsize

    hsize = ysize/2

    if xsize % 2 == 0:
        x = np.arange(-hsize,hsize)
    else:
        x = np.arange(-hsize,hsize)

    xx, yy = np.meshgrid(x,x, sparse=True)
    psf_ar[(xx**2+yy**2)>hsize**2] = 0

    if clean_fact != 0:
        mask_clean = (psf_ar < (np.amax(psf_ar)*clean_fact))
        psf_ar[mask_clean]=0

    psf_ar_norm = psf_ar / np.sum(psf_ar)

    return(psf_ar_norm)
#####################################################################################################################


#####################################################################################################################
def get_psf(image):
    """
    Finds the psf using PSFex.
    """

    r = rdfits(image)
    sexcat = image.replace('.fits', '_PSFCAT.fits')
    sexcat = sexcat.replace('./output/','')
    print('Making PSF Catalog') 
    subprocess.call([sex, '-c', root_config+'/default.sex', r, '-CATALOG_NAME', sexcat,
                    '-PARAMETERS_NAME', root_config+'/default.param', '-STARNNW_NAME', 
                     root_config+'/default.nnw','-FILTER_NAME', root_config+'/default.conv'
                     ], stderr=subprocess.DEVNULL)

    outcat = sexcat.replace('_PSFCAT.fits', '.psfexcat')

    talk2 = ['psfex', sexcat, '-c', root_config+'/psfex.conf', '-OUTCAT_NAME', outcat]
    print('Modelling PSF '+image)
    subprocess.call(talk2, stderr=subprocess.DEVNULL)
    with fits.open(sexcat.replace('.fits','.psf')) as hdulist:
        header = hdulist[1].header
        data = hdulist[1].data
    dat = data[0][0][:]
    subprocess.call(['rm', r])
    return(dat, header, sexcat, outcat)

#####################################################################################################################


#####################################################################################################################
def psf_map(dat, header, const, xl, yl, xc, yc, slices, NUMBER):
    """
    Takes a slice of the data and maps the psf to a kernal for that slice (allows for a varying psf)
    """
    polzero1 = header['POLZERO1']
    polzero2 = header['POLZERO2']
    polscal1 = header['POLSCAL1']
    polscal2 = header['POLSCAL2']
    poldeg = header['POLDEG1']
    psf_samp = header['PSF_SAMP']


    psf_size_config = header['PSFAXIS1']
    psf_size = np.int(np.ceil(psf_size_config * psf_samp))
    if psf_size % 2 == 0:
        psf_size += 1
    psf_samp_update = float(psf_size) / float(psf_size_config)

    ysize_fft = yl
    xsize_fft = xl

    xcenter_fft, ycenter_fft = xsize_fft/2, ysize_fft/2


    psf_ima_center = np.zeros((ysize_fft,xsize_fft), dtype='float32')
    # [psf_ima_shift] is [psf_ima_center] shifted - this is
     # the input PSF image needed in the zogy function
    psf_ima_shift = np.zeros((ysize_fft,xsize_fft), dtype='float32')


    x = (xc - polzero1) / polscal1
    y = (yc - polzero2) / polscal2

    if slices == 1:
        psf = dat[0]
    else:
        if poldeg == 2:
            psf = dat[0] + dat[1] * x + dat[2] * x**2 + dat[3] * y + dat[4] * x * y + dat[5] * y**2
        elif poldeg == 3:
            psf = dat[0] + dat[1] * x + dat[2] * x**2 + dat[3] * x**3 + \
            dat[4] * y + dat[5] * x * y + dat[6] * x**2 * y + \
            dat[7] * y**2 + dat[8] * x * y**2 + dat[9] * y**3

    psf_ima_resized = ndimage.zoom(psf, psf_samp_update)
    psf_ima_resized_norm = clean_norm_psf(psf_ima_resized, const)
    psf_hsize = math.floor(psf_size/2)

    ind = [slice(int(ycenter_fft-psf_hsize), int(ycenter_fft+psf_hsize+1)),
           slice(int(xcenter_fft-psf_hsize), int(xcenter_fft+psf_hsize+1))]

    psf_ima_center[tuple(ind)] = psf_ima_resized_norm

    # perform fft shift
    psf_ima_shift = fft.fftshift(psf_ima_center)
    return(psf_ima_shift)
#####################################################################################################################


#####################################################################################################################
def data_chunks(data, xslice, yslice):
    """
    slices image into sub_images returns centre co-ords as well
    """
    sub_img = []
    centres = []
    x = np.linspace(0, data.shape[0], xslice+1)
    y = np.linspace(0, data.shape[1], yslice+1)

    for i in range(len(x)-1):
        for j in range(len(y)-1):
            sub_img.append(data[int(x[i]):int(x[i+1]), int(y[j]):int(y[j+1])])
            centres.append([(int(x[i])+int(x[i+1]-1))/2, (int(y[j]) + int(y[j+1]))/2])

    return(sub_img, centres)
#####################################################################################################################


#####################################################################################################################
def restitcher(data, new_cut_data, xslice, yslice):
    """
    stitch the image back together
    """
    #Empty array to fill with new cut data
    data_empty = np.zeros((data.shape[0], data.shape[1]), dtype='float32') 

    x = np.linspace(0, data.shape[0], xslice+1)
    y = np.linspace(0, data.shape[1], yslice+1)
    for i in range(len(x)-1):
        for j in range(len(y)-1):
            data_empty[int(x[i]):int(x[i+1]), int(y[j]):int(y[j+1])] =  new_cut_data[(len(y)-1)*i + j]

    return(data_empty)
#####################################################################################################################


#####################################################################################################################
def chop_kern(data, psf_dat, psf_hed, xslice, yslice, clean_const=0.65):
    """
    A handy third party function to find the PSF for specific segments of the input image returns the sub images and 
    corresponding PSF
    """
    slices = xslice * yslice
    psf = []
    sub_img, cents = data_chunks(data, xslice, yslice)
    for i in range(len(sub_img)):
        psf.append(psf_map(psf_dat, psf_hed, clean_const, sub_img[i].shape[1], sub_img[i].shape[0], cents[i][1], cents[i][0], slices, i))

    return(sub_img, psf)
#####################################################################################################################

#####################################################################################################################
def coad_func(j,base_psfcat, base_sexcat):
    # j is the jth coadded image
    psf_dat, psf_hed, sexcat2, psfcat2 = get_psf(j)
    x_fratio, y_fratio, fratio, dx, dy = get_fratio(base_psfcat, psfcat2, base_sexcat, sexcat2) 
    R_j = fits.getdata(j)

    R_j_hat = fft.fft2(R_j)

    _, PSF = chop_kern(R_j, psf_dat, psf_hed,1,1)
    P_j =  PSF[0]
    #print(P_j.shape)
    P_j_hat = fft.fft2(P_j)
    P_j_hat_bar = np.conj(P_j_hat)
    
    F_j = 1.0/np.median(fratio) #Zeropoint flux
    sig_j_2 = np.std(R_j)**2

    Nomin = (F_j/sig_j_2) * P_j_hat_bar * R_j_hat
    Denom = (F_j**2/sig_j_2) * np.abs(P_j_hat**2)
    subprocess.call(['rm', sexcat2, psfcat2, sexcat2.replace('.fits', '.psf')])
    return(Nomin, Denom)
#####################################################################################################################


#####################################################################################################################
def prop_coad(ref_dir, make_fits=False):
    """
    Proper coaddition function
    """
 
    if len(ref_dir)==1:
        F  = glob.glob(ref_dir[0]+'/*.fits') # collect images you want to coadd
        print(F)
    else:
        F = ref_dir
        print(F)

    psf_dat, psf_hed, sexcat1, psfcat1 = get_psf(F[0])  
    pool = multiprocessing.Pool(len(F)-1)
    tmp_array = pool.starmap(coad_func, [(F[i+1] , psfcat1, sexcat1) for i in range(len(F)-1)])

    Nomin = sum(x[0] for x  in tmp_array)
    Denom = sum(x[1] for x in tmp_array)

    Denom = np.sqrt(Denom)
    if np.any(Denom==0):
        print('ZEROS')

    R_hat = Nomin/Denom
    R = np.real(fft.ifft2(R_hat))
    subprocess.call(['rm', sexcat1, psfcat1, sexcat1.replace('.fits', '.psf')])
    if make_fits == True:
        hed = fits.getheader(F[0])
        hed['COMMENT'] = 'ZO coaddition from ZiP'
        hed['COMMENT'] = 'List of stacked fits'
        hed['COMMENT'] =  ', '.join(F)
        fits.writeto(F[0].replace('.fits', '_COAD.fits'), R, hed, overwrite = True)
        return(F[0].replace('.fits', '_COAD.fits'), R)
    else:
        return(R)
#####################################################################################################################
def med_comb(ref_dir, make_fits=False):
 
    if len(ref_dir)==1:
        F  = glob.glob(ref_dir[0]+'/*.fits') # collect images you want to coadd
        print(F)
    else:
        F = ref_dir
        print(F)
    data = [] 
    for fil in F:
        n = rdfits(fil)
        hed = fits.getheader(n)
        data_fil = fits.getdata(n)
        data.append(data_fil)
        subprocess.call(['rm', n])

    stack = np.median(data, axis=0)
    if make_fits == True:
        hed['COMMENT'] = 'Median Combined from ZiP'
        hed['COMMENT'] = 'List of stacked fits'
        hed['COMMENT'] =  ', '.join(F)
        fits.writeto(F[0].replace('.fits', '_MED_COMB.fits'), stack, hed, overwrite = True)
        return(F[0].replace('.fits', '_MED_COMB.fits'), stack)  
    else:
        return(stack)
#####################################################################################################################
