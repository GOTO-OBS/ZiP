import math
import glob
import sep
import os
import shutil
import multiprocessing
import multiprocessing.pool
from multiprocessing import Pool
import sys
import inspect
import time
from random import randint
import subprocess
from itertools import product

from astropy.nddata.utils import Cutout2D as cut
from astropy.convolution import convolve
from astropy.io import fits
from astropy.io import ascii
from astropy.wcs import WCS
from numpy import ndarray
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from scipy import stats
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

from zogyp.spali2 import spalipy

#use the correct call for sextractor
try:
    sex = 'sextractor'
    subprocess.call(sex, stderr=subprocess.DEVNULL)
except:
    sex = 'sex'
    subprocess.call(sex, stderr=subprocess.DEVNULL) 



root_config = os.path.dirname(os.path.realpath(__file__))+'/configfls' #path to config files

def config_loc(): #Quick function to find config files
    print(root_config)
    return(root_config)

#####################################################################################################################
class NoDaemonProcess(multiprocessing.Process):
    """ 
    These two classes allow the daemons to have children... how nice
    Basically, allows a function in a pool to use its own pool
    """
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)


# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class NoDaemonProcessPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess
#####################################################################################################################


#####################################################################################################################
def rdfits(image):
    """ 
    Compressed fits images are not compatible with SExtractor
    This function remakes the fits files so they are
    """
    #Depending on compression type
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

    sexcat = image.replace('.fits', '_PSFCAT.fits')
    sexcat = sexcat.replace('./Zoutput/','') 
    talk = [sex ,image,'-c',root_config+'/default.sex' , '-CATALOG_NAME' , sexcat,
           '-PARAMETERS_NAME', root_config+'/default.param', '-STARNNW_NAME', 
           root_config+'/default.nnw', '-FILTER_NAME', root_config+'/default.conv']

    print('Making PSF catalog')
    subprocess.call(talk, stderr=subprocess.DEVNULL)

    outcat = sexcat.replace('_PSFCAT.fits', '.psfexcat')
 
    print('Modelling PSF for '+image)
    talk2 = ['psfex', sexcat, '-c', root_config+'/psfex.conf', '-OUTCAT_NAME', outcat]
    subprocess.call(talk2,stderr=subprocess.DEVNULL)
    with fits.open(sexcat.replace('.fits','.psf')) as hdulist:
        header = hdulist[1].header
        data = hdulist[1].data

    dat = data[0][0][:]
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
    # [psf_ima_shift] is [psf_ima_center] FFT shifted - this is
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
    data_empty = np.zeros((data.shape[0], data.shape[1]), dtype='float32') #empty array to fill with new cut data
    x = np.linspace(0, data.shape[0], xslice+1)
    y = np.linspace(0, data.shape[1], yslice+1)
    for i in range(len(x)-1):
        for j in range(len(y)-1):
            data_empty[int(x[i]):int(x[i+1]), int(y[j]):int(y[j+1])] =  new_cut_data[(len(y)-1)*i + j]

    return(data_empty)
#####################################################################################################################


#####################################################################################################################
def chop_kern(data, psf_dat, psf_hed, xslice, yslice, clean_const=0.5):
   """
   A handy third party function to find the PSF for specific segments of the input image returns the sub images and 
   corresponding PSF
   """
   slices = xslice * yslice
   psf = []
   sub_img, cents = data_chunks(data, xslice, yslice)
   for i in range(len(sub_img)):
       psf.append(psf_map(psf_dat, psf_hed, clean_const, sub_img[i].shape[1], 
                  sub_img[i].shape[0], cents[i][1], cents[i][0], slices, i))

   return(sub_img, psf)

#####################################################################################################################


####################################################################################################################
def get_fratio(psfcat_sci, psfcat_ref, sexcat_sci, sexcat_ref):
    """
    Function that takes in output catalogs of stars used for the PSFex runs, and 
    returns the arrays with pixel coordinates x, y and fratios (flux ratio not 
    the statistic) for the matching stars. 

    In addition, it provides the difference in stars' RAs and DECs in arcseconds 
    between the two catalogs. (This step is broken currently!) It also determines 
    the error in position due to the Ful Width HalF Maximum (FWHM).
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
        return(np.array(ra), np.array(dec), np.array(X_pos), 
               np.array(Y_pos), np.array(FWHM), np.array(ELON))

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

####################################################################################################################


####################################################################################################################
def imprep(sci_im, ref_im, chopx, chopy, inject = False, lineup = False):
    """
    If not aligned, will realign sci to match ref (as Ref should have better astometry)
    Will chop up the images
    Can inject fake sources into the science image (for testing)
    Returns fits files of the aligned, chopped, and source injected images
    """
    # F = glob.glob('./Zoutput/*')
    # for fil in F:
    #  subprocess.call(['rm', fil])

    # This part formats the files to guarentee the other programs can use them#
    name1 = rdfits(sci_im)
    headerA = fits.getheader(name1)
    name2 = rdfits(ref_im)

    hdulist = fits.open(name2)
    ref_data = hdulist[0].data
    header= hdulist[0].header
    try:
        W = WCS(header) #World coord system
    except:
        W = 0
        print('no compatible WCS')

    if lineup == True: #Register step
        name_new = spalipy(name1, name2) #Spalipy remap
        subprocess.call(['rm', name1, name2])
        hdulist2 = fits.open(name_new)
    else:
        hdulist2 = fits.open(name1)
        subprocess.call(['rm', name1, name2])

    sci_data = hdulist2[0].data

    X = int(header['NAXIS1'])
    Xfrac = X/chopx
    xcuts = math.ceil(X/Xfrac) # number of cuts in the x plane
    Y = int(header['NAXIS2'])
    Yfrac = Y/chopy
    ycuts= math.ceil(Y/Yfrac) #number of cuts in the y plane

    CC = [] #Centre Co-ordinates of each cut

    Xlin = np.linspace(0, X, xcuts+1)
    XC = [] #Centre of X points

    for i in range(len(Xlin)-1):
        XC.append( int((Xlin[i] + Xlin[i+1])/2) )

    Ylin = np.linspace(0, Y, ycuts+1)
    YC = [] #Centre of Y points

    for i in range(len(Ylin)-1):
        YC.append( int((Ylin[i] + Ylin[i+1])/2) )

    for i in range (len(XC)):
        for j in range(len(YC)):
            CC.append([XC[i], YC[j]]) #Centre Co-ords for each slice

    for i in range(len(CC)):
        if inject == True:
            datainj = sci_data
            print('make comprehensive injection function')
            sci_data = datainj

        if W != 0:
            #Includes WCS from ref, should be same for both images if 
             #they are alligned
            CUT = cut(ref_data, (CC[i][0], CC[i][1]), (Yfrac, Xfrac), W)
            CUT2 = cut(sci_data, (CC[i][0], CC[i][1]), (Yfrac, Xfrac), W)
        else:
            CUT = cut(ref_data, (CC[i][0], CC[i][1]), (Yfrac, Xfrac))
            CUT2 = cut(sci_data, (CC[i][0], CC[i][1]), (Yfrac, Xfrac))

        hdu = fits.PrimaryHDU(CUT.data, header=CUT.wcs.to_header())
        hdu.writeto('Zoutput/ref_cut%s.fits' %(i+1), overwrite=True)
        hdu = fits.PrimaryHDU(CUT2.data, header=CUT2.wcs.to_header())
        try:
            hdu.header['DATE2'] = headerA['DATE-OBS']
            hdu.header['ORIG-ID']= headerA['RUN-ID']
            hdu.header['INST']= headerA['INSTRUME']
            hdu.header['REF-ID'] = header['RUN-ID']
            hdu.header['INST2']= header['INSTRUME']
            hdu.writeto('Zoutput/sci_cut%s.fits' %(i+1), overwrite=True)
        except:
            hdu.writeto('Zoutput/sci_cut%s.fits' %(i+1), overwrite=True)
    return(None)
#####################################################################################################################


#####################################################################################################################
def ZOGY(R,N,Pr,Pn,sr,sn,fr,fn,Vr,Vn,dx,dy):
    """
    Optimal image subtraction in a pythonic layout! 
    Where the magic happens. Algorithm from
    Zackay et al. (2016)
    """
    R_hat = fft.fft2(R)

    N_hat = fft.fft2(N)

    Pn_hat = fft.fft2(Pn)

    Pn_hat2_abs = np.abs(Pn_hat**2)

    Pr_hat = fft.fft2(Pr)

    Pr_hat2_abs = np.abs(Pr_hat**2)

    sn2 = sn**2

    sr2 = sr**2

    fn2 = fn**2

    fr2 = fr**2

    fD = fr*fn / np.sqrt(sn2*fr2+sr2*fn2)

    denom = sn2*fr2*Pr_hat2_abs + sr2*fn2*Pn_hat2_abs
    if np.any(denom==0):
        print('There are zeros!')

    D_hat = (fr*Pr_hat*N_hat - fn*Pn_hat*R_hat) / np.sqrt(denom)
    D = np.real(fft.ifft2(D_hat)) / fD

    P_D_hat = (fr*fn/fD) * (Pr_hat*Pn_hat) / np.sqrt(denom)

    S_hat = fD*D_hat*np.conj(P_D_hat)
    S = np.real(fft.ifft2(S_hat))

    kr_hat = fr*fn2*np.conj(Pr_hat)*Pn_hat2_abs / denom
    kr = np.real(fft.ifft2(kr_hat))
    kr2 = kr**2
    kr2_hat = fft.fft2(kr2)

    kn_hat = fn*fr2*np.conj(Pn_hat)*Pr_hat2_abs / denom
    kn = np.real(fft.ifft2(kn_hat))
    kn2 = kn**2
    kn2_hat = fft.fft2(kn2)

    Vr_hat = fft.fft2(Vr)
    Vn_hat = fft.fft2(Vn)

    VSr = np.real(fft.ifft2(Vr_hat*kr2_hat))
    VSn = np.real(fft.ifft2(Vn_hat*kn2_hat))

    dx2 = dx**2
    dy2 = dy**2
    # and calculate astrometric variance
    Sn = np.real(fft.ifft2(kn_hat*N_hat))
    dSndy = Sn - np.roll(Sn,1,axis=0)
    dSndx = Sn - np.roll(Sn,1,axis=1)
    VSn_ast = dx2 * dSndx**2 + dy2 * dSndy**2

    Sr = np.real(fft.ifft2(kr_hat*R_hat))
    dSrdy = Sr - np.roll(Sr,1,axis=0)
    dSrdx = Sr - np.roll(Sr,1,axis=1)
    VSr_ast = dx2 * dSrdx**2 + dy2 * dSrdy**2

    # and finally S_corr
    V_S = VSr + VSn
    V_ast = VSr_ast + VSn_ast
    V = V_S + V_ast

    # make sure there's no division by zero
    S_corr = np.copy(S)
    S_corr[V>0] /= np.sqrt(V[V>0])

    F_S =  np.sum((fn2*Pn_hat2_abs*fr2*Pr_hat2_abs) / denom)
    F_S /= R.size

    alpha = S / F_S
    alpha_std = np.zeros(alpha.shape)
    alpha_std[V_S>0] = np.sqrt(V_S[V_S>0]) / F_S

    return(D, S, S_corr, alpha, alpha_std)
####################################################################################################################

###################################################################################################################
def speedup_func(i, cdat2, cdat,dat, dat2,  psf2, psf, stdb2, stdb, 
                 f_ref, f_new, dx_full, dy_full):

    """
    A function to apply ZOGY in parallel
    """

    var_sci = (abs(cdat - np.median(dat))**2)
    var_ref = (abs(cdat2 - np.median(dat2))**2)

    data_D, data_S, data_Sc, fpsf, fpsf_std  = ZOGY(cdat2, cdat,  psf2, psf, np.median(stdb2), 
                                                    np.median(stdb), f_ref, f_new, var_ref, var_sci, 
                                                    dx_full, dy_full)

    return([data_D, data_S, data_Sc, i])



###################################################################################################################

####################################################################################################################
def finp(image, name, xslice, yslice, clean_sci, clean_ref, blackout):
    """ 
    Using all of the above, this function will find the psf of background subtracted data.
    Will then find the F-ratio, and pixel properties. The subtraction occurs producing D, S, 
    and Scorr images.
    """

    sci = image.replace('ref', 'sci')

    #Parallell PSF modelling
    with multiprocessing.Pool(2) as p:
        V1, V2 = p.map(get_psf, [sci, image])
    psf_dat, psf_hed, sexcat1, psfcat1 = V1 
    psf2_dat, psf2_hed, sexcat2, psfcat2 = V2 
    p.close()

    x_fratio, y_fratio, fratio, dx, dy = get_fratio(psfcat1, psfcat2, sexcat1, sexcat2)
    FM = np.median(fratio)

    fnum = image.replace('./Zoutput/ref_cut','')
    fnum2= fnum.replace('.fits','')

    f_new = 1.0
    f_ref = f_new/FM

    dx_full = np.median(dx)
    dy_full = np.median(dy)

    ######## Science data ########

    hdu = fits.open(sci)
    dat = hdu[0].data
    dat = dat.byteswap().newbyteorder()
    head = hdu[0].header
    W = WCS(head)
    bkg = sep.Background(dat, bw = 16, bh =16)
    stdb = bkg.rms()

    sub_dat = dat - bkg #bkg subtracted data

    ######### Ref data ########
    hdu2 = fits.open(image)
    dat2 = hdu2[0].data
    dat2 = dat2.byteswap().newbyteorder()
 
    if blackout == True: #remove data where there is no overlap
        dat2[dat==0] = 0
  	    		 
    head2 = hdu2[0].header
    W2 = WCS(head)
    bkg2 = sep.Background(dat2, bw =16, bh = 16)
    stdb2 = bkg2.rms()

    sub_dat2 = dat2 - bkg2 #bkg subtracted data
    

    cdat, psf = chop_kern(sub_dat, psf_dat, psf_hed, xslice, yslice, clean_sci)
    cdat2, psf2 = chop_kern(sub_dat2, psf2_dat, psf2_hed, xslice, yslice, clean_ref)


    data_D = [0]*len(cdat) #empty frames
    data_S = [0]*len(cdat)
    data_Sc = [0]*len(cdat)

    if multiprocessing.cpu_count() > 3: #faster to do in serial with fewer cores 
        with multiprocessing.Pool(len(cdat)) as p:
            NV_in = [stdb2, stdb, f_ref, f_new, dx_full, dy_full] # Non-varrying inputs
            OUT = p.starmap(speedup_func, [([i, cdat2[i], cdat[i], dat[i], dat2[i], psf2[i], 
                            psf[i]] + NV_in) for i in range(len(cdat))])

        j = 0
        for i in range(len(cdat)):
            if OUT[i][3] == j: #stitch data to empty frame
                data_D[j] = OUT[i][0]
                data_S[j] = OUT[i][1]
                data_Sc[j] = OUT[i][2]
                j += 1

    else:
        for i in range(len(cdat)):
            var_sci = (abs(cdat[i] - np.median(dat[i]))**2)
            var_ref = (abs(cdat2[i] - np.median(cdat2[i]))**2)
            data_D[i], data_S[i], data_Sc[i], fpsf, fpsf_std  = ZOGY(cdat2[i], cdat[i], psf2[i], psf[i], 
                                                                     np.median(stdb2), np.median(stdb), 
                                                                     f_ref, f_new, var_ref, var_sci, dx_full, 
                                                                     dy_full)


    D_img = restitcher(sub_dat, data_D, xslice, yslice)
    S_img = restitcher(sub_dat, data_S, xslice, yslice)
    Scorr_img = restitcher(sub_dat, data_Sc, xslice, yslice)

    hdu = fits.PrimaryHDU(D_img, header=head )
    hdu.writeto('./Zoutput/'+name+'_D'+fnum ,overwrite=True)

    hdu = fits.PrimaryHDU(S_img, header=head)
    hdu.writeto('./Zoutput/'+name+'_S'+fnum ,overwrite=True)

    hdu = fits.PrimaryHDU(Scorr_img, header=head)
    hdu.writeto('./Zoutput/'+name+'_Scorr'+fnum ,overwrite=True)

    subprocess.call(['rm', 'sci_cut%s_PSFCAT.psf' %(fnum2), 
                     'sci_cut%s.psfexcat' %(fnum2), 'sci_cut%s_PSFCAT.fits' %(fnum2)])
    subprocess.call(['rm', 'ref_cut%s_PSFCAT.psf' %(fnum2), 
                    'ref_cut%s.psfexcat' %(fnum2), 'ref_cut%s_PSFCAT.fits' %(fnum2)])

    return(None)
####################################################################################################################


###################################################################################################################
def run_ZOGY(sci_im, ref_im, outname = 'data', xslice = 1, yslice = 1, 
             align = False, clean_sci = 0.75,clean_ref = 0.75, 
             sub_imagex = 1, sub_imagey =1, blackout = False):


    """
    A compact function to make ZOGY callable from in a single line

    Only create sub_image if PSF variation over the field can't be modelled with a 3rd order polynomial
    or you are trying to reduce computation time. (The PSF model will degrade if this is used)
    """
 
     #      Make directory suitable       #
    #######################################
    if os.path.isdir('./Zoutput') == False:
        os.makedirs('./Zoutput')
        print('Output directory made!')
    else:
    	  for fil in glob.glob('./Zoutput/sci_cut*'):
    	      subprocess.call(['rm', fil])
    	  for fil in glob.glob('./Zoutput/ref_cut*'):
    	      subprocess.call(['rm', fil])
    ########################################

    x = multiprocessing.cpu_count()
    if x < 6:
        print('Serial version')
        ncores = x
        imprep(sci_im, ref_im, sub_imagex, sub_imagey, lineup = align)
        refs = glob.glob('./Zoutput/ref_cut*.fits')
        for SLICE in refs:
            finp(SLICE, outname, xslice, yslice, clean_sci, clean_ref, blackout)
    else:
        if x > 45:
            ncores = 15
        else:
            ncores = (int(x/3))
        print('Parallell version, using %s cores' %(ncores*3))
        imprep(sci_im, ref_im, sub_imagex, sub_imagey, lineup = align)
        refs = glob.glob('./Zoutput/ref_cut*.fits')
        p = NoDaemonProcessPool(processes = ncores)
        p.starmap(finp, [([R] + [outname, xslice, yslice, clean_sci, clean_ref, blackout])for R in refs])
        p.close()

    return(None)

