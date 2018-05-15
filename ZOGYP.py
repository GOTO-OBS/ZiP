from astropy.nddata.utils import Cutout2D as cut
from numpy import ndarray
from multiprocessing import Pool
from ZOGYfunc import fin #Incase you don't have enough cores
from astropy.convolution import convolve
from astropy.io import fits
from astropy.io import ascii
from astropy.wcs import WCS
from random import randint
import math
import glob
import numpy as np
import sep
import os
import subprocess
import scipy.interpolate as interp
import pyfftw.interfaces.numpy_fft as fft
#import numpy.fft as fft  #Use this if you don't have pyfftw (it's a little slower)
from scipy import ndimage
import shutil
import multiprocessing
import multiprocessing.pool
from itertools import product
import sys
import time

#####################################################################################################################
"""
Find transients inputs image number and if you want to find variables
"""

def fitra(NUM, var=False):

 #Make directories keep them neat#
 if os.path.exists('./tfinds') == False:
  os.makedirs('./tfinds')
 if os.path.exists('./tfinds/list'+str(NUM)):
  shutil.rmtree('./tfinds/list'+str(NUM))
  os.makedirs('./tfinds/list'+str(NUM))
 else:
  os.makedirs('./tfinds/list'+str(NUM))   #Make sure path is empty

 subprocess.call(['sex', './output/data_D'+str(NUM)+'.fits', '-c', './configfls/check.sex','-CATALOG_NAME',  str(NUM)+'.cat'])
 f = open(str(NUM)+'.cat')
 tl = open('./tfinds/list'+str(NUM)+'/transientlist.cat', 'w')
 tl.write('# Source Num'+'\n')
 tl.write('# RA'+'\n')
 tl.write('# DEC'+'\n')
 tl.write('# Ref flux'+'\n')
 tl.write('# Sci flux'+'\n')
 tl.write('\n')
 tl.write('\n')
 counter = 0 #counts number of detections

 dat1 = fits.getdata('./output/sci_cut'+str(NUM)+'.fits')
 dat1 = dat1.byteswap().newbyteorder()
 bkg = sep.Background(dat1)
 dat1 = dat1 - bkg

 dat2 = fits.getdata('./output/ref_cut'+str(NUM)+'.fits')
 dat2 = dat2.byteswap().newbyteorder()
 bkg = sep.Background(dat2)
 datB2 = dat2 - bkg #Need to check dat2 for 0

 dat3 = fits.getdata('./output/data_D'+str(NUM)+'.fits')
 dat3 = dat3.byteswap().newbyteorder()
 dat4 = fits.getdata('./output/data_Scorr'+str(NUM)+'.fits')
 dat4 = dat4.byteswap().newbyteorder()
 hed = fits.getheader('./output/data_Scorr'+str(NUM)+'.fits')
 World  = WCS(hed)

 if var == True: #finds all variables
  for line in f:
   fn = 0
   for i in line.split(' '):
    if i != '' and fn == 0:
     fn+= 1
     flux = float(i)
    elif i != '' and fn == 1 and abs(flux) > 450:# and 1106 < float(i) < 1108 :
     fn+=1
     x = i
    elif i != '' and fn == 2 and abs(flux) > 450: #and 436 < float(i) < 438 :
     fn+=1
     y = i
     f1 , ferr1, flag1 = sep.sum_circle(dat1, x ,y, 3.0)
     f2 , ferr2, flag2 = sep.sum_circle(dat2, x ,y, 3.0)
     fB2 , ferr2, flag2 = sep.sum_circle(dat2, x ,y, 3.0)
     f3 , ferr3, flag3 = sep.sum_circle(dat4, x ,y, 3.0)
     #print(max(f1, f2), f2, f3)

     if max(f1,fB2) > min(f1,fB2)*3  and f3 > 0.3 and f2 != 0:
      counter += 1
      print(f1, f2, f3)
      N = World.wcs_pix2world(float(x),float(y),1) #New coord
      tl.write(str(counter)+' '+str(N[0])+' '+str(N[1])+' '+str(f1)+'\n')
      os.makedirs('./tfinds/list'+str(NUM)+'/F%s' %(counter))

      CUT = cut(dat1, (N[0],N[1]), (50, 50), World)
      hdu = fits.PrimaryHDU(CUT.data, header=CUT.wcs.to_header())
      hdu.writeto('./tfinds/list'+str(NUM)+'/F%s/ref.fits' %(counter), overwrite=True)

      CUT = cut(datB2, (N[0],N[1]), (50, 50), World)
      hdu = fits.PrimaryHDU(CUT.data, header=CUT.wcs.to_header())
      hdu.writeto('./tfinds/list'+str(NUM)+'/F%s/sci.fits' %(counter), overwrite=True)

      CUT = cut(dat3, (N[0],N[1]), (50, 50), World)
      hdu = fits.PrimaryHDU(CUT.data, header=CUT.wcs.to_header())
      hdu.writeto('./tfinds/list'+str(NUM)+'/F%s/D.fits' %(counter), overwrite=True)

      CUT = cut(dat4, (N[0],N[1]), (50, 50), World)
      hdu = fits.PrimaryHDU(CUT.data, header=CUT.wcs.to_header())
      hdu.writeto('./tfinds/list'+str(NUM)+'/F%s/Scorr.fits' %(counter), overwrite=True)

 else: #find transients
  for line in f:
   fn = 0
   for i in line.split(' '):
    if i != '' and fn == 0:
     fn+= 1
     flux = float(i)
    elif i != '' and fn == 1 and flux > 450:# and 1106 < float(i) < 1108 :
     fn+=1
     x = float(i)

    elif i != '' and fn == 2 and flux > 450: #and 436 < float(i) < 438 :
     fn+=1
     y = float(i)
     f1 , ferr1, flag1 = sep.sum_circle(dat1, x ,y, 3.0)
     f2 , ferr2, flag2 = sep.sum_circle(dat2, x ,y, 3.0)
     fB2 , ferr2, flag2 = sep.sum_circle(datB2, x ,y, 3.0)
     f3 , ferr3, flag3 = sep.sum_circle(dat4, x ,y, 3.0)

     if f1 > fB2*2  and f3 > 0.3 and f2 != 0:
      counter += 1
      #print(x, y)
      N = World.wcs_pix2world(float(x),float(y),1) #New coord
      tl.write(str(counter)+' '+str(N[0])+' '+str(N[1])+' '+str(f1)+'\n')
      os.makedirs('./tfinds/list'+str(NUM)+'/F%s' %(counter))

      CUT = cut(dat1, (x,y), (50, 50), World)
      hdu = fits.PrimaryHDU(CUT.data, header=CUT.wcs.to_header())
      hdu.writeto('./tfinds/list'+str(NUM)+'/F%s/ref.fits' %(counter), overwrite=True)

      CUT = cut(datB2, (x,y), (50, 50), World)
      hdu = fits.PrimaryHDU(CUT.data, header=CUT.wcs.to_header())
      hdu.writeto('./tfinds/list'+str(NUM)+'/F%s/sci.fits' %(counter), overwrite=True)

      CUT = cut(dat3, (x,y), (50, 50), World)
      hdu = fits.PrimaryHDU(CUT.data, header=CUT.wcs.to_header())
      hdu.writeto('./tfinds/list'+str(NUM)+'/F%s/D.fits' %(counter), overwrite=True)

      CUT = cut(dat4, (x,y), (50, 50), World)
      hdu = fits.PrimaryHDU(CUT.data, header=CUT.wcs.to_header())
      hdu.writeto('./tfinds/list'+str(NUM)+'/F%s/Scorr.fits' %(counter), overwrite=True)
 subprocess.call(['rm', str(NUM)+'.cat'])
#####################################################################################################################


#####################################################################################################################
"""
Selects most appropriate ref tile
"""

def refsel(sci):
 print('finding ref tile')
 hdulist = fits.open(sci)
 header = hdulist[1]._header

 Tile_list = []


 F = glob.glob('./ref-tile/*') ###CHANGE HERE FOR REF TILE DIRECTORY###
 for i in F:
  hdulist2 = fits.open(i)
  header2 = hdulist2[1]._header
  #Group files that have the same filter object and UT as sci image
  if header['FILTER'] == header2['FILTER'] and header['OBJECT'] == header2['OBJECT'] and header['UT'] == header2['UT']:
   Tile_list.append(i)

 num = -1
 for i in Tile_list: #select most recent file
  hdulist2 = fits.open(i)
  header2 = hdulist2[1]._header
  num+=1
  if num == 0:
   date = float(header2['JD'])
   pos = num #file position in array
  else:
   if date < float(header2['JD']):
    date = float(header2['JD'])
    pos = num #file position in array

 return(Tile_list[pos])
#####################################################################################################################


#####################################################################################################################
""" 
These two classes allow the daemons to have children... how nice
Basically, a function in a pool can use its own pool
"""

class NoDaemonProcess(multiprocessing.Process):
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
""" 
compressed photo images are not compatible with SExtractor and SWARP
This function remakes the fits files so they are
"""
def rdfits(image):	
    hdulist =  fits.open(image)

    header= hdulist[1]._header

    data = hdulist[1].data
    image_RD = image.replace('.fits', '_RD.fits')

    fits.writeto(image_RD, data, header, overwrite=True)
    return(image_RD)

#####################################################################################################################

#####################################################################################################################
"""
Normalises the psf required for ZOGY, will remove any values smaller than the clean factor 
to avoid clutter by precision
"""

def clean_norm_psf (psf_ar, clean_fact = 0.25):
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

""" 
Finds the psf using PSFex (CITE). Uses the polynomial and expands it into a kernal 
the size of the subtraction images.
"""
 
def get_psf(image, xl, yl, const):
    sexcat = image.replace('.fits', '_PSFCAT.fits')
    sexcat = sexcat.replace('./output/','') 
    talk = ['sex' ,image,'-c','./configfls/default.sex' , '-CATALOG_NAME' , sexcat]
    print(sexcat)
    subprocess.call(talk)

    outcat = sexcat.replace('_PSFCAT.fits', '.psfexcat')

    talk2 = ['psfex', sexcat, '-c', './configfls/psfex.conf', '-OUTCAT_NAME', outcat]

    subprocess.call(talk2)

    with fits.open(sexcat.replace('.fits','.psf')) as hdulist:
     header = hdulist[1].header
     data = hdulist[1].data
    polzero1 = header['POLZERO1']
    polzero2 = header['POLZERO2']
    polscal1 = header['POLSCAL1']
    polscal2 = header['POLSCAL2']
    poldeg = header['POLDEG1']
    psf_samp = header['PSF_SAMP']

###!!!! INSERT IF STATEMENT FOR LOWER POLYNOMIAL !!!!####

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

    x = (xcenter_fft - polzero1) / polscal1
    y = (ycenter_fft - polzero2) / polscal2
    

    dat = data[0][0][:]
    if poldeg == 2:
     psf = dat[0] + dat[1] * x + dat[2] * x**2 + dat[3] * y + dat[4] * x * y + dat[5] * y**2
    elif poldeg == 3:
     psf = dat[0] + dat[1] * x + dat[2] * x**2 + dat[3] * x**3 + \
    dat[4] * y + dat[5] * x * y + dat[6] * x**2 * y + \
    dat[7] * y**2 + dat[8] * x * y**2 + dat[9] * y**3

    psf_ima_resized = ndimage.zoom(psf, psf_samp_update)
    psf_ima_resized_norm = clean_norm_psf(psf_ima_resized, const)
    psf_hsize = math.floor(psf_size/2)

    
    index = [slice(int(ycenter_fft-psf_hsize), int(ycenter_fft+psf_hsize+1)), 
     slice(int(xcenter_fft-psf_hsize), int(xcenter_fft+psf_hsize+1))]

    psf_ima_center[index] = psf_ima_resized_norm

    # perform fft shift
    psf_ima_shift = fft.fftshift(psf_ima_center)
    return(psf_ima_shift, sexcat, outcat)

#####################################################################################################################

#####################################################################################################################
"""
Using SWARP (cite) remap the reference image so it alligns with the co-ords of the science image.

NOTE: This is quite a cumbersome way of remapping, it's woth looking into a cleaner/faster way of doing it 
"""

def register(sci_im, ref_im):

    im_out = ref_im.replace('.fits','_REMAP.fits')

    hdulist_sci = fits.open(sci_im)
    header_sci = hdulist_sci[0].header

    hdulist_ref = fits.open(ref_im)
    header_ref = hdulist_ref[0].header

    head_out = header_sci[:]

    sci_xl = header_sci['ZNAXIS1']
    sci_yl = header_sci['ZNAXIS2']

    for i in ['NAXIS1', 'NAXIS2']:
     head_out[i] = header_sci['Z'+i]

    with open(im_out.replace('.fits','.head'),'w') as new_head:
     for card in head_out.cards:
      new_head.write(str(card)+'\n')


    TALK = ['swarp', ref_im, '-IMAGEOUT_NAME', im_out,  '-RESAMPLING_TYPE', 'LANCZOS3']

    subprocess.call(TALK)
    return(im_out)

#####################################################################################################################

####################################################################################################################
"""
Function that takes in output catalogs of stars used in the PSFex runs on the new and the ref image, and 
returns the arrays with pixel coordinates (!) x, y (in the new frame) and fratios for the matching stars. 
In addition, it provides the difference in stars' RAs and DECs in arcseconds between the two catalogs.

"""

def get_fratio(psfcat_sci, psfcat_ref, sexcat_sci, sexcat_ref):
    
    def readcat (psfcat):
        table = ascii.read(psfcat, format='sextractor')
        number = table['SOURCE_NUMBER']
        x = table['X_IMAGE']
        y = table['Y_IMAGE']
        norm = table['NORM_PSF']
        return number, x, y, norm
        
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
        # record in ra, dec for each source
        ra = []
        dec = []
        for n in number:
            ra.append(ra_sex[n-1])
            dec.append(dec_sex[n-1])      
        return np.array(ra), np.array(dec)
    
    # get ra, dec in pixel coords
    ra_sci, dec_sci = xy2radec(number_sci, sexcat_sci)
    ra_ref, dec_ref = xy2radec(number_ref, sexcat_ref)

    # now find matching sources, this step needs improving
    x_sci_match = []
    y_ref_match = []
    dra_match = []
    ddec_match = []
    fratio = []
    nmatch = 0
    for i_sci in range(len(x_sci)):
        # calculate distance to ref objects
        dra = 3600.*(ra_sci[i_sci]-ra_ref)*np.cos(dec_sci[i_sci]*np.pi/180.)
        ddec = 3600.*(dec_sci[i_sci]-dec_ref)
        dist = np.sqrt(dra**2 + ddec**2)
        # minimum distance and its index
        dist_min, i_ref = np.amin(dist), np.argmin(dist)
        #print(dist_min, 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
        if dist_min <250.: #This min distance is dependant on your registrtion. The less confident you are in your registration the bigger it needs to be.
            nmatch += 1
            x_sci_match.append(x_sci[i_sci])
            y_ref_match.append(y_sci[i_sci])
            dra_match.append(dra[i_ref])
            ddec_match.append(ddec[i_ref])
            # append ratio of normalized counts to fratios
            fratio.append(norm_sci[i_sci] / norm_ref[i_ref])

    return np.array(x_sci_match), np.array(y_ref_match), np.array(fratio), \
        np.array(dra_match), np.array(ddec_match)

####################################################################################################################

####################################################################################################################
"""
Reads in science image and cuts it into sub images.
Reads in ref image and remaps. 
Finally cuts the ref image in the same places
"""

def imprep(sci_im, ref_im, inject = False):
    F = glob.glob('./output/*')
    for fil in F:
     subprocess.call(['rm', fil])
    # This part formats the files so the other programs can use them#
    name1 = rdfits(sci_im)
    #name2 = rdfits(sci_im)
    #nameT = register(name2, name1)
    hdulist = fits.open(name1)
    data = hdulist[0].data
    header= hdulist[0].header

    W = WCS(header) #World coord system

    ##############################################
    # remap and stuff  #
    name2 = rdfits(ref_im)
    #name2 = rdfits(sci_im)
    #name3 = register(name2, name1)
    #subprocess.call(['rm', name1])
    #subprocess.call(['rm', name2])
    try:
     shutil.rmtree('alipy_out')
    except:
     print('no folder')
    subprocess.call(['python', 'ali.py' , name2, name1]) #alipy remap
    subprocess.call(['rm', name1, name2])

    name_new = glob.glob('./alipy_out/*')[0]

    hdulist2 = fits.open(name_new)
    data2 = hdulist2[0].data

    header2= hdulist[0].header
    W2 = WCS(header2)
    
    ############################################## 

    X = int(header['NAXIS1'])
    xcuts = math.ceil(X/2000) # number of cuts in the x plane
    Y = int(header['NAXIS2'])
    ycuts= math.ceil(Y/4000) #number of cuts in the y plane

    #print(XCORD)

    CC = [] #Centre Co-ordinates of each cut

    Xlin = np.linspace(0, X, xcuts)
    XC = [] #Centre of X points

    for i in range(len(Xlin)-1):
     XC.append( int((Xlin[i] + Xlin[i+1])/2) )

    Ylin = np.linspace(0, Y, ycuts)
    YC = [] #Centre of Y points

    for i in range(len(Ylin)-1):
     YC.append( int((Ylin[i] + Ylin[i+1])/2) )

    for i in range (len(XC)):
     for j in range(len(YC)):
      CC.append([XC[i], YC[j]])

    for i in range(len(CC)):
     CUT = cut(data, (CC[i][0],CC[i][1]), (2000, 4000), W) #Create a subimage with centre co-ords CC
     if inject == True:
      datainj = CUT.data
      for inj in range(randint(0,3)):
       Xinj = randint(1, CUT.data.shape[0])
       Yinj = randint(1, CUT.data.shape[1])
       datainj[Xinj][Yinj] = 10000
       datainj[Xinj][Yinj+1] = 10000
       datainj[Xinj][Yinj-1] = 10000
       datainj[Xinj+1][Yinj] = 10000
       datainj[Xinj-1][Yinj] = 10000
       datainj[Xinj+1][Yinj-1] = 10000
       datainj[Xinj+1][Yinj+1] = 10000
       datainj[Xinj-1][Yinj-1] = 10000
       datainj[Xinj-1][Yinj+1] = 10000
       datainj[Xinj+2][Yinj-2] = 10000
       datainj[Xinj+2][Yinj-1] = 10000
       datainj[Xinj+2][Yinj] = 10000
       datainj[Xinj+2][Yinj+1] = 10000
       datainj[Xinj+2][Yinj+2] = 10000
       datainj[Xinj-2][Yinj-1] = 10000
       datainj[Xinj-2][Yinj-2] = 10000
       datainj[Xinj-2][Yinj] = 10000
       datainj[Xinj-2][Yinj+1] = 10000
       datainj[Xinj-2][Yinj+2] = 10000
       datainj[Xinj][Yinj+2] = 10000
       datainj[Xinj][Yinj-2] = 10000
       datainj[Xinj-1][Yinj+2] = 10000
       datainj[Xinj+1][Yinj-2] = 10000
       datainj[Xinj+1][Yinj+2] = 10000
       datainj[Xinj-1][Yinj-2] = 10000



       print(i, inj, Xinj, Yinj)


     OC = W.wcs_pix2world(CC[i][0],CC[i][1],1) #original coords
     NC = W2.wcs_world2pix(OC[0],OC[1],1) #New coords

     CUT2 = cut(data2, (NC[0], NC[1]), (2000, 4000), W2)

     hdu = fits.PrimaryHDU(CUT2.data, header=CUT2.wcs.to_header())
     hdu.writeto('output/ref_cut%s.fits' %(i+1), overwrite=True)
     if inject == True:
      hdu = fits.PrimaryHDU(datainj, header=CUT.wcs.to_header())
      hdu.writeto('output/sci_cut%s.fits' %(i+1), overwrite=True)
     else:
      hdu = fits.PrimaryHDU(CUT.data, header=CUT.wcs.to_header())
      hdu.writeto('output/sci_cut%s.fits' %(i+1), overwrite=True)

#####################################################################################################################

#####################################################################################################################
"""
Optimal image subtraction in a pythonic layout! 
Where the magic happens
"""
def ZOGY(R,N,Pr,Pn,sr,sn,fr,fn,Vr,Vn,dx,dy):

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

    denominator = sn2*fr2*Pr_hat2_abs + sr2*fn2*Pn_hat2_abs
    if np.any(denominator==0):
        print('There are zeros!')

    D_hat = (fr*Pr_hat*N_hat - fn*Pn_hat*R_hat) / np.sqrt(denominator)

    D = np.real(fft.ifft2(D_hat)) / fD
    
    P_D_hat = (fr*fn/fD) * (Pr_hat*Pn_hat) / np.sqrt(denominator)
    
    S_hat = fD*D_hat*np.conj(P_D_hat)
    S = np.real(fft.ifft2(S_hat))

    kr_hat = fr*fn2*np.conj(Pr_hat)*Pn_hat2_abs / denominator
    kr = np.real(fft.ifft2(kr_hat))
    kr2 = kr**2
    kr2_hat = fft.fft2(kr2)

    kn_hat = fn*fr2*np.conj(Pn_hat)*Pr_hat2_abs / denominator
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

    F_S =  np.sum((fn2*Pn_hat2_abs*fr2*Pr_hat2_abs) / denominator)
    F_S /= R.size

    alpha = S / F_S
    alpha_std = np.zeros(alpha.shape)
    alpha_std[V_S>=0] = np.sqrt(V_S[V_S>=0]) / F_S

    return D, S, S_corr, alpha, alpha_std
####################################################################################################################

####################################################################################################################
""" 
Using all of the above, this function will find the psf of background subtracted data.
After will find the F-ratio, and pixel properties. The subtraction occurs producing D, S, and Scorr images.
"""

def finp(image):
 Data_Df = []
 Data_Sf = []
 Data_Cf = []
 sci = image.replace('ref', 'sci')

 SH = fits.getdata(sci).shape
 #Parallell step
 with multiprocessing.Pool(2) as p:
  V1, V2 = p.starmap(get_psf, [(sci, SH[1], SH[0], 0.95), (image, SH[1], SH[0], 0.95)])
 psf, sexcat1, psfcat1 = V1#get_psf(sci, SH[1], SH[0],0.95)
 psf2, sexcat2, psfcat2 = V2#get_psf(image, SH[1], SH[0],0.95)
 x_fratio, y_fratio, fratio, dra, ddec = get_fratio(psfcat1, psfcat2, sexcat1, sexcat2)
 FM = np.median(fratio)

 fnum = image.replace('./output/ref_cut','')
 fnum2= fnum.replace('.fits','')

 ######## Science data ########

 hdu = fits.open(sci)
 dat = hdu[0].data
 dat = dat.byteswap().newbyteorder()
 head = hdu[0].header
 W = WCS(head)
 bkg = sep.Background(dat)
 std = bkg.rms()

 sub_dat = dat - bkg #bkg subtracted data

 ######### Ref data ########
 hdu2 = fits.open(image)
 dat2 = hdu2[0].data
 dat2 = dat2.byteswap().newbyteorder()
 head2 = hdu2[0].header
 W2 = WCS(head)
 bkg2 = sep.Background(dat2)
 std2 = bkg2.rms()

 sub_dat2 = dat2 - bkg2 #bkg subtracted data

 if FM > 2.5 or FM < 0.09 or math.isnan(FM)==True:
  FM = 1.0 # Fratio can be unstable. If it's to high or too low subtraction breaks
  print('FM CHANGED')

 f_new = 1.0
 f_ref = f_new/FM

 dx = dra / 1.24
 dy = ddec / 1.24

 dx_full = np.sqrt(np.median(dx)**2 + np.std(dx)**2)
 dy_full = np.sqrt(np.median(dy)**2 + np.std(dy)**2) 

 var_sci = sub_dat + 1.2**2 
 var_ref = sub_dat2 + 1.2**2 

 data_D, data_S, data_Scorr, fpsf, fpsf_std  = ZOGY(sub_dat2, sub_dat,  psf2, psf, np.median(std2), np.median(std), f_ref, f_new, var_ref, var_sci, dx_full, dy_full)

 hdu = fits.PrimaryHDU(data_D, header= head)
 hdu.writeto('./output/data_D'+fnum,overwrite=True)

 hdu = fits.PrimaryHDU(data_S, header= head)
 hdu.writeto('./output/data_S'+fnum,overwrite=True)

 hdu = fits.PrimaryHDU(data_Scorr, header= head)
 hdu.writeto('./output/data_Scorr'+fnum,overwrite=True)

 Data_Df.append('./output/data_D'+fnum)
 Data_Sf.append('./output/data_S'+fnum)
 Data_Cf.append('./output/data_Scorr'+fnum)


 subprocess.call(['rm', 'sci_cut%s_PSFCAT.psf' %(fnum2), 'sci_cut%s.psfexcat' %(fnum2), 'sci_cut%s_PSFCAT.fits' %(fnum2)])
 subprocess.call(['rm', 'ref_cut%s_PSFCAT.psf' %(fnum2), 'ref_cut%s.psfexcat' %(fnum2), 'ref_cut%s_PSFCAT.fits' %(fnum2)])
####################################################################################################################


"""The Program"""

t0 = time.time()
x = multiprocessing.cpu_count()
print('you have', x, 'cores')

if len(sys.argv) == 1:
 print(' ')
 print(' ')
 print('This is ZiP, the image subtraction software! This particular version has been focused towards GOTO data')
 print('--------------------------------------------------------------------------------------------------------------')
 print(' ')
 print(' ')
 print('To use this simply type [python3 ZOGYP.py sci_im ref_im]')
 print('where images are fits files you want subtracting')
 print('or if you want to see if the software is working [python3 ZOGYP.py test]')
 print(' ')
 print(' ')
 print('If you have directory with a selection of ref tiles, just submit the sci image')
 print('and the ref selctor can find the most fitting ref tile')

elif sys.argv[1] == 'test':
 if x < 6:
  print('Serial version')
  ncores = x
  fin('./test/2.fits', './test/1.fits')
 else:
  if x>45:
   ncores = 15
  else:
   ncores = (int(x/3))
  print('Parallell version, using %s cores' %(ncores*3))
  imprep('./test/2.fits', './test/1.fits')
  refs = glob.glob('./output/ref_cut*.fits')
  p = NoDaemonProcessPool(processes = ncores)
  p.starmap(finp, product(refs, repeat=1))

else:
 if x < 6:
  print('Serial version')
  ncores = x
  if len(sys.argv) == 3:
   fin(sys.argv[1], sys.argv[2])
  else:
   ref = refsel(sys.argv[1])
   fin(sys.argv[1], ref)
 else:
  if x>45:
   ncores = 15
  else:
   ncores = (int(x/3))
  print('Parallell version, using %s cores' %(ncores*3))
  if len(sys.argv) == 3:
   imprep(sys.argv[1], sys.argv[2])
   ref = sys.argv[2] #removing stuff
   refs = glob.glob('./output/ref_cut*.fits')
   p = NoDaemonProcessPool(processes = ncores)
   p.starmap(finp, product(refs, repeat=1))
  else:
   ref = refsel(sys.argv[1]) #Select ref image
   imprep(sys.argv[1], ref)
   refs = glob.glob('./output/ref_cut*.fits')
   p = NoDaemonProcessPool(processes = ncores)
   p.starmap(finp, product(refs, repeat=1))
 #subprocess.call(['rm', ref.replace('.fits','_RD_REMAP.fits'),ref.replace('.fits','_RD_REMAP.head')])
 #subprocess.call(['rm', sys.argv[1].replace('.fits','_RD_REMAP.fits'),sys.argv[1].replace('.fits','_RD_REMAP.head')])

NUMB = len(glob.glob('./output/*Scorr*'))
if NUMB > ncores:
 p = Pool(ncores)
 p.map(fitra, range(1, NUMB+1))
else:
 p = Pool(NUMB)
 p.map(fitra, range(1, NUMB+1))

print(glob.glob('./tfinds/*/F*'))
t1 = time.time()
print((t1 -t0)/60 , 'minutes')
