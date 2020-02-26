Spali2
======

A robust image registration technique for wide fields of view 
--------------------------------------------------------------

*spali2 is the successor to the detection based image registration algorithm*  `spalipy <https://github.com/Lyalpha/spalipy>`_ *. Here, this section will cover the functionality and how to use spali2*

The Algorithm
-------------

Firstly, two catalogs are built using *sextractor*. These catalogs will detail source positions and fluxes. Quads are built in both catalogs ( `Hogg et al 2010 <https://iopscience.iop.org/article/10.1088/0004-6256/139/5/1782/pdf>`_ ) and are then compared to one another. Using the position, rotation, and scaling differences, an affine transform is built (`scipy <https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.ndimage.interpolation.affine_transform.html>`_). **!!! WARNING: the affine transformation uses scipy.ndimage.interpolation.affine_transform which does not handle NaNs properly. Replace all NaN values in the input image prior to running spalipy !!!**

Second, once the affine transform is done, the remaining residuals are compared to see where higher order transfroms are needed. This builds a 2D-spline surface to represent the residuals in x and y axes. The spline surfaces are then implemeneted to properly register the images.

Using Spali2
------------

As usual, import stuff
::
   
   import glob
   import ntpath
   import time
   import shutil
   import subprocess
   import matplotlib.pyplot as plt
   from astropy.io import fits
   import numpy as np

   #ZOGY in Parallel routines
   from zogyp.zip import run_ZOGY
   from zogyp.zip import rdfits
   from zogyp.zip import config_loc

   #Image alignment
   from zogyp.spali2 import spalipy

   #Stacking
   #from zogyp.zo_coadd import med_comb
   #from zogyp.zo_coadd import prop_coad


To align image1 to image2 
::

   Output = spalipy(image1, image2, name='SPLINE.fits')

**If your images are not sextractor compliant, you can remake the fits files and try again**
::
   
   R = rdfits(image1) 
   D = rdfits(image2)
   Output = spalipy(R, D, name='SPLINE.fits')
   subprocess.call(['rm', R, D])
   
Finally, if your images are small and not prone to wide distorions or you are really pressed for time, you can take out the spline step.
::

   ali = spalipy(image1, image2, spline=False)

