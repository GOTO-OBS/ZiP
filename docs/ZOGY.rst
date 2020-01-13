ZOGY
====

The Maths
---------

TO DO

Running ZOGY Tutorial
---------------------

Install the standards ::
    
        import glob
        import ntpath
        import time
        import shutil 
        import matplotlib.pyplot as plt
        from astropy.io import fits
        import numpy as np

        #ZOGY in Parallel routines
        from zogyp.zip import run_ZOGY
        from zogyp.zip import rdfits
        from zogyp.zip import config_loc

        #Image alignment
        #from zogyp.spali2 import spalipy

        #Stacking
        #from zogyp.zo_coadd import med_comb
        #from zogyp.zo_coadd import prop_coad
      
**Find the directory of the configuration files, these will need to be edited for best peformance** ::

       print(config_loc())
The two files of interest are default.sex and psfex.conf

Get the test images :: 

    t_f = LOC.replace('configfls','test')
    T = [i for i in glob.glob(t_f+'/*')]
    
Basic Subtraction :: 

   run_ZOGY(T[0], T[1], outname='TEST_SUB') #ZOGY image subtraction T[0] - T[1]
   
Let's plot the subtractions ::

   file_names = glob.glob('Zoutput/*')
   fig, axs = plt.subplots(1, 5, figsize=(20,45))
   for i in range(len(file_names)):
       D = fits.getdata(file_names[i])
       axs[i].imshow(D, cmap = 'gray_r', vmin = -np.mean(D)*20, vmax= np.mean(D)*20)
       axs[i].set_xticks([] , [])
       axs[i].set_yticks([] , [])
       axs[i].set_title(ntpath.basename(file_names[i]))
   plt.show()
   
..:image:: 

Tips and Tricks
---------------

To Do

