# ZOGYP.run_ZOGY

`ZOGYP.run_ZOGY`(sci_im, ref_im,  xslice=1, yslice=1, align = False, Ex = 'N', figs = False, sub_imagex = 1, sub_imagey =1)
  
  Using the ZOGY algorithm, subtracts ref_im from sci_im
  
  ---
  
  ## Parameters
  
   * **sci_im, ref_im** - Input fits files
  
   * **xslice, yslice** - Will split the image into (xslice by yslice) and apply a varrying PSF to each segment
  
   * **align** - True/False. If true, Alipy will be used to register the sci_im to the ref_im. If False the images are assumed to be alligned already.
  
   * **Ex** - N/V/T. Extractor run. If N, no extrator run. If V, looks for all variables. If T, looks just for transients (sources that get brighter)
  
   * **figs** - Create thumbnails of extracted sources
  
   * **sub-imagex, sub-imagey** - Only use for big images. Chops the image into sub images. Used for images that use up too much memory or with a complex PSF. (This can also be used to speed up the parallel version, however can result in lower quality PSF estimation.)
  
  ## Output
  
   * **D-image** - Raw subtraction
   
   * **S-image** - PSF Convolution of the raw image
   
   * **Scorr** - A statistical map in image space. Where bright sources are more probable residuals.
   
  
  
  
