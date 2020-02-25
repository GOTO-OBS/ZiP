# ZiP 

ZOGY in Parallell (ZiP) is a fast(ish) computation of proper image subtraction [B.Zackay, E.Ofek, A.Gal-Yam (2016)](http://iopscience.iop.org/article/10.3847/0004-637X/830/1/27/pdf). Inspired by [Ofek 2014](http://adsabs.harvard.edu/abs/2014ascl.soft07005O) and [pmvreeswijk](https://github.com/pmvreeswijk/ZOGY). ZiP offers a faster subtraction at the expense of a more comprehensive input. I.e. The program should be tailored for one telescope or input of images. This code has a parallell function, however it requires 6+ cores to operate. This particular Case is optimised for the Gravitational-Wave Optical Transient Observer ([GOTO](https://goto-observatory.org/)) However, simple fudging of the parameters should make it possible to make this work for other telescopes.

An internal version of [spalipy](https://github.com/GOTO-OBS/spalipy) has been added as the alignment algorithm. This uses sextractor to find source locations in two images and then aligns them with an affine transform. Residuals are then used to build a 2D spline surface to correct for warping due to large field distortions.

Finally, a parallel version of [proper coadition](https://arxiv.org/abs/1512.06879) is used for stacking images. It still isn't increadibly fast for on the spot coaddition; so a meidian combine tool is also included.

---

In Serial the program takes ~ 2:06 per subtraction [for a field 8000 X 6000 pixels big]

In Parallell it takes ~ 33s per subtraction [for a field 8000 X 6000 pixels big]


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ryan Cutter 
V1.4.00 (25/02/2019)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
