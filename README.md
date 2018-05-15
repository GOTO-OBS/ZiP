# ZiP 

ZOGY in Parallell (ZiP) is a fast(ish) computation of proper image subtraction [B.Zackay, E.Ofek, A.Gal-Yam 2016](http://iopscience.iop.org/article/10.3847/0004-637X/830/1/27/pdf). Inspired by [Ofek 2014](http://adsabs.harvard.edu/abs/2014ascl.soft07005O) and [pmvreeswijk](https://github.com/pmvreeswijk/ZOGY). ZiP offers a faster subtraction at the expense of a more comprehensive input. I.e. The program should be tailored for one telescope or input of images. This code has a parallell function, however it requires 6+ cores to operate. This particular Case is the Gravitational-Wave Optical Transient Observer ([GOTO](https://goto-observatory.org/)) However, simple fudging of the parameters should make it possible to make this work for other telescopes.

In Serial the program takes ~ 1:45 per subtraction

In Parallell it takes ~ 38s per subtraction

## Other needed programs + modules:
* Python - 3.6.4 (+ 2.4 for alipy)
* [SExtractor](https://www.astromatic.net/software/sextractor) - 2.19.5
* ~~[Swarp](https://www.astromatic.net/software/swarp) - 2.38.0~~
* [PSFex](https://www.astromatic.net/software/psfex) - 3.17.1
* [pyfftw](https://hgomersall.github.io/pyFFTW/) - 0.10.3
* [SEP](http://sewpy.readthedocs.io/en/latest/) - 1.0.1
* [astropy](http://www.astropy.org/) - 2.0.4
* [numpy](http://www.numpy.org/) - 1.14.0
* [scipy](https://www.scipy.org/) - 1.0.0
* [alipy](https://obswww.unige.ch/~tewes/alipy/)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ryan Cutter 
V1.3.5 (15/05/2018)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
