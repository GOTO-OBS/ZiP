#!/usr/bin/env python

"""
This script is to make the program executable in the shell
"""

import sys
import glob
import ntpath
import time
import shutil
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np

#ZOGY in Parallel routines
from zogyp.zip import run_ZOGY

def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1", "True")

def main():
    skip = 0 #lazy way of making errors easy to flag
    cleanref = 0.25
    cleansci = 0.25
    xsl = 1
    ysl = 1
    psfy = 1
    psfx = 1
    BO = False
    al = False
    ONM = 'ZOGY_sub'

    if len(sys.argv) == 1:
        print('~~~~~  ZOGY in Parallel ~~~~~')
        print('~           V1.6.3          ~')
        print('-----------------------------')
        print(' ')
        print('          Ry Cutter          ')
        print(' ')
        print(' ')
        print('Type: "ZOGY -help" for list of inputs')
        print(' ')

    elif sys.argv[1] == '-help':
        print('~~~~~  ZOGY in Parallel ~~~~~')
        print('~           V1.6.3          ~')
        print('-----------------------------')
        print(' ')
        print('          Ry Cutter          ')
        print(' ')
        print(' ')
        print('To use, type:')
        print(' ZOGY sci.fits ref.fits')
        print(' ')
        print(' Other Parameters: \n \n -sub_imagex [int] #number of x sub image slices \n -sub_imagey [int] #number of y sub image slices')
        print(' -xslice [int] #number of x PSF slices \n -yslice [int] #number of y PSF slices \n -blackout [True/False] #remove data where images do not overlap')
        print(' -clean_sci [0-1] #handles overfitting of PSF model for sci.fits \n -clean_ref [0-1] #handles overfitting of PSF model for ref.fits')
        print(' align [True/False] #aligns sci to ref image befroe subtrcation \n -outname [str] #filenames for output')
        print('------- \n \n')
        print('example:\n')
        print('ZOGY sci.fits ref.fits -sub_imagex 2 -sub_imagey 2 -clean_ref 0.29')
        

   
    else:
        for i in range(len(sys.argv)):
            if i == 0:
                continue
            elif i == 1:
                sci_T = sys.argv[1]
            elif i == 2:
                ref_T = sys.argv[2]
            elif skip == 1:
                skip = 0
                continue
            elif sys.argv[i] == '-sub_imagex':
                xsl = int(sys.argv[i+1])
                skip = 1
                if 'XSL_v' in locals():
                    raise ValueError('-sub_imagex assigned multiple times')
                else:
                    XSL_v = True

            elif sys.argv[i] == '-sub_imagey':
                ysl = int(sys.argv[i+1])
                skip = 1
                if 'YSL_v' in locals():
                    raise ValueError('-sub_imagey assigned multiple times')
                else:
                    YSL_v = True

            elif sys.argv[i] == '-xslice':
                psfx = int(sys.argv[i+1])
                skip = 1
                if 'PXSL_v' in locals():
                    raise ValueError('-xslice assigned multiple times')
                else:
                    PXSL_v = True

            elif sys.argv[i] == '-yslice':
                psfy = int(sys.argv[i+1])
                skip = 1
                if 'PYSL_v' in locals():
                    raise ValueError('-yslice assigned multiple times')
                else:
                    PYSL_v = True

            elif sys.argv[i] == '-clean_ref':
                cleanref = float(sys.argv[i+1])
                skip = 1
                if 'ACr_v' in locals():
                    raise ValueError('-clean_ref assigned multiple times')
                else:
                    ACr_v = True

            elif sys.argv[i] == '-clean_sci':
                cleansci = float(sys.argv[i+1])
                skip = 1
                if 'ACs_v' in locals():
                    raise ValueError('-clean_sci assigned multiple times')
                else:
                    ACs_v = True

            elif sys.argv[i] == '-outname':
                ONM = sys.argv[i+1]
                skip = 1
                if 'ONM_v' in locals():
                    raise ValueError('-outname assigned multiple times')
                else:
                    ONM_v = True

            elif sys.argv[i] == '-blackout':
                BO = str2bool(sys.argv[i+1])
                skip = 1
                if 'BO_v' in locals():
                    raise ValueError('-blackout assigned multiple times')
                else:
                    BO_v = True

            elif sys.argv[i] == '-align':
                al = str2bool(sys.argv[i+1])
                skip = 1
                if 'BOAL_v' in locals():
                    raise ValueError('-align assigned multiple times')
                else:
                    BOAL_v = True

            else:
                raise ValueError(sys.argv[i]+' is not a valid input parameter \n --- use "ZOGY -help" to get list of valid inputs') 


        run_ZOGY(sci_T, ref_T, clean_ref=cleanref, clean_sci=cleansci,
                 sub_imagex = xsl, sub_imagey = ysl, xslice = psfx, yslice = psfy,
                 blackout = BO, align = al, outname = ONM)








