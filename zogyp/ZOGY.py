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



run_ZOGY(sys.argv[1], sys.argv[2])
