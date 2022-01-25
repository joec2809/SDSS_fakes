import matplotlib.pyplot as plt

import Hirogen_Functions
import mods_functions

import sys
import os
import math

import numpy as np
import scipy.constants as constants
import mysql.connector
import decimal

from astropy.io import fits
from astropy import units as u
from mysql.connector import errorcode
from astropy.convolution import convolve, Box1DKernel




filename_1 = "./dr16/sdss/spectro/redux/26/spectra/0423/spec-0423-51821-0013.fits"
filename_2 = "./dr16/sdss/spectro/redux/26/spectra/0423/spec-0423-51821-0013-fake.fits"

spectrum = Hirogen_Functions.spec_reader(filepath=filename_2)
spec_header = spectrum[0]

data = spectrum[1]

lamb_observed = data['loglam'] * u.AA  # Observed wavelength

flux_in = data['flux'] * u.Unit('erg cm-2 s-1 AA-1')  # Flux
error = data['ivar']  # Flux error

flags = data['and_mask']  # Bitmask and mask map


print(flux_in)
print(len(lamb_observed), len(flux_in), len(error), len(flags))