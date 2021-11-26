import sys
import os
import math

import numpy as np
import scipy.constants as constants
import Hirogen_Functions
from astropy.io import fits


def flux_scaler(flux, wavelengths, line_names, scale_value):
    scale_array = np.ones(len(flux))
    lines = Hirogen_Functions.lines_for_analysis()
    for i, name in enumerate(line_names):
        line = lines[name]
        line_location = line[0]
        for i, x in enumerate(scale_array):
            if wavelengths[i].value >= line_location - 12 and wavelengths[i].value <= line_location + 12 and flux[i] >= 0: 
                scale_array[i] = scale_value
            elif wavelengths[i].value >= line_location - 12 and wavelengths[i].value <= line_location + 12 and flux[i] < 0: 
                scale_array[i] = 1/scale_value   
    scaled_flux = flux*scale_array
    return scaled_flux

def fits_file_gen(wave, flux, flux_err):
    hdu = fits.BinTableHDU.from_columns(
        [fits.Column(name='flux', format='E', array = flux),
        fits.Column(name='loglam', format='E', array = wave),
        fits.Column(name = 'ivar', format = 'E', array = flux_err)]
    )

    hdu.writeto('scaled_spec-0907-52373-0419.fits', overwrite = True)
    
    return None