import sys
import os
import math

import numpy as np
import scipy.constants as constants
import Hirogen_Functions
from astropy.io import fits
from astropy import units as u

c = constants.value('speed of light in vacuum') / 1000


def flux_scaler(flux, wavelengths, line_names, scale_value):
    object_name = ''
    flux_without_peaks = np.zeros(len(flux))
    scale_array = np.ones(len(flux))
    lines = Hirogen_Functions.lines_for_analysis()
    for ii, item in enumerate(lines):
        if item in line_names:
            line_location = lines[item][0]
            shift = (((np.array(wavelengths) * c) / line_location) - c)
            shift_region, wave_region, flux_region, error_region = Hirogen_Functions.region_cutter(
                    shift, 
                    wavelengths,
                    flux*u.Unit('erg cm-2 s-1 AA-1'),
                    line_location - 40,
                    line_location + 40,
                    mode='Wavelength'
                )
            continuum = Hirogen_Functions.continua_maker(
                    wave_region,
                    flux_region,
                    # shift=Shift_Region,
                    item,
                    line_location,
                    object_name
                )
            peak_start = find_nearest(wavelengths, continuum[2][0])
            flux_without_peaks[peak_start:peak_start+len(continuum[0])] = continuum[0]

    flux = flux_without_peaks

    for i, name in enumerate(line_names):
        line = lines[name]
        line_location = line[0]
        for i, x in enumerate(scale_array):
            if wavelengths[i].value >= line_location - 12 and wavelengths[i].value <= line_location + 12 and flux[i] >= 0: 
                scale_array[i] = scale_value
            elif wavelengths[i].value >= line_location - 12 and wavelengths[i].value <= line_location + 12 and flux[i] < 0: 
                scale_array[i] = 1/scale_value
    scaled_flux = flux*scale_array
    scaled_flux += flux_without_peaks
    flux += flux_without_peaks
    return scaled_flux

def fits_file_gen(wave, flux, flux_err, filename):
    hdu = fits.BinTableHDU.from_columns(
        [fits.Column(name='flux', format='E', array = flux),
        fits.Column(name='loglam', format='E', array = wave),
        fits.Column(name = 'ivar', format = 'E', array = flux_err)]
    )

    hdu.writeto(f"scaled_{filename}", overwrite = True)
    
    return None

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx