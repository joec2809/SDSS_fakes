import sys
import os
import math

import numpy as np
import scipy.constants as constants
import Hirogen_Functions

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

def continua_maker(spec_region, flux, line_name, line_loc, object_name):
    """This version of the continua generator and corrector uses a wavelength region rather than a relative velocity"""

    if line_name.lower() in {'halpha', '[nii]6548', '[nii]6584'}:

        # These are the values Decker used for her paper with Or.
        # They will include parts of broad features so I am tweaking.
        """
        blue_start = 6492.8
        blue_end = 6532.8

        red_start = 6592.8
        red_end = 6632.8
        """
        blue_start = 6465
        blue_end = 6470

        red_start = 6620
        red_end = 6635

    elif line_name.lower() in {'hbeta', '[oiii]5007'}:

        blue_start = line_loc - 35
        blue_end = line_loc - 30

        red_start = line_loc + 30
        red_end = line_loc + 35

    elif line_name.lower() in {'[fexiv]5304', 'heii4686'}:

        blue_start = line_loc - 30
        blue_end = line_loc - 25

        red_start = line_loc + 25
        red_end = line_loc + 30

    elif line_name.lower() in {'[oiii]4959'}:

        blue_start = line_loc - 20
        blue_end = line_loc - 15

        red_start = line_loc + 15
        red_end = line_loc + 20

    else:
        blue_start = line_loc - 50
        blue_end = line_loc - 30

        red_start = line_loc + 30
        red_end = line_loc + 50

    continuum_regions = [blue_start, blue_end, red_start, red_end]

    blue_middle = (blue_start + blue_end) / 2
    red_middle = (red_start + red_end) / 2

    blue_flux = []
    red_flux = []

    for xx, wavelength_point in enumerate(spec_region):
        if blue_start < spec_region[xx] < blue_end:
            blue_flux.append(flux[xx])
        if red_start < spec_region[xx] < red_end:
            red_flux.append(flux[xx])

    try:
        m = (np.nanmean(red_flux) - np.nanmean(blue_flux)) / (red_middle - blue_middle)
        c = np.nanmean(blue_flux) - (m * blue_middle)

    except FloatingPointError:
        print("Something has gone wrong with the selection of part of the continuum - likely no valid points"
              "\nSkipping object for now")
        print(f"{object_name}\t{line_name}")

        print(spec_region)

        return [], []

    continuum = []

    for xx, velocity in enumerate(spec_region):
        continuum.append((m * spec_region[xx]) + c)

    try:
        scaled_flux = np.array(flux) / np.array(continuum)
    except:
        print("Continuum Scaling Failure")
        print(flux)
        print(continuum)

        scaled_flux = np.array(flux)

    return continuum, scaled_flux, continuum_regions