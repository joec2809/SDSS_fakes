import sys
import os
import math

import numpy as np
import scipy.constants as constants
import Hirogen_Functions

def peak_scaler(flux, wavelengths, line_names, scale_value):
    scale_array = np.ones(len(flux))
    lines = Hirogen_Functions.lines_for_analysis()
    for i, name in enumerate(line_names):
        line = lines[name]
        line_location = line[0]
        for i, x in enumerate(scale_array):
            if wavelengths[i].value >= line_location - 12 and wavelengths[i].value <= line_location + 12: 
                scale_array[i] = scale_value
    scaled_peaks = flux*scale_array
    return scaled_peaks