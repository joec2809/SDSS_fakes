import matplotlib.pyplot as plt
import mpl_style
import Hirogen_Functions
import mods_functions
plt.style.use(mpl_style.scientific)

import sys
import os
import math
import spectres

import numpy as np
import scipy.constants as constants
import mysql.connector

from astropy.io import fits
from astropy import units as u
from astropy.convolution import convolve, Box1DKernel
from scipy.signal import medfilt
from scipy.ndimage import gaussian_filter1d
from astropy.nddata import InverseVariance
from dust_extinction.parameter_averages import F99
from numpy.fft import *
from mysql.connector import errorcode



np.set_printoptions(threshold=sys.maxsize)
c = constants.value('speed of light in vacuum') / 1000

def find_nearest_idx(arr, val):
    idx = (np.abs(arr - val)).argmin()
    return idx

smoothing_boxcar = 5


User = 'Joe'
User_Config = Hirogen_Functions.user_config(user=User)

#TableID = "SDSS_Test_QSO_Sample"
#TableID = "SDSS_Test_Sample"
TableID = "SDSS_Confirmed_Objects"

config_parameters = Hirogen_Functions.main_config()  # Draws from centralised parameter declarations

Lower_Wave = config_parameters[0]
Upper_Wave = config_parameters[1]

Lower_Shift = config_parameters[2]
Upper_Shift = config_parameters[3]

# Any object with a candidate score equal to or exceeding this value will be treated as an ECLE candidate
Candidate_Threshold = config_parameters[4]

# Maximum possible candidate score
Candidate_Score_Max = config_parameters[5]

# To qualify as a line detection for candidate selection the overall feature pEQW must exceed this value:
LineDetection_pEQW_Threshold = config_parameters[6]
LineDetection_Max_Threshold = config_parameters[7]

LineDetection_Max_Above_Av_Continua_Threshold = config_parameters[8]

# To qualify as a line detection for candidate selection the maximum flux point of the feature must occur
# within this +_ kms^-1 of the zero velocity point
LineDetection_Peak_Tolerance = config_parameters[9]
Line_Peak_Location_Region = config_parameters[10]

Line_Peak_Region_Minima_Threshold = config_parameters[11]

Strong_EQW_Threshold = config_parameters[12]
Strong_Peak_Max_Threshold = config_parameters[13]

# Search for 'QUERY' to find the actual database access location where parameters can be adjusted

Database = User_Config[0]
Database_User = User_Config[1]
Database_Password = User_Config[2]
Main_Spectra_Path = User_Config[3]

Data = Hirogen_Functions.database_connection(user=Database_User, password=Database_Password, database=Database)

print("\n")
cursor = Data.cursor()



# Can be set to only run the analysis for specific IDs given in the following format
IDs = 1021306429161629696

# Manually_Inspected_Flag = -10 are spectra with external issues that prevent their use

# QUERY
cursor.execute(
    f"SELECT DR16_Spectroscopic_ID, spec_Plate, spec_MJD, spec_FiberID, z_SDSS_spec, generalised_extinction, "
    f"survey, run2d, Standard_Inclusion,Path_Override_Flag, Path_Override, Follow_Up_ID, Smoothing_Override,"
    f"z_corrected_flag, extinction_corrected_flag "
    f"FROM `{Database}`.`{TableID}`"
    #f"WHERE Manually_Inspected_Flag != -10 AND lin_con_LineFlux_Halpha is NULL"
    #f"WHERE lin_con_pEQW_Halpha is NULL"
    #f"WHERE z_SDSS_spec >= 0.2"
    #f"WHERE Follow_Up_ID = 0"
    # f"WHERE Follow_Up_ID in (1,2,3)"
    #f"WHERE Standard_Inclusion = 1"
    f"WHERE DR16_Spectroscopic_ID = {IDs}"
)

Candidate_Data = cursor.fetchall()
# Candidate_Data = cursor.fetchmany(200)

#print(Candidate_Data)



Object_Name_List = Candidate_Data[0][0]
Plate_List = Candidate_Data[0]
MJD_List = Candidate_Data[0][2]
FiberID_List = Candidate_Data[0][3]
Redshift_List = Candidate_Data[0][4]
Extinction_List = Candidate_Data[0][5]
Survey_List = Candidate_Data[0][6]
run2d_List = Candidate_Data[0][7]
Standard_Inclusion_List = Candidate_Data[0][8]
Path_Override_Flag_List = Candidate_Data[0][9]
Path_Override_List = Candidate_Data[0][10]
Follow_Up_ID_List = Candidate_Data[0][11]
Smoothing_Override_List = Candidate_Data[0][12]
z_Correction_Override_List = Candidate_Data[0][13]
Extinction_Correction_Override_List = Candidate_Data[0][14]

print(Standard_Inclusion_List)

extinction = Extinction_List
z = Redshift_List

filename = f"spec-{Plate_List[1]}-{MJD_List}-{FiberID_List}.fits"

with fits.open(filename, comments='#') as hdul:
    header = hdul[0].header
    spectrum_data = hdul[1].data

"""Assumes the spectrum data is in the format:
Wavelength, Flux, FluxError"""

wave = np.array(spectrum_data['loglam'])
flux = np.array(spectrum_data['flux'])
flux_err = np.array(spectrum_data['ivar'])


lamb_observed = 10**wave * u.AA

spec_res = (lamb_observed[1] - lamb_observed[0])  # Assumes that the file has a fixed wavelength resolution

###########
# Redshift Correction
###########

lamb_rest = Hirogen_Functions.rest_wavelength_converter(observer_frame_wave=lamb_observed.value, z=z) * u.AA

flux_mean = np.mean(flux)
flux -= flux_mean

lines_to_scale = Hirogen_Functions.lines_for_scoring()

s_flux = mods_functions.flux_scaler(flux, lamb_rest, lines_to_scale, 0.5)
s_flux_err = mods_functions.flux_scaler(flux_err, lamb_rest, lines_to_scale, 1)

flux += flux_mean
s_flux += flux_mean

mods_functions.fits_file_gen(wave, s_flux, s_flux_err, filename)

flux_in = flux * u.Unit('erg cm-2 s-1 AA-1')
error = flux_err

s_flux_in = s_flux * u.Unit('erg cm-2 s-1 AA-1')
s_error = s_flux_err

all_flag_percentage = -999
counted_flag_percentage = -999

flux = flux_in

s_flux = s_flux_in

###########
# Smoothing
###########

# If active runs a boxcar smooth, defaults to a width of 5 but can be set via the function call
smoothed_flux = convolve(flux, Box1DKernel(smoothing_boxcar))
flux = smoothed_flux

s_smoothed_flux = convolve(s_flux, Box1DKernel(smoothing_boxcar))
s_flux = s_smoothed_flux

###########
# Extinction Correction
###########

# SDSS Extinction comes from the mean value of the per band extinctions converted using the constants given in the
# 2002 Early Data release paper - this uses the Schlegel, Finkbeiner and Davis dust maps

# New version of the extinction correction using the "dust_extinction" package

extinction_model = F99(Rv=3.1)
flux_extinction_corrected = flux / extinction_model.extinguish(
    lamb_observed, Ebv=extinction
)
flux = flux_extinction_corrected

s_flux_extinction_corrected = s_flux / extinction_model.extinguish(
    lamb_observed, Ebv=extinction
)
s_flux = s_flux_extinction_corrected

fig, ax = plt.subplots(figsize = (14,8))
ax.step(lamb_rest, flux, 'r')
ax.step(lamb_rest, s_flux, 'k')
ax.set_xlabel(r"Rest Wavelength ($\AA$)")
ax.set_ylabel(r"Flux (10$^{-17}$ erg/cm$^{2}$/s/$\AA$)")

plt.show()


