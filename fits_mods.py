import matplotlib.pyplot as plt
import mpl_style
import Hirogen_Functions
plt.style.use(mpl_style.style1)

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

extinction = Extinction_List
z = Redshift_List

hdul = fits.open('spec-0907-52373-0419.fits', comments='#')

header = hdul[0].header
spectrum_data = hdul[1].data

"""Assumes the spectrum data is in the format:
Wavelength, Flux, FluxError"""

wave = np.array(spectrum_data['loglam'])
flux = np.array(spectrum_data['flux'])
flux_err = np.array(spectrum_data['ivar'])

lamb_observed = 10**wave * u.AA

spec_res = (lamb_observed[1] - lamb_observed[0])  # Assumes that the file has a fixed wavelength resolution

flux_in = flux * u.Unit('erg cm-2 s-1 AA-1')
error = flux_err

all_flag_percentage = -999
counted_flag_percentage = -999

flux = flux_in



###########
# Smoothing
###########

# If active runs a boxcar smooth, defaults to a width of 5 but can be set via the function call
smoothed_flux = convolve(flux, Box1DKernel(smoothing_boxcar))
flux = smoothed_flux

###########
# Median Filter
###########

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

###########
# Redshift Correction
###########

lamb_rest = Hirogen_Functions.rest_wavelength_converter(observer_frame_wave=lamb_observed.value, z=z) * u.AA

fig, ax = plt.subplots(figsize = (14,8))
ax.step(lamb_rest, flux, 'k')
ax.set_xlabel(r"Rest Wavelength ($\AA$)")
ax.set_ylabel(r"Flux (10$^{-17}$ erg/cm$^{2}$/s/$\AA$)")

plt.show()