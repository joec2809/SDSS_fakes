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
from tqdm.notebook import tqdm

# This code takes the 7 confirmed ECLE objects, produces 9 scaled copies of each and saves them as fits files
# It also adds the corresponding info to an SQL table that contains other galaxy spectra

np.set_printoptions(threshold=sys.maxsize)
c = constants.value('speed of light in vacuum') / 1000

smoothing_boxcar = 5

scale_factors = np.around(np.arange(0.1, 1.1, 0.1), 1)

User = 'Joe'
User_Config = Hirogen_Functions.user_config(user=User)

Objects_TableID = "SDSS_Confirmed_Objects"
Spectra_TableID = 'SDSS_Galaxy_Spectra'

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

coronal_lines = Hirogen_Functions.coronal_lines()


Data = Hirogen_Functions.database_connection(user=Database_User, password=Database_Password, database=Database)

print("\n")
cursor = Data.cursor()

cursor.execute(
        f"SELECT DR16_Spectroscopic_ID, spec_Plate, spec_MJD, spec_FiberID, z_SDSS_spec, generalised_extinction, "
        f"survey, run2d, Standard_Inclusion,Path_Override_Flag, Path_Override, Follow_Up_ID, Smoothing_Override,"
        f"z_corrected_flag, extinction_corrected_flag "
        f"FROM `{Database}`.`{Spectra_TableID}`"
        #f"WHERE Manually_Inspected_Flag != -10 AND lin_con_LineFlux_Halpha is NULL"
        #f"WHERE lin_con_pEQW_Halpha is NULL"
        #f"WHERE z_SDSS_spec >= 0.2"
        #f"WHERE Follow_Up_ID = 0"
        #f"WHERE Follow_Up_ID in (1,2,3)"
        #f"WHERE Nickname = 'Leafeon' AND Follow_Up_ID in (2,3)"
        #f"WHERE Standard_Inclusion = 1"
        #f"WHERE DR16_Spectroscopic_ID IN {IDs} and Standard_Inclusion = 1"
        #f"WHERE DR16_Spectroscopic_ID = {IDs}"
    )

galaxy_data = cursor.fetchall()
# galaxy_data = cursor.fetchmany(200)

#print(galaxy_data)

if len(galaxy_data) >= 1:

    galaxy_Object_Name_List = [item[0] for item in galaxy_data]
    galaxy_Plate_List = [item[1] for item in galaxy_data]
    galaxy_MJD_List = [item[2] for item in galaxy_data]
    galaxy_FiberID_List = [item[3] for item in galaxy_data]
    galaxy_Redshift_List = [item[4] for item in galaxy_data]
    galaxy_Extinction_List = [item[5] for item in galaxy_data]
    galaxy_Survey_List = [item[6] for item in galaxy_data]
    galaxy_run2d_List = [item[7] for item in galaxy_data]
    galaxy_Standard_Inclusion_List = [item[8] for item in galaxy_data]
    galaxy_Path_Override_Flag_List = [item[9] for item in galaxy_data]
    galaxy_Path_Override_List = [item[10] for item in galaxy_data]
    galaxy_Follow_Up_ID_List = [item[11] for item in galaxy_data]
    galaxy_Smoothing_Override_List = [item[12] for item in galaxy_data]
    galaxy_z_Correction_Override_List = [item[13] for item in galaxy_data]
    galaxy_Extinction_Correction_Override_List = [item[14] for item in galaxy_data]

else:
    print("galaxy_data Length error: Check and try again.")

    sys.exit()

galaxy_FilePaths = Hirogen_Functions.sdss_spectra_file_path_generator(
    Main_Spectra_Path, galaxy_Plate_List, galaxy_MJD_List, galaxy_FiberID_List, galaxy_Survey_List, galaxy_run2d_List, 
    galaxy_Path_Override_Flag_List, galaxy_Path_Override_List
)

cursor.execute(
    f"SELECT DR16_Spectroscopic_ID, spec_Plate, spec_MJD, spec_FiberID, z_SDSS_spec, generalised_extinction, "
    f"survey, run2d, Standard_Inclusion,Path_Override_Flag, Path_Override, Follow_Up_ID, Smoothing_Override,"
    f"z_corrected_flag, extinction_corrected_flag "
    f"FROM `{Database}`.`{Objects_TableID}`"
    #f"WHERE Manually_Inspected_Flag != -10 AND lin_con_LineFlux_Halpha is NULL"
    #f"WHERE lin_con_pEQW_Halpha is NULL"
    #f"WHERE z_SDSS_spec >= 0.2"
    #f"WHERE Follow_Up_ID = 0"
    #f"WHERE Follow_Up_ID in (1,2,3)"
    #f"WHERE Nickname = 'Leafeon' AND Follow_Up_ID in (2,3)"
    #f"WHERE Standard_Inclusion = 1"
    #f"WHERE DR16_Spectroscopic_ID IN {IDs} and Standard_Inclusion = 1"
    #f"WHERE DR16_Spectroscopic_ID = {IDs}"
)

ecle_data = cursor.fetchall()
# ecle_data = cursor.fetchmany(200)

#print(ecle_data)

if len(ecle_data) >= 1:

    ecle_Object_Name_List = [item[0] for item in ecle_data]
    ecle_Plate_List = [item[1] for item in ecle_data]
    ecle_MJD_List = [item[2] for item in ecle_data]
    ecle_FiberID_List = [item[3] for item in ecle_data]
    ecle_Redshift_List = [item[4] for item in ecle_data]
    ecle_Extinction_List = [item[5] for item in ecle_data]
    ecle_Survey_List = [item[6] for item in ecle_data]
    ecle_run2d_List = [item[7] for item in ecle_data]
    ecle_Standard_Inclusion_List = [item[8] for item in ecle_data]
    ecle_Path_Override_Flag_List = [item[9] for item in ecle_data]
    ecle_Path_Override_List = [item[10] for item in ecle_data]
    ecle_Follow_Up_ID_List = [item[11] for item in ecle_data]
    ecle_Smoothing_Override_List = [item[12] for item in ecle_data]
    ecle_z_Correction_Override_List = [item[13] for item in ecle_data]
    ecle_Extinction_Correction_Override_List = [item[14] for item in ecle_data]

else:
    print("ecle_data Length error: Check and try again.")

    sys.exit()

peaks_FilePaths = Hirogen_Functions.sdss_peaks_spectra_file_path_generator(
    Main_Spectra_Path, ecle_Plate_List, ecle_MJD_List, ecle_FiberID_List, ecle_Survey_List, ecle_run2d_List, 
    ecle_Path_Override_Flag_List, ecle_Path_Override_List
)

print(len(peaks_FilePaths))

for i, spectrum in tqdm(enumerate(galaxy_Object_Name_List)):

    with fits.open(galaxy_FilePaths[i], comments='#') as hdul:
            galaxy_header = hdul[0].header
            galaxy_spectrum_data = hdul[1].data

    #Assumes the spectrum data is in the format:
    #Wavelength, Flux, FluxError

    galaxy_wave = np.array(galaxy_spectrum_data['loglam'])
    galaxy_flux = np.array(galaxy_spectrum_data['flux'])
    galaxy_flux_err = np.array(galaxy_spectrum_data['ivar'])
    galaxy_flags = np.array(galaxy_spectrum_data['and_mask'])

    lamb_observed = 10**galaxy_wave * u.AA

    spec_res = (lamb_observed[1] - lamb_observed[0])  # Assumes that the file has a fixed wavelength resolution

    ###########
    # Redshift Correction
    ###########

    galaxy_lamb_rest = Hirogen_Functions.rest_wavelength_converter(observer_frame_wave=lamb_observed.value, z=galaxy_Redshift_List[i]) * u.AA

    galaxy_flux_with_continua = mods_functions.peak_remover(galaxy_flux, galaxy_lamb_rest, coronal_lines)

    peak_file_idx = np.random.randint(0,7)

    peak_filepath = peaks_FilePaths[peak_file_idx]

    with fits.open(peak_filepath, comments='#') as hdul:
            peak_header = hdul[0].header
            peak_spectrum_data = hdul[1].data

    #Assumes the spectrum data is in the format:
    #Wavelength, Flux, FluxError

    peak_wave = np.array(peak_spectrum_data['loglam'])
    peak_flux = np.array(peak_spectrum_data['flux'])
    peak_flux_err = np.array(peak_spectrum_data['ivar'])
    peak_flags = np.array(peak_spectrum_data['and_mask'])

    lamb_observed = 10**peak_wave * u.AA

    scaled_peak_flux = peak_flux * np.random.random()

    spec_res = (lamb_observed[1] - lamb_observed[0])  # Assumes that the file has a fixed wavelength resolution

    ###########
    # Redshift Correction
    ###########

    peak_lamb_rest = Hirogen_Functions.rest_wavelength_converter(observer_frame_wave=lamb_observed.value, z=ecle_Redshift_List[peak_file_idx]) * u.AA

    peak_wavelengths, peak_flux, galaxy_wavelengths, galaxy_flux = mods_functions.spectra_resizer(peak_lamb_rest, scaled_peak_flux,
                                                                                                galaxy_lamb_rest, galaxy_flux_with_continua)

    fake_flux = peak_flux + galaxy_flux

    new_flux_err = np.sqrt(galaxy_flux**2 + peak_flux**2)

    save_wave = np.log10((1+ecle_Redshift_List[peak_file_idx])*galaxy_wavelengths.value)

    mods_functions.fits_file_gen(save_wave, fake_flux, new_flux_err, galaxy_flags, Main_Spectra_Path, galaxy_Plate_List[i], galaxy_MJD_List[i], galaxy_FiberID_List[i],
                                galaxy_Survey_List[i], galaxy_run2d_List[i], galaxy_Path_Override_Flag_List[i], galaxy_Path_Override_List[i], 'fakes')