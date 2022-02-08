import Hirogen_Functions
import mods_functions

import sys
import os
import math

import numpy as np
import scipy.constants as constants
import mysql.connector
import decimal
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy import units as u
from mysql.connector import errorcode

# This code takes the 7 confirmed ECLE objects, produces 9 scaled copies of each and saves them as fits files
# It also adds the corresponding info to an SQL table that contains other galaxy spectra

np.set_printoptions(threshold=sys.maxsize)
c = constants.value('speed of light in vacuum') / 1000

smoothing_boxcar = 5

User = 'Joe'
User_Config = Hirogen_Functions.user_config(user=User)

TableID = "SDSS_Confirmed_Spectra"

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

cursor.execute(
        f"SELECT DR16_Spectroscopic_ID, spec_Plate, spec_MJD, spec_FiberID, z_SDSS_spec, generalised_extinction, "
        f"survey, run2d, Standard_Inclusion,Path_Override_Flag, Path_Override, Follow_Up_ID, Smoothing_Override,"
        f"z_corrected_flag, extinction_corrected_flag "
        f"FROM `{Database}`.`{TableID}`"
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

Candidate_Data = cursor.fetchall()
# Candidate_Data = cursor.fetchmany(200)

#print(Candidate_Data)

if len(Candidate_Data) >= 1:

    Object_Name_List = [item[0] for item in Candidate_Data]
    Plate_List = [item[1] for item in Candidate_Data]
    MJD_List = [item[2] for item in Candidate_Data]
    FiberID_List = [item[3] for item in Candidate_Data]
    Redshift_List = [item[4] for item in Candidate_Data]
    Extinction_List = [item[5] for item in Candidate_Data]
    Survey_List = [item[6] for item in Candidate_Data]
    run2d_List = [item[7] for item in Candidate_Data]
    Standard_Inclusion_List = [item[8] for item in Candidate_Data]
    Path_Override_Flag_List = [item[9] for item in Candidate_Data]
    Path_Override_List = [item[10] for item in Candidate_Data]
    Follow_Up_ID_List = [item[11] for item in Candidate_Data]
    Smoothing_Override_List = [item[12] for item in Candidate_Data]
    z_Correction_Override_List = [item[13] for item in Candidate_Data]
    Extinction_Correction_Override_List = [item[14] for item in Candidate_Data]

else:
    print("Candidate_Data Length error: Check and try again.")

    sys.exit()


FilePaths = Hirogen_Functions.sdss_spectra_file_path_generator(
    Main_Spectra_Path, Plate_List, MJD_List, FiberID_List, Survey_List, run2d_List,
    override_path_flag_list=Path_Override_Flag_List, override_path_list=Path_Override_List
)

for i, item in enumerate(Candidate_Data):

    with fits.open(FilePaths[i], comments='#') as hdul:
            header = hdul[0].header
            spectrum_data = hdul[1].data

    #Assumes the spectrum data is in the format:
    #Wavelength, Flux, FluxError

    wave = np.array(spectrum_data['loglam'])
    flux = np.array(spectrum_data['flux'])
    flux_err = np.array(spectrum_data['ivar'])
    flags = np.array(spectrum_data['and_mask'])

    lamb_observed = 10**wave * u.AA

    spec_res = (lamb_observed[1] - lamb_observed[0])  # Assumes that the file has a fixed wavelength resolution

    ###########
    # Redshift Correction
    ###########

    lamb_rest = Hirogen_Functions.rest_wavelength_converter(observer_frame_wave=lamb_observed.value, z=Redshift_List[i]) * u.AA

    coronal_lines = Hirogen_Functions.coronal_lines()

    # Set most of spectra flux to zero, leaving the peaks only

    peaks_only = mods_functions.peak_selector(flux, lamb_rest, coronal_lines)

    # Error propagation of subtracting flux

    new_flux_err = np.sqrt(2)*flux_err

    # Save as new fits file

    mods_functions.fits_file_gen(wave, peaks_only, new_flux_err, flags, Main_Spectra_Path,
        Plate_List[i], MJD_List[i], FiberID_List[i], Survey_List[i], run2d_List[i],
        Path_Override_Flag_List[i], Path_Override_List[i], 'peaks'
    )