import matplotlib.pyplot as plt
import mpl_style
import Hirogen_Functions
import mods_functions
plt.style.use(mpl_style.scientific)

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



np.set_printoptions(threshold=sys.maxsize)
c = constants.value('speed of light in vacuum') / 1000

smoothing_boxcar = 5

scale_factors = np.arange(0.1, 1, 0.1)

User = 'Joe'
User_Config = Hirogen_Functions.user_config(user=User)

#TableID = "SDSS_Scaled_Objects"
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


Data = Hirogen_Functions.database_connection(user=Database_User, password=Database_Password, database=Database)

print("\n")
cursor = Data.cursor()


# Can be set to only run the analysis for specific IDs given in the following format
IDs = 1021306429161629696

# Manually_Inspected_Flag = -10 are spectra with external issues that prevent their use

# QUERY
cursor.execute(
    f"SELECT DR16_Spectroscopic_ID, Right_Ascension, Declination, RA_HMS, DEC_DMS, SDSS_ShortName, DR16_ParentID, DR16_Photometric_ID,"
    f"SDSS_DR16_Explore_Link, SDSS_DR16_Navigate_Link, survey, run2d, field, spec_Plate, spec_MJD, spec_FiberID, SDSS_Type, spec_class_SDSS,"
    f"spec_subclass_SDSS, spec_Human_Comments, z_SDSS_spec, z_err_SDSS_spec, z_warning_SDSS_spec, Distance_MPC, Max_Distance_MPC, Min_Distance_MPC,"
    f"Distance_Modulus, Distance_Modulus_Err, median_SNR_SDSS_spec, u_extinction, g_extinction, r_extinction, i_extinction, z_extinction, generalised_extinction,"
    f"u_petro_mag, u_petro_mag_err, u_petro_AB_mag, u_petro_AB_mag_err, u_petro_rad, u_petro_rad_err, g_petro_mag, g_petro_mag_err, g_petro_AB_mag, g_petro_AB_mag_err,"
    f"g_petro_rad, g_petro_rad_err, r_petro_mag, r_petro_mag_err, r_petro_AB_mag, r_petro_AB_mag_err, r_petro_rad, r_petro_rad_err, i_petro_mag, i_petro_mag_err,"
    f"i_petro_AB_mag, i_petro_AB_mag_err, i_petro_rad, i_petro_rad_err, z_petro_mag, z_petro_mag_err, z_petro_AB_mag, z_petro_AB_mag_err, z_petro_rad, z_petro_rad_err,"
    f"Petro_Total_Host_Mass_Estimate, Petro_Total_Host_Mass_Estimate_err, Model_Total_Host_Mass_Estimate, Model_Total_Host_Mass_Estimate_err, u_model_mag, u_model_mag_err,"
    f"u_model_AB_mag, u_model_AB_mag_err, g_model_mag, g_model_mag_err, g_model_AB_mag, g_model_AB_mag_err, r_model_mag, r_model_mag_err, r_model_AB_mag, r_model_AB_mag_err,"
    f"i_model_mag, i_model_mag_err, i_model_AB_mag, i_model_AB_mag_err, z_model_mag, z_model_mag_err, z_model_AB_mag, z_model_AB_mag_err, u_minus_r_petro, u_minus_r_petro_observation_err,"
    f"g_minus_r_petro, g_minus_r_petro_observation_err, g_minus_i_petro, g_minus_i_petro_observation_err, u_minus_r_model, u_minus_r_model_observation_err, g_minus_r_model,"
    f"g_minus_r_model_observation_err, g_minus_i_model, g_minus_i_model_observation_err, flags, clean, ABS_Petro_i_Expected_Mag_Range, Manually_Inspected_Flag, Follow_Up_ID,"
    f"Standard_Inclusion, Path_Override_Flag, Path_Override, Smoothing_Override, z_corrected_flag, extinction_corrected_flag "
    f"FROM `{Database}`.`{Objects_TableID}`"
    #f"WHERE Manually_Inspected_Flag != -10 AND lin_con_LineFlux_Halpha is NULL"
    #f"WHERE lin_con_pEQW_Halpha is NULL"
    #f"WHERE z_SDSS_spec >= 0.2"
    #f"WHERE Follow_Up_ID = 0"
    # f"WHERE Follow_Up_ID in (1,2,3)"
    #f"WHERE Standard_Inclusion = 1"
    #f"WHERE DR16_Spectroscopic_ID = {IDs}"
)

Candidate_Data = cursor.fetchall()
# Candidate_Data = cursor.fetchmany(200)

if len(Candidate_Data) >= 1:

    DR16_SpectroscopicID = [item[0] for item in Candidate_Data]
    RA = [item[1] for item in Candidate_Data]
    DEC = [item[2] for item in Candidate_Data]
    RA_HMS = [item[3] for item in Candidate_Data]
    DEC_DMS = [item[4] for item in Candidate_Data]
    SDSS_Shortnames = [item[5] for item in Candidate_Data]
    DR16_ParentID = [item[6] for item in Candidate_Data]
    DR16_PhotometricID = [item[7] for item in Candidate_Data]
    SDSS_Explore_Links = [item[8] for item in Candidate_Data]
    SDSS_Navigate_Links = [item[9] for item in Candidate_Data]
    Survey = [item[10] for item in Candidate_Data]
    Run2D = [item[11] for item in Candidate_Data]
    Field = [item[12] for item in Candidate_Data]
    Plate = [item[13] for item in Candidate_Data]
    MJD = [item[14] for item in Candidate_Data]
    FiberID = [item[15] for item in Candidate_Data]
    SDSS_Type = [item[16] for item in Candidate_Data]
    Spec_Classification = [item[17] for item in Candidate_Data]
    Spec_Sub_Classification = [item[18] for item in Candidate_Data]
    Spec_Human_Comments = [item[19] for item in Candidate_Data]
    Spec_z = [item[20] for item in Candidate_Data]
    Spec_z_err = [item[21] for item in Candidate_Data]
    Spec_z_warning = [item[22] for item in Candidate_Data]
    Distance = [item[23] for item in Candidate_Data]
    MaxDistance = [item[24] for item in Candidate_Data]
    MinDistance = [item[25] for item in Candidate_Data]
    Distance_Modulus = [item[26] for item in Candidate_Data]
    Distance_Modulus_Err = [item[27] for item in Candidate_Data]
    Median_SNR = [item[28] for item in Candidate_Data]
    Extinction_u = [item[29] for item in Candidate_Data]
    Extinction_g = [item[30] for item in Candidate_Data]
    Extinction_r = [item[31] for item in Candidate_Data]
    Extinction_i = [item[32] for item in Candidate_Data]
    Extinction_z = [item[33] for item in Candidate_Data]
    Mean_E_BminusV = [item[34] for item in Candidate_Data]
    Petrosian_u = [item[35] for item in Candidate_Data]
    Petrosian_u_err = [item[36] for item in Candidate_Data]
    Petrosian_u_AB = [item[37] for item in Candidate_Data]
    Petrosian_u_AB_err = [item[38] for item in Candidate_Data]
    Petrosian_u_rad = [item[39] for item in Candidate_Data]
    Petrosian_u_rad_err = [item[40] for item in Candidate_Data]
    Petrosian_g = [item[41] for item in Candidate_Data]
    Petrosian_g_err = [item[42] for item in Candidate_Data]
    Petrosian_g_AB = [item[43] for item in Candidate_Data]
    Petrosian_g_AB_err = [item[44] for item in Candidate_Data]
    Petrosian_g_rad = [item[45] for item in Candidate_Data]
    Petrosian_g_rad_err = [item[46] for item in Candidate_Data]
    Petrosian_r = [item[47] for item in Candidate_Data]
    Petrosian_r_err = [item[48] for item in Candidate_Data]
    Petrosian_r_AB = [item[49] for item in Candidate_Data]
    Petrosian_r_AB_err = [item[50] for item in Candidate_Data]
    Petrosian_r_rad = [item[51] for item in Candidate_Data]
    Petrosian_r_rad_err = [item[52] for item in Candidate_Data]
    Petrosian_i = [item[53] for item in Candidate_Data]
    Petrosian_i_err = [item[54] for item in Candidate_Data]
    Petrosian_i_AB = [item[55] for item in Candidate_Data]
    Petrosian_i_AB_err = [item[56] for item in Candidate_Data]
    Petrosian_i_rad = [item[57] for item in Candidate_Data]
    Petrosian_i_rad_err = [item[58] for item in Candidate_Data]
    Petrosian_z = [item[59] for item in Candidate_Data]
    Petrosian_z_err = [item[60] for item in Candidate_Data]
    Petrosian_z_AB = [item[61] for item in Candidate_Data]
    Petrosian_z_AB_err = [item[62] for item in Candidate_Data]
    Petrosian_z_rad = [item[63] for item in Candidate_Data]
    Petrosian_z_rad_err = [item[64] for item in Candidate_Data]
    TotalHostMassEstimate_Petro = [item[65] for item in Candidate_Data]
    TotalHostMassEstimate_Err_Petro = [item[66] for item in Candidate_Data]
    TotalHostMassEstimate_Model = [item[67] for item in Candidate_Data]
    TotalHostMassEstimate_Err_Model = [item[68] for item in Candidate_Data]
    Model_u = [item[69] for item in Candidate_Data]
    Model_u_err = [item[70] for item in Candidate_Data]
    Model_u_AB = [item[71] for item in Candidate_Data]
    Model_u_AB_err = [item[72] for item in Candidate_Data]
    Model_g = [item[73] for item in Candidate_Data]
    Model_g_err = [item[74] for item in Candidate_Data]
    Model_g_AB = [item[75] for item in Candidate_Data]
    Model_g_AB_err = [item[76] for item in Candidate_Data]
    Model_r = [item[77] for item in Candidate_Data]
    Model_r_err = [item[78] for item in Candidate_Data]
    Model_r_AB = [item[79] for item in Candidate_Data]
    Model_r_AB_err = [item[80] for item in Candidate_Data]
    Model_i = [item[81] for item in Candidate_Data]
    Model_i_err = [item[82] for item in Candidate_Data]
    Model_i_AB = [item[83] for item in Candidate_Data]
    Model_i_AB_err = [item[84] for item in Candidate_Data]
    Model_z = [item[85] for item in Candidate_Data]
    Model_z_err = [item[86] for item in Candidate_Data]
    Model_z_AB = [item[87] for item in Candidate_Data]
    Model_z_AB_err = [item[88] for item in Candidate_Data]
    Petrosian_u_minus_r = [item[89] for item in Candidate_Data]
    Petrosian_u_minus_r_observation_err = [item[90] for item in Candidate_Data]
    Petrosian_g_minus_r = [item[91] for item in Candidate_Data]
    Petrosian_g_minus_r_observation_err = [item[92] for item in Candidate_Data]
    Petrosian_g_minus_i = [item[93] for item in Candidate_Data]
    Petrosian_g_minus_i_observation_err = [item[94] for item in Candidate_Data]
    Model_u_minus_r = [item[95] for item in Candidate_Data]
    Model_u_minus_r_observation_err = [item[96] for item in Candidate_Data]
    Model_g_minus_r = [item[97] for item in Candidate_Data]
    Model_g_minus_r_observation_err = [item[98] for item in Candidate_Data]
    Model_g_minus_i = [item[99] for item in Candidate_Data]
    Model_g_minus_i_observation_err = [item[100] for item in Candidate_Data]
    Spec_Flags = [item[101] for item in Candidate_Data]
    Clean_Flag = [item[102] for item in Candidate_Data]
    Petro_AB_i_Check = [item[103] for item in Candidate_Data]
    Manual_Flag = [item[104] for item in Candidate_Data]
    Follow_UpID = [item[105] for item in Candidate_Data]
    Standard_Inclusion = [item[106] for item in Candidate_Data]
    Path_Override_Flag = [item[107] for item in Candidate_Data]
    Path_Override = [item[108] for item in Candidate_Data]
    Smoothing_Override = [item[109] for item in Candidate_Data]
    z_Corrected_Flag = [item[110] for item in Candidate_Data]
    Extinction_Corrected_Flag = [item[111] for item in Candidate_Data]

cursor.execute(f"SELECT * FROM {Database}.{Spectra_TableID}")
rows = cursor.fetchall()

print(f'Number of rows prior to row entry is: {len(rows)}')

filenames = Hirogen_Functions.sdss_spectra_file_path_generator(
    Main_Spectra_Path, Plate, MJD, FiberID, Survey, Run2D, 
    Path_Override_Flag, Path_Override
)

for i, object in enumerate(Candidate_Data):

    with fits.open(filenames[i], comments='#') as hdul:
            header = hdul[0].header
            spectrum_data = hdul[1].data

    #Assumes the spectrum data is in the format:
    #Wavelength, Flux, FluxError

    wave = np.array(spectrum_data['loglam'])
    flux = np.array(spectrum_data['flux'])
    flux_err = np.array(spectrum_data['ivar'])

    flux_to_scale = flux

    lamb_observed = 10**wave * u.AA

    spec_res = (lamb_observed[1] - lamb_observed[0])  # Assumes that the file has a fixed wavelength resolution

    ###########
    # Redshift Correction
    ###########

    lamb_rest = Hirogen_Functions.rest_wavelength_converter(observer_frame_wave=lamb_observed.value, z=Spec_z[i]) * u.AA

    lines_to_scale = Hirogen_Functions.lines_for_scaling()

    for j, scale_factor in enumerate(scale_factors):

        s_flux = mods_functions.flux_scaler(flux_to_scale, lamb_rest, lines_to_scale, scale_factor)
        s_flux_err = mods_functions.flux_scaler(flux_err, lamb_rest, lines_to_scale, scale_factor)

        save_name = f"spec-{Plate[i]}-{MJD[i]}-{FiberID[i]}-{int(scale_factor*10)}.fits"

        
        scale_ID = DR16_SpectroscopicID[i] + decimal.Decimal(scale_factor)
        

        #mods_functions.fits_file_gen(wave, s_flux, s_flux_err, save_name, scale_factor)

        cmd = f"INSERT IGNORE INTO `{Database}`.`{Spectra_TableID}` (DR16_Spectroscopic_ID) " \
                f"VALUES ({scale_ID});"

        cursor.execute(cmd)

        cmd = f"UPDATE {Database}.{Spectra_TableID} " \
            f"set Right_Ascension = '{RA[i]}'," \
            f"Declination = '{DEC[i]}'," \
            f"RA_HMS = '{RA_HMS[i]}'," \
            f"DEC_DMS = '{DEC_DMS[i]}'," \
            f"SDSS_ShortName = '{SDSS_Shortnames[i]}'," \
            f"DR16_ParentID = '{DR16_ParentID[i]}'," \
            f"DR16_Photometric_ID = '{DR16_PhotometricID[i]}'," \
            f"SDSS_DR16_Explore_Link ='{SDSS_Explore_Links[i]}'," \
            f"SDSS_DR16_Navigate_Link ='{SDSS_Navigate_Links[i]}'," \
            f"survey = '{Survey[i]}'," \
            f"run2d = '{Run2D[i]}'," \
            f"field = '{Field[i]}'," \
            f"spec_Plate = '{Plate[i]}'," \
            f"spec_MJD = '{MJD[i]}' ," \
            f"spec_FiberID = '{FiberID[i]}'," \
            f"SDSS_Type = '{SDSS_Type[i]}'," \
            f"spec_class_SDSS = '{Spec_Classification[i]}'," \
            f"spec_subclass_SDSS = '{Spec_Sub_Classification[i]}'," \
            f"spec_Human_Comments = '{Spec_Human_Comments[i]}'," \
            f"z_SDSS_spec = '{Spec_z[i]}'," \
            f"z_err_SDSS_spec = '{Spec_z_err[i]}'," \
            f"z_warning_SDSS_spec = '{Spec_z_warning[i]}'," \
            f"Distance_MPC = '{Distance[i]}'," \
            f"Max_Distance_MPC = '{MaxDistance[i]}'," \
            f"Min_Distance_MPC = '{MinDistance[i]}'," \
            f"Distance_Modulus='{Distance_Modulus[i]}'," \
            f"Distance_Modulus_Err='{Distance_Modulus_Err[i]}'," \
            f"median_SNR_SDSS_spec = '{Median_SNR[i]}'," \
            f"u_extinction = '{Extinction_u[i]}'," \
            f"g_extinction = '{Extinction_g[i]}'," \
            f"r_extinction = '{Extinction_r[i]}'," \
            f"i_extinction = '{Extinction_i[i]}'," \
            f"z_extinction = '{Extinction_z[i]}'," \
            f"generalised_extinction = '{Mean_E_BminusV[i]}'," \
            f"u_petro_mag = '{Petrosian_u[i]}'," \
            f"u_petro_mag_err = '{Petrosian_u_err[i]}'," \
            f"u_petro_AB_mag = '{Petrosian_u_AB[i]}'," \
            f"u_petro_AB_mag_err = '{Petrosian_u_AB_err[i]}'," \
            f"u_petro_rad= '{Petrosian_u_rad[i]}'," \
            f"u_petro_rad_err= '{Petrosian_u_rad_err[i]}'," \
            f"g_petro_mag = '{Petrosian_g[i]}'," \
            f"g_petro_mag_err = '{Petrosian_g_err[i]}'," \
            f"g_petro_AB_mag = '{Petrosian_g_AB[i]}'," \
            f"g_petro_AB_mag_err = '{Petrosian_g_AB_err[i]}'," \
            f"g_petro_rad= '{Petrosian_g_rad[i]}'," \
            f"g_petro_rad_err= '{Petrosian_g_rad_err[i]}'," \
            f"r_petro_mag = '{Petrosian_r[i]}'," \
            f"r_petro_mag_err = '{Petrosian_r_err[i]}'," \
            f"r_petro_AB_mag = '{Petrosian_r_AB[i]}'," \
            f"r_petro_AB_mag_err = '{Petrosian_r_AB_err[i]}'," \
            f"r_petro_rad= '{Petrosian_r_rad[i]}'," \
            f"r_petro_rad_err= '{Petrosian_r_rad_err[i]}'," \
            f"i_petro_mag = '{Petrosian_i[i]}'," \
            f"i_petro_mag_err = '{Petrosian_i_err[i]}'," \
            f"i_petro_AB_mag = '{Petrosian_i_AB[i]}'," \
            f"i_petro_AB_mag_err = '{Petrosian_i_AB_err[i]}'," \
            f"i_petro_rad= '{Petrosian_i_rad[i]}'," \
            f"i_petro_rad_err= '{Petrosian_i_rad_err[i]}'," \
            f"z_petro_mag = '{Petrosian_z[i]}'," \
            f"z_petro_mag_err = '{Petrosian_z_err[i]}'," \
            f"z_petro_AB_mag = '{Petrosian_z_AB[i]}'," \
            f"z_petro_AB_mag_err = '{Petrosian_z_AB_err[i]}'," \
            f"z_petro_rad= '{Petrosian_z_rad[i]}'," \
            f"z_petro_rad_err= '{Petrosian_z_rad_err[i]}'," \
            f"Petro_Total_Host_Mass_Estimate = '{TotalHostMassEstimate_Petro[i]}'," \
            f"Petro_Total_Host_Mass_Estimate_err = '{TotalHostMassEstimate_Err_Petro[i]}'," \
            f"Model_Total_Host_Mass_Estimate = '{TotalHostMassEstimate_Model[i]}'," \
            f"Model_Total_Host_Mass_Estimate_err = '{TotalHostMassEstimate_Err_Model[i]}'," \
            f"u_model_mag = '{Model_u[i]}'," \
            f"u_model_mag_err = '{Model_u_err[i]}'," \
            f"u_model_AB_mag = '{Model_u_AB[i]}'," \
            f"u_model_AB_mag_err = '{Model_u_AB_err[i]}'," \
            f"g_model_mag = '{Model_g[i]}'," \
            f"g_model_mag_err = '{Model_g_err[i]}'," \
            f"g_model_AB_mag = '{Model_g_AB[i]}'," \
            f"g_model_AB_mag_err = '{Model_g_AB_err[i]}'," \
            f"r_model_mag = '{Model_r[i]}'," \
            f"r_model_mag_err = '{Model_r_err[i]}'," \
            f"r_model_AB_mag = '{Model_r_AB[i]}'," \
            f"r_model_AB_mag_err = '{Model_r_AB_err[i]}'," \
            f"i_model_mag = '{Model_i[i]}'," \
            f"i_model_mag_err = '{Model_i_err[i]}'," \
            f"i_model_AB_mag = '{Model_i_AB[i]}'," \
            f"i_model_AB_mag_err = '{Model_i_AB_err[i]}'," \
            f"z_model_mag = '{Model_z[i]}'," \
            f"z_model_mag_err = '{Model_z_err[i]}'," \
            f"z_model_AB_mag = '{Model_z_AB[i]}'," \
            f"z_model_AB_mag_err = '{Model_z_AB_err[i]}'," \
            f"u_minus_r_petro = '{Petrosian_u_minus_r[i]}'," \
            f"u_minus_r_petro_observation_err = '{Petrosian_u_minus_r_observation_err[i]}'," \
            f"g_minus_r_petro = '{Petrosian_g_minus_r[i]}'," \
            f"g_minus_r_petro_observation_err = '{Petrosian_g_minus_r_observation_err[i]}'," \
            f"g_minus_i_petro = '{Petrosian_g_minus_i[i]}'," \
            f"g_minus_i_petro_observation_err = '{Petrosian_g_minus_i_observation_err[i]}'," \
            f"u_minus_r_model = '{Model_u_minus_r[i]}'," \
            f"u_minus_r_model_observation_err = '{Model_u_minus_r_observation_err[i]}'," \
            f"g_minus_r_model = '{Model_g_minus_r[i]}'," \
            f"g_minus_r_model_observation_err = '{Model_g_minus_r_observation_err[i]}'," \
            f"g_minus_i_model = '{Model_g_minus_i[i]}'," \
            f"g_minus_i_model_observation_err = '{Model_g_minus_i_observation_err[i]}'," \
            f"flags = '{Spec_Flags[i]}'," \
            f"clean = '{Clean_Flag[i]}'," \
            f"ABS_Petro_i_Expected_Mag_Range = '{Petro_AB_i_Check[i]}'," \
            f"Manually_Inspected_Flag = '{Manual_Flag[i]}'," \
            f"Follow_Up_ID = '{Follow_UpID[i]}'," \
            f"Standard_Inclusion = '{Standard_Inclusion[i]}'," \
            f"Path_Override_Flag= '{Path_Override_Flag[i]}'," \
            f"Path_Override = '{Path_Override[i]}'," \
            f"Smoothing_Override = '{Smoothing_Override[i]}'," \
            f"z_corrected_flag = '{z_Corrected_Flag[i]}'," \
            f"extinction_corrected_flag = '{Extinction_Corrected_Flag[i]}'," \
            f"Scale_Factor = '{scale_factor}'" \
            f"WHERE DR16_Spectroscopic_ID = '{scale_ID}' "

        cursor.execute(cmd)

cursor.execute(f"SELECT DR16_Spectroscopic_ID FROM {Database}.{Spectra_TableID}")
rows = cursor.fetchall()
print(f'Number of rows is now: {len(rows)}')

Hirogen_Functions.commit_and_close(Data, cursor)