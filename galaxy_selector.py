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
from random import randint

User = 'Joe'
Input_TableID = "SDSS_Galaxy_Spectra"
Output_TableID = 'SDSS_Fake_Spectra'

number_to_output = 50

User_Config = Hirogen_Functions.user_config(user=User)

Database = User_Config[0]
Database_User = User_Config[1]
Database_Password = User_Config[2]

Data = Hirogen_Functions.database_connection(user=Database_User, password=Database_Password, database=Database)

print("\n")
cursor = Data.cursor()

cursor.execute(
        f"SELECT DR16_Spectroscopic_ID, SDSS_ShortName, SDSS_DR16_Explore_Link, SDSS_DR16_Navigate_Link, survey, Right_Ascension, Declination, RA_HMS, DEC_DMS, "
        f"DR16_ParentID, DR7_Photometric_ID, DR16_Photometric_ID, median_SNR_SDSS_spec, z_SDSS_spec, z_err_SDSS_spec, z_warning_SDSS_spec, Distance_MPC, Max_Distance_MPC, "
        f"Min_Distance_MPC, Distance_Modulus, Distance_Modulus_Err, run2d, field, spec_MJD, spec_Plate, spec_FiberID, SDSS_Type, spec_class_SDSS, "
        f"spec_subclass_SDSS, spec_Human_comments, u_extinction, g_extinction, r_extinction, i_extinction, z_extinction, generalised_extinction, "
        f"clean, flags, u_petro_mag, u_petro_mag_err, u_petro_AB_mag, u_petro_AB_mag_err, u_petro_rad, u_petro_rad_err, g_petro_mag, g_petro_mag_err, g_petro_AB_mag, g_petro_AB_mag_err, g_petro_rad, g_petro_rad_err, "
        f"r_petro_mag, r_petro_mag_err, r_petro_AB_mag, r_petro_AB_mag_err, r_petro_rad, r_petro_rad_err, i_petro_mag, i_petro_mag_err, i_petro_AB_mag, "
        f"i_petro_AB_mag_err, i_petro_rad, i_petro_rad_err, z_petro_mag, z_petro_mag_err, z_petro_AB_mag, z_petro_AB_mag_err, z_petro_rad, z_petro_rad_err, "
        f"u_minus_r_petro, u_minus_r_petro_observation_err, g_minus_r_petro, g_minus_r_petro_observation_err, g_minus_i_petro, g_minus_i_petro_observation_err, "
        f"u_model_mag, u_model_mag_err, u_model_AB_mag, u_model_AB_mag_err, g_model_mag, g_model_mag_err, g_model_AB_mag, g_model_AB_mag_err, r_model_mag, "
        f"r_model_mag_err, r_model_AB_mag, r_model_AB_mag_err, i_model_mag, i_model_mag_err, i_model_AB_mag, i_model_AB_mag_err, z_model_mag, z_model_mag_err, "
        f"z_model_AB_mag, z_model_AB_mag_err, u_minus_r_model, u_minus_r_model_observation_err, g_minus_r_model, g_minus_r_model_observation_err, g_minus_i_model, "
        f"g_minus_i_model_observation_err, Petro_Total_Host_Mass_Estimate, Petro_Total_Host_Mass_Estimate_err, Model_Total_Host_Mass_Estimate, Model_Total_Host_Mass_Estimate_err "
        f"FROM `{Database}`.`{Input_TableID}`"
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
    SDSS_Shortnames = [item[1] for item in Candidate_Data]
    SDSS_Explore_Links = [item[2] for item in Candidate_Data]
    SDSS_Navigate_Links = [item[3] for item in Candidate_Data]
    Survey = [item[4] for item in Candidate_Data]
    RA = [item[5] for item in Candidate_Data]
    DEC = [item[6] for item in Candidate_Data]
    RA_HMS = [item[7] for item in Candidate_Data]
    DEC_HMS = [item[8] for item in Candidate_Data]
    DR16_ParentID = [item[9] for item in Candidate_Data]
    DR16_PhotometricID = [item[10] for item in Candidate_Data]
    Median_SNR = [item[11] for item in Candidate_Data]
    Spec_z = [item[12] for item in Candidate_Data]
    Spec_z_err = [item[13] for item in Candidate_Data]
    Spec_z_warning = [item[14] for item in Candidate_Data]
    Distance = [item[15] for item in Candidate_Data]
    MaxDistance = [item[16] for item in Candidate_Data]
    MinDistance = [item[17] for item in Candidate_Data]
    Distance_Modulus = [item[18] for item in Candidate_Data]
    Distance_Modulus_Err = [item[19] for item in Candidate_Data]
    Run2D = [item[20] for item in Candidate_Data]
    Field = [item[21] for item in Candidate_Data]
    MJD = [item[22] for item in Candidate_Data]
    Plate = [item[23] for item in Candidate_Data]
    FiberID = [item[24] for item in Candidate_Data]
    SDSS_Type = [item[25] for item in Candidate_Data]
    Spec_Classification = [item[26] for item in Candidate_Data]
    Spec_Sub_Classification = [item[27] for item in Candidate_Data]
    Spec_Human_Comments = [item[28] for item in Candidate_Data]
    Extinction_u = [item[29] for item in Candidate_Data]
    Extinction_g = [item[30] for item in Candidate_Data]
    Extinction_r = [item[31] for item in Candidate_Data]
    Extinction_i = [item[32] for item in Candidate_Data]
    Extinction_z = [item[33] for item in Candidate_Data]
    Mean_E_BminusV = [item[34] for item in Candidate_Data]
    Clean_Flag = [item[35] for item in Candidate_Data]
    Spec_Flags = [item[36] for item in Candidate_Data]
    Petrosian_u = [item[37] for item in Candidate_Data]
    Petrosian_u_err = [item[38] for item in Candidate_Data]
    Petrosian_u_AB = [item[39] for item in Candidate_Data]
    Petrosian_u_AB_err = [item[40] for item in Candidate_Data]
    Petrosian_u_rad = [item[41] for item in Candidate_Data]
    Petrosian_u_rad_err = [item[42] for item in Candidate_Data]
    Petrosian_g = [item[43] for item in Candidate_Data]
    Petrosian_g_err = [item[44] for item in Candidate_Data]
    Petrosian_g_AB = [item[45] for item in Candidate_Data]
    Petrosian_g_AB_err = [item[46] for item in Candidate_Data]
    Petrosian_g_rad = [item[47] for item in Candidate_Data]
    Petrosian_g_rad_err = [item[48] for item in Candidate_Data]
    Petrosian_r = [item[49] for item in Candidate_Data]
    Petrosian_r_err = [item[50] for item in Candidate_Data]
    Petrosian_r_AB = [item[51] for item in Candidate_Data]
    Petrosian_r_AB_err = [item[52] for item in Candidate_Data]
    Petrosian_r_rad = [item[53] for item in Candidate_Data]
    Petrosian_r_rad_err = [item[54] for item in Candidate_Data]
    Petrosian_i = [item[55] for item in Candidate_Data]
    Petrosian_i_err = [item[56] for item in Candidate_Data]
    Petrosian_i_AB = [item[57] for item in Candidate_Data]
    Petrosian_i_AB_err = [item[58] for item in Candidate_Data]
    Petrosian_i_rad = [item[59] for item in Candidate_Data]
    Petrosian_i_rad_err = [item[60] for item in Candidate_Data]
    Petrosian_z = [item[61] for item in Candidate_Data]
    Petrosian_z_err = [item[62] for item in Candidate_Data]
    Petrosian_z_AB = [item[63] for item in Candidate_Data]
    Petrosian_z_AB_err = [item[64] for item in Candidate_Data]
    Petrosian_z_rad = [item[65] for item in Candidate_Data]
    Petrosian_z_rad_err = [item[66] for item in Candidate_Data]
    Petrosian_u_minus_r = [item[67] for item in Candidate_Data]
    Petrosian_u_minus_r_observation_err = [item[68] for item in Candidate_Data]
    Petrosian_g_minus_r = [item[69] for item in Candidate_Data]
    Petrosian_g_minus_r_observation_err = [item[70] for item in Candidate_Data]
    Petrosian_g_minus_i = [item[71] for item in Candidate_Data]
    Petrosian_g_minus_i_observation_err = [item[72] for item in Candidate_Data]
    Model_u = [item[73] for item in Candidate_Data]
    Model_u_err = [item[74] for item in Candidate_Data]
    Model_u_AB = [item[75] for item in Candidate_Data]
    Model_u_AB_err = [item[76] for item in Candidate_Data]
    Model_g = [item[77] for item in Candidate_Data]
    Model_g_err = [item[78] for item in Candidate_Data]
    Model_g_AB = [item[79] for item in Candidate_Data]
    Model_g_AB_err = [item[80] for item in Candidate_Data]
    Model_r = [item[81] for item in Candidate_Data]
    Model_r_err = [item[82] for item in Candidate_Data]
    Model_r_AB = [item[83] for item in Candidate_Data]
    Model_r_AB_err = [item[84] for item in Candidate_Data]
    Model_i = [item[85] for item in Candidate_Data]
    Model_i_err = [item[86] for item in Candidate_Data]
    Model_i_AB = [item[87] for item in Candidate_Data]
    Model_i_AB_err = [item[88] for item in Candidate_Data]
    Model_z = [item[89] for item in Candidate_Data]
    Model_z_err = [item[90] for item in Candidate_Data]
    Model_z_AB = [item[91] for item in Candidate_Data]
    Model_z_AB_err = [item[92] for item in Candidate_Data]
    Model_u_minus_r= [item[93] for item in Candidate_Data]
    Model_u_minus_r_observation_err= [item[94] for item in Candidate_Data]
    Model_g_minus_r= [item[95] for item in Candidate_Data]
    Model_g_minus_r_observation_err= [item[96] for item in Candidate_Data]
    Model_g_minus_i= [item[97] for item in Candidate_Data]
    Model_g_minus_i_observation_err= [item[98] for item in Candidate_Data]
    TotalHostMassEstimate_Petro = [item[99] for item in Candidate_Data]
    TotalHostMassEstimate_Err_Petro = [item[100] for item in Candidate_Data]
    TotalHostMassEstimate_Model = [item[101] for item in Candidate_Data]
    TotalHostMassEstimate_Err_Model = [item[102] for item in Candidate_Data]

else:
    print("Candidate_Data Length error: Check and try again.")

    sys.exit()

used_spectra = []

while len(used_spectra) < number_to_output:
    j = randint(0,len(Candidate_Data))
    if j in used_spectra:
        continue
    used_spectra.append(j)

    cmd = f"INSERT IGNORE INTO `{Database}`.`{Output_TableID}` (DR16_Spectroscopic_ID) " \
            f"VALUES ({Object_Name_List[j]});"

    cursor.execute(cmd)

for i, j in enumerate(used_spectra):

    cmd = f"update {Database}.{Output_TableID} " \
        f"set Right_Ascension = '{RA[j]}'," \
        f"Declination = '{DEC[j]}'," \
        f"RA_HMS = '{RA_HMS[j]}'," \
        f"DEC_DMS = '{DEC_HMS[j]}'," \
        f"SDSS_ShortName = '{SDSS_Shortnames[j]}'," \
        f"DR16_ParentID = '{DR16_ParentID[j]}'," \
        f"DR16_Photometric_ID = '{DR16_PhotometricID[j]}'," \
        f"SDSS_DR16_Explore_Link ='{SDSS_Explore_Links[j]}'," \
        f"SDSS_DR16_Navigate_Link ='{SDSS_Navigate_Links[j]}'," \
        f"survey = '{Survey[j]}'," \
        f"run2d = '{Run2D[j]}'," \
        f"field = '{Field[j]}'," \
        f"spec_Plate = '{Plate[j]}'," \
        f"spec_MJD = '{MJD[j]}' ," \
        f"spec_FiberID = '{FiberID[j]}'," \
        f"SDSS_Type = '{SDSS_Type[j]}'," \
        f"spec_class_SDSS = '{Spec_Classification[j]}'," \
        f"spec_subclass_SDSS = '{Spec_Sub_Classification[j]}'," \
        f"spec_Human_Comments = '{Spec_Human_Comments[j]}'," \
        f"z_SDSS_spec = '{Spec_z[j]}'," \
        f"z_err_SDSS_spec = '{Spec_z_err[j]}'," \
        f"z_warning_SDSS_spec = '{Spec_z_warning[j]}'," \
        f"Distance_MPC = '{Distance[j]}'," \
        f"Max_Distance_MPC = '{MaxDistance[j]}'," \
        f"Min_Distance_MPC = '{MinDistance[j]}'," \
        f"Distance_Modulus='{Distance_Modulus[j]}'," \
        f"Distance_Modulus_Err='{Distance_Modulus_Err[j]}'," \
        f"median_SNR_SDSS_spec = '{Median_SNR[j]}'," \
        f"u_extinction = '{Extinction_u[j]}'," \
        f"g_extinction = '{Extinction_g[j]}'," \
        f"r_extinction = '{Extinction_r[j]}'," \
        f"i_extinction = '{Extinction_i[j]}'," \
        f"z_extinction = '{Extinction_z[j]}'," \
        f"generalised_extinction = '{Mean_E_BminusV[j]}'," \
        f"u_petro_mag = '{Petrosian_u[j]}'," \
        f"u_petro_mag_err = '{Petrosian_u_err[j]}'," \
        f"u_petro_AB_mag = '{Petrosian_u_AB[j]}'," \
        f"u_petro_AB_mag_err = '{Petrosian_u_AB_err[j]}'," \
        f"u_petro_rad= '{Petrosian_u_rad[j]}'," \
        f"u_petro_rad_err= '{Petrosian_u_rad_err[j]}'," \
        f"g_petro_mag = '{Petrosian_g[j]}'," \
        f"g_petro_mag_err = '{Petrosian_g_err[j]}'," \
        f"g_petro_AB_mag = '{Petrosian_g_AB[j]}'," \
        f"g_petro_AB_mag_err = '{Petrosian_g_AB_err[j]}'," \
        f"g_petro_rad= '{Petrosian_g_rad[j]}'," \
        f"g_petro_rad_err= '{Petrosian_g_rad_err[j]}'," \
        f"r_petro_mag = '{Petrosian_r[j]}'," \
        f"r_petro_mag_err = '{Petrosian_r_err[j]}'," \
        f"r_petro_AB_mag = '{Petrosian_r_AB[j]}'," \
        f"r_petro_AB_mag_err = '{Petrosian_r_AB_err[j]}'," \
        f"r_petro_rad= '{Petrosian_r_rad[j]}'," \
        f"r_petro_rad_err= '{Petrosian_r_rad_err[j]}'," \
        f"i_petro_mag = '{Petrosian_i[j]}'," \
        f"i_petro_mag_err = '{Petrosian_i_err[j]}'," \
        f"i_petro_AB_mag = '{Petrosian_i_AB[j]}'," \
        f"i_petro_AB_mag_err = '{Petrosian_i_AB_err[j]}'," \
        f"i_petro_rad= '{Petrosian_i_rad[j]}'," \
        f"i_petro_rad_err= '{Petrosian_i_rad_err[j]}'," \
        f"z_petro_mag = '{Petrosian_z[j]}'," \
        f"z_petro_mag_err = '{Petrosian_z_err[j]}'," \
        f"z_petro_AB_mag = '{Petrosian_z_AB[j]}'," \
        f"z_petro_AB_mag_err = '{Petrosian_z_AB_err[j]}'," \
        f"z_petro_rad= '{Petrosian_z_rad[j]}'," \
        f"z_petro_rad_err= '{Petrosian_z_rad_err[j]}'," \
        f"Petro_Total_Host_Mass_Estimate = '{TotalHostMassEstimate_Petro[j]}'," \
        f"Petro_Total_Host_Mass_Estimate_err = '{TotalHostMassEstimate_Err_Petro[j]}'," \
        f"Model_Total_Host_Mass_Estimate = '{TotalHostMassEstimate_Model[j]}'," \
        f"Model_Total_Host_Mass_Estimate_err = '{TotalHostMassEstimate_Err_Model[j]}'," \
        f"u_model_mag = '{Model_u[j]}'," \
        f"u_model_mag_err = '{Model_u_err[j]}'," \
        f"u_model_AB_mag = '{Model_u_AB[j]}'," \
        f"u_model_AB_mag_err = '{Model_u_AB_err[j]}'," \
        f"g_model_mag = '{Model_g[j]}'," \
        f"g_model_mag_err = '{Model_g_err[j]}'," \
        f"g_model_AB_mag = '{Model_g_AB[j]}'," \
        f"g_model_AB_mag_err = '{Model_g_AB_err[j]}'," \
        f"r_model_mag = '{Model_r[j]}'," \
        f"r_model_mag_err = '{Model_r_err[j]}'," \
        f"r_model_AB_mag = '{Model_r_AB[j]}'," \
        f"r_model_AB_mag_err = '{Model_r_AB_err[j]}'," \
        f"i_model_mag = '{Model_i[j]}'," \
        f"i_model_mag_err = '{Model_i_err[j]}'," \
        f"i_model_AB_mag = '{Model_i_AB[j]}'," \
        f"i_model_AB_mag_err = '{Model_i_AB_err[j]}'," \
        f"z_model_mag = '{Model_z[j]}'," \
        f"z_model_mag_err = '{Model_z_err[j]}'," \
        f"z_model_AB_mag = '{Model_z_AB[j]}'," \
        f"z_model_AB_mag_err = '{Model_z_AB_err[j]}'," \
        f"u_minus_r_petro = '{Petrosian_u_minus_r[j]}'," \
        f"u_minus_r_petro_observation_err = '{Petrosian_u_minus_r_observation_err[j]}'," \
        f"g_minus_r_petro = '{Petrosian_g_minus_r[j]}'," \
        f"g_minus_r_petro_observation_err = '{Petrosian_g_minus_r_observation_err[j]}'," \
        f"g_minus_i_petro = '{Petrosian_g_minus_i[j]}'," \
        f"g_minus_i_petro_observation_err = '{Petrosian_g_minus_i_observation_err[j]}'," \
        f"u_minus_r_model = '{Model_u_minus_r[j]}'," \
        f"u_minus_r_model_observation_err = '{Model_u_minus_r_observation_err[j]}'," \
        f"g_minus_r_model = '{Model_g_minus_r[j]}'," \
        f"g_minus_r_model_observation_err = '{Model_g_minus_r_observation_err[j]}'," \
        f"g_minus_i_model = '{Model_g_minus_i[j]}'," \
        f"g_minus_i_model_observation_err = '{Model_g_minus_i_observation_err[j]}'," \
        f"flags = '{Spec_Flags[j]}'," \
        f"clean = '{Clean_Flag[j]}'," \
        f" where DR16_Spectroscopic_ID = '{Object_Name_List[j]}' "

    cursor.execute(cmd)

Hirogen_Functions.commit(Data)
