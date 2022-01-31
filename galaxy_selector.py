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
from random import randrange

User = 'Joe'
Input_TableID = "SDSS_Galaxy_Spectra"
Output_TableID = 'SDSS_Fake_Spectra'

number_to_output = 5000

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

used_spectra = []

while len(used_spectra) < number_to_output:
    j = randrange(len(Candidate_Data))
    if j in used_spectra:
        continue
    used_spectra.append(j)

    cmd = f"INSERT INTO {Database}.{Output_TableID} SELECT * FROM {Database}.{Input_TableID} WHERE DR16_Spectroscopic_ID = {Object_Name_List[j]};"

    cursor.execute(cmd)

Hirogen_Functions.commit(Data)
