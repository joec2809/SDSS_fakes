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

User = 'Joe'
Input_TableID = "SDSS_Galaxy_Spectra"
Output_TableID = 'SDSS_Non_FeVII_Fake_Spectra'

number_to_output = 10000

TDE_type_dist = True

if TDE_type_dist:
    no_psb = (5/41)*number_to_output
    no_qbs = (8/41)*number_to_output
    no_quiscent = (13/41)*number_to_output
    no_sf = (15/41)*number_to_output

    gal_types = np.array([no_psb, no_qbs, no_quiscent, no_sf])
    gal_types_rounded = np.zeros(4)

    no_floored = 0
    no_ceiled = 0
    i = 0

    while no_floored < 2 or no_ceiled < 2:
        randint = np.random.randint(0,2)
        if randint == 0:
            if no_floored != 2:
                gal_types_rounded[i] = np.floor(gal_types[i])
                no_floored += 1
            else:
                gal_types_rounded[i] = np.ceil(gal_types[i])
                no_ceiled += 1
        elif randint == 1:
            if no_ceiled != 2:
                gal_types_rounded[i] = np.ceil(gal_types[i])
                no_ceiled += 1
            else:
                gal_types_rounded[i] = np.floor(gal_types[i])
                no_floored += 1
        i += 1

    psb_cond = "WHERE (lin_con_pEQW_Halpha > -3) AND (Lick_HDelta_Index - Lick_HDelta_Index_Err > 4)"
    qbs_cond = "WHERE (lin_con_pEQW_Halpha > -3) AND (Lick_HDelta_Index - Lick_HDelta_Index_Err BETWEEN 1.31 AND 4)"
    quiscent_cond = "WHERE (lin_con_pEQW_Halpha > -3) AND (Lick_HDelta_Index - Lick_HDelta_Index_Err < 1.31)"
    sf_cond = "WHERE (lin_con_pEQW_Halpha < -3)"

    conditions = [psb_cond, qbs_cond, quiscent_cond, sf_cond]

    User_Config = Hirogen_Functions.user_config(user=User)

    Database = User_Config[0]
    Database_User = User_Config[1]
    Database_Password = User_Config[2]

    Data = Hirogen_Functions.database_connection(user=Database_User, password=Database_Password, database=Database)

    for k, condition in enumerate(conditions):

        print("\n")
        cursor = Data.cursor()

        cursor.execute(
                f"SELECT DR16_Spectroscopic_ID "
                f"FROM `{Database}`.`{Input_TableID}`"
                #f"WHERE Manually_Inspected_Flag != -10 AND lin_con_LineFlux_Halpha is NULL"
                f"{condition}"
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

        #print(len(Candidate_Data))

        if len(Candidate_Data) >= 1:

            Object_Name_List = [item[0] for item in Candidate_Data]

        used_spectra = []
        while len(used_spectra) < gal_types_rounded[k]:
            j = np.random.randint(len(Candidate_Data))
            if j in used_spectra:
                continue
            used_spectra.append(j)

            cmd = f"INSERT INTO {Database}.{Output_TableID} SELECT * FROM {Database}.{Input_TableID} WHERE DR16_Spectroscopic_ID = {Object_Name_List[j]};"

            cursor.execute(cmd)

        Hirogen_Functions.commit(Data)

else:
    User_Config = Hirogen_Functions.user_config(user=User)

    Database = User_Config[0]
    Database_User = User_Config[1]
    Database_Password = User_Config[2]

    Data = Hirogen_Functions.database_connection(user=Database_User, password=Database_Password, database=Database)

    print("\n")
    cursor = Data.cursor()

    cursor.execute(
            f"SELECT DR16_Spectroscopic_ID "
            f"FROM `{Database}`.`{Input_TableID}`"
            #f"WHERE Manually_Inspected_Flag != -10 AND lin_con_LineFlux_Halpha is NULL"
            f"WHERE lin_con_LineFlux_Halpha is not NULL"
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

    #print(len(Candidate_Data))

    if len(Candidate_Data) >= 1:

        Object_Name_List = [item[0] for item in Candidate_Data]

    used_spectra = []
    while len(used_spectra) < number_to_output:
        j = np.random.randint(len(Candidate_Data))
        if j in used_spectra:
            continue
        used_spectra.append(j)

        cmd = f"INSERT INTO {Database}.{Output_TableID} SELECT * FROM {Database}.{Input_TableID} WHERE DR16_Spectroscopic_ID = {Object_Name_List[j]};"

        cursor.execute(cmd)

    Hirogen_Functions.commit(Data)