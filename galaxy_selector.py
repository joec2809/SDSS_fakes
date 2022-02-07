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

number_to_output = 10000

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
