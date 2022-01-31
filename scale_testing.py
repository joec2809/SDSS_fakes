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
from astropy.convolution import convolve, Box1DKernel
from random import uniform

User = 'Joe'
User_Config = Hirogen_Functions.user_config(user=User)
Main_Config = Hirogen_Functions.main_config()

TableID = "SDSS_Fake_Spectra"

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
        f"survey, run2d, Standard_Inclusion,Path_Override_Flag, Path_Override, Follow_Up_ID, Smoothing_Override, "
        f"z_corrected_flag, extinction_corrected_flag, lin_con_pEQW_OIII5007, lin_con_LineFlux_OIII5007 "
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
    OIII_pEQW = [item[15] for item in galaxy_data]
    OIII_flux = [item[16] for item in galaxy_data]

else:
    print("galaxy_data Length error: Check and try again.")

    sys.exit()

galaxy_FilePaths = Hirogen_Functions.sdss_peaks_spectra_file_path_generator(
    Main_Spectra_Path, galaxy_Plate_List, galaxy_MJD_List, galaxy_FiberID_List, galaxy_Survey_List, galaxy_run2d_List, 
    galaxy_Path_Override_Flag_List, galaxy_Path_Override_List, Mode = 'fakes'
)

file = np.random.randint(0,len(galaxy_FilePaths))

spectra = Hirogen_Functions.sdss_spectrum_reader(galaxy_FilePaths[file], galaxy_Redshift_List[file], galaxy_Extinction_List[file]) 

lamb_rest = spectra[0]
flux = spectra[1]

lines = Hirogen_Functions.lines_for_analysis()

primary = Hirogen_Functions.primary_lines()
labels = Hirogen_Functions.primary_line_labels()

fig, ax = plt.subplots(figsize = (20,10))

for i, line in enumerate(lines):
    if line in primary:
        ax.axvline(lines[line][0], 0, 1, color = 'k')
        ax.annotate(line, (lines[line][0], 30))

ax.plot(lamb_rest, flux)

if OIII_pEQW[file] < Main_Config[6]:
    print("yes")


plt.show()