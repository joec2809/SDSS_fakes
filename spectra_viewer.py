import sys
import time
import warnings
import math
import os
import Hirogen_Functions
import mods_functions

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy import units as u

User = 'Joe'
User_Config = Hirogen_Functions.user_config(user=User)

TableID = "SDSS_FeVII_Fake_Spectra"

# Search for 'QUERY' to find the actual database access location where parameters can be adjusted

Database = User_Config[0]
Database_User = User_Config[1]
Database_Password = User_Config[2]
Main_Spectra_Path = User_Config[3]

Data = Hirogen_Functions.database_connection(user=Database_User, password=Database_Password, database=Database)


print("\n")
cursor = Data.cursor()

cursor.execute(
    f"SELECT DR16_Spectroscopic_ID "
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

ID_Data = cursor.fetchall()

if len(ID_Data) >= 1:
    ID_List = [item[0] for item in ID_Data]

rand_idx = np.random.randint(0,len(ID_List))

IDs = (ID_List[rand_idx])

#IDs = (494429519780800512)

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
    f"WHERE DR16_Spectroscopic_ID = {IDs}"
)

Candidate_Data = cursor.fetchall()

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

filepath = Hirogen_Functions.sdss_peaks_spectra_file_path_generator(Main_Spectra_Path, Plate_List, MJD_List, FiberID_List, Survey_List, run2d_List, Path_Override_Flag_List, Path_Override_List, Mode = 'fakes')

spectra = Hirogen_Functions.sdss_spectrum_reader(filepath[0], Redshift_List[0], Extinction_List[0])

fake_wave, fake_flux = spectra[0], spectra[1]

filepath = Hirogen_Functions.sdss_spectra_file_path_generator(Main_Spectra_Path, Plate_List, MJD_List, FiberID_List, Survey_List, run2d_List, Path_Override_Flag_List, Path_Override_List)

print(filepath)

spectra = Hirogen_Functions.sdss_spectrum_reader(filepath[0], Redshift_List[0], Extinction_List[0])

wave, flux = spectra[0], spectra[1]

fig, ax = plt.subplots()

ax.plot(wave, flux, 'r')

ax.plot(fake_wave, fake_flux, 'k')


plt.show()