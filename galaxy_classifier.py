import sys
import time
import warnings
import math
import os

import Hirogen_Functions
import mods_functions

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.constants as constants

from astropy.cosmology import Planck18
from astropy import units as u

User = 'Joe'
User_Config = Hirogen_Functions.user_config(user=User)

Database = User_Config[0]
Database_User = User_Config[1]
Database_Password = User_Config[2]
Main_Spectra_Path = User_Config[3]

TableID = "SDSS_FeVII_Fake_Spectra"

Data = Hirogen_Functions.database_connection(user=Database_User, password=Database_Password, database=Database)

print("\n")
cursor = Data.cursor()

cursor.execute(
    f"SELECT DR16_Spectroscopic_ID, lin_con_pEQW_Halpha, Lick_HDelta_Index, Lick_HDelta_Index_Err " 
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

if len(Candidate_Data) >= 1:

    Object_Name_List = [item[0] for item in Candidate_Data]
    h_alpha_EW = [item[1] for item in Candidate_Data]
    h_delta_lick = [item[2] for item in Candidate_Data]
    h_delta_lick_err = [item[3] for item in Candidate_Data]

h_alpha_EW = -np.array(h_alpha_EW)
h_delta_lick = np.array(h_delta_lick)
h_delta_lick_err = np.array(h_delta_lick_err)

non_nans = (h_alpha_EW > -1000) & (h_delta_lick != -999) & (h_delta_lick_err > 0) & (h_alpha_EW < 999)

non_nan_EW = h_alpha_EW[non_nans]
non_nan_lick = h_delta_lick[non_nans]
non_nan_lick_err = h_delta_lick_err[non_nans]

psb = 0
mod_bal = 0
quiscent = 0
sf = 0

star_forming_EW = non_nan_EW[non_nan_EW >= 3]
star_forming_lick = non_nan_lick[non_nan_EW >= 3]

quiscent_EW = non_nan_EW[non_nan_EW < 3]
quiscent_lick = non_nan_lick[non_nan_EW < 3]

for i, object in enumerate(non_nan_EW):
    if non_nan_EW[i] < 3:
        if non_nan_lick[i] - non_nan_lick_err[i] > 4:
            psb += 1
        elif non_nan_lick[i] - non_nan_lick_err[i] > 1.31:
            mod_bal += 1
        else:
            quiscent += 1
    else:
        sf +=1

psb_percent = psb/len(non_nan_EW)*100
mod_bal_percent = mod_bal/len(non_nan_EW)*100
quiscent_percent = (quiscent/len(non_nan_EW))*100
sf_percent = (sf/len(non_nan_EW))*100

print(f"{psb_percent}% of the sample are post-starburst, that's {psb} galaxies")
print(f"{mod_bal_percent}% of the sample are moderately balmer strong, that's {mod_bal} galaxies")
print(f"{quiscent_percent}% of the sample are other quiscent, that's {quiscent} galaxies")#
print(f"{sf_percent}% of the sample are star forming, that's {sf} galaxies")


fig, (ax1, ax2) = plt.subplots(2,1, sharex = True)

ax1.scatter(star_forming_lick, star_forming_EW, color = 'b', marker = '.')
ax1.set_ylim(5, 100)


ax2.scatter(star_forming_lick, star_forming_EW, color = 'b', marker = '.')
ax2.scatter(quiscent_lick, quiscent_EW, color = 'r', marker = '.')
ax2.vlines(1.31, -10, 3, 'k')
ax2.vlines(4, -5, 3, 'k', 'dashed')
ax2.hlines(3, 1.31, 12, 'k')

ax2.set_xlim(-6,10)
ax2.set_ylim(-2.5,4)



plt.show()