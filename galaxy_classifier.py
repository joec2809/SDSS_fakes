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

TableID = "SDSS_Galaxy_Spectra"

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
others = 0

psb_idxs = []
mod_bal_idxs = []
others_idxs = []

for i, object in enumerate(non_nan_EW):
    if non_nan_EW[i] < 3:
        if non_nan_lick[i] - non_nan_lick_err[i] > 4:
            psb += 1
            psb_idxs.append(i)
        elif non_nan_lick[i] - non_nan_lick_err[i] > 1.31:
            mod_bal += 1
            mod_bal_idxs.append(i)
        else:
            others += 1
            others_idxs.append(i)
    else:
        others +=1
        others_idxs.append(i)

psb_percent = psb/len(non_nan_EW)*100
mod_bal_percent = mod_bal/len(non_nan_EW)*100
others_percent = (others/len(non_nan_EW))*100

print(f"{psb_percent}% of the sample are post-starburst")
print(f"{mod_bal_percent}% of the sample are moderately balmer strong")
print(f"{others_percent}% of the sample are other types")

psb_EW = non_nan_EW[psb_idxs]
psb_lick = non_nan_lick[psb_idxs]

mod_bal_EW = non_nan_EW[mod_bal_idxs]
mod_bal_lick = non_nan_lick[mod_bal_idxs]

others_EW = non_nan_EW[others_idxs]
others_lick = non_nan_lick[others_idxs]


fig, ax = plt.subplots()
ax.scatter(non_nan_lick, non_nan_EW)
ax.vlines(1.31, -10, 3, 'k')
ax.vlines(4, -5, 3, 'k', 'dashed')
ax.hlines(3, 1.31, 12, 'k')

plt.show()