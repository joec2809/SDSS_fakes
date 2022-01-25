import sys
import time
import warnings
import math
import os
import Hirogen_Functions
import mods_functions

import numpy as np
import scipy.constants as constants
import matplotlib.pyplot as plt
import pandas as pd

from astropy import units as u
from astropy.visualization import hist
from datetime import datetime, date

User = 'Joe'
User_Config = Hirogen_Functions.user_config(user=User)

TableID = "SDSS_Galaxy_Spectra"


# Search for 'QUERY' to find the actual database access location where parameters can be adjusted

Database = User_Config[0]
Database_User = User_Config[1]
Database_Password = User_Config[2]
Main_Spectra_Path = User_Config[3]

Data = Hirogen_Functions.database_connection(user=Database_User, password=Database_Password, database=Database)

print("\n")
cursor = Data.cursor()

cursor.execute(
            f"SELECT DR16_Spectroscopic_ID, lin_con_pEQW_FeVII6008, lin_con_pEQW_FeX6376, lin_con_pEQW_FeXI7894, lin_con_pEQW_FeXIV5304, "
            f"Strong_FeVII_Flag, Strong_FeX_Flag, Strong_FeXI_Flag, Strong_FeXIV_Flag, ECLE_Candidate_Score "
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
        fevii_pEQW = [item[1] for item in Candidate_Data]
        fex_pEQW = [item[2] for item in Candidate_Data]
        fexi_pEQW = [item[3] for item in Candidate_Data]
        fexiv_pEQW = [item[4] for item in Candidate_Data]
        fevii_flag = [item[5] for item in Candidate_Data]
        fex_flag = [item[6] for item in Candidate_Data]
        fexi_flag = [item[7] for item in Candidate_Data]
        fexiv_flag = [item[8] for item in Candidate_Data]
        ecle_candidate_score = [item[9] for item in Candidate_Data]
        

else:
    print("Candidate_Data Length error: Check and try again.")

    sys.exit()


ecle_candidate = np.zeros(len(ecle_candidate_score))
for i, score in enumerate(ecle_candidate_score):
    if score >= 7 or (fevii_flag[i]+fex_flag[i]+fexi_flag[i]+fexiv_flag[i])> 0:
        ecle_candidate[i] = 1

# FeVII
fevii_pEQW = np.array(fevii_pEQW)
fevii_scores = np.array(ecle_candidate)

nans = np.argwhere(fevii_pEQW <= -999)
new_fevii_pEQW = np.delete(fevii_pEQW, nans)
new_fevii_scores = np.delete(fevii_scores, nans)

fevii_sort_ind = new_fevii_pEQW.argsort()
fevii_pEQW_sort = fevii_pEQW[fevii_sort_ind]
fevii_scores_sort = fevii_scores[fevii_sort_ind]

fevii_bins_template = mods_functions.create_bins(fevii_pEQW_sort, 1)
fevii_bins = np.array(pd.cut(fevii_pEQW_sort, bins = fevii_bins_template, right = True).value_counts())
fevii_binned = np.array(np.split(fevii_scores_sort, fevii_bins), dtype = object)

print(fevii_binned)

fevii_scores_binned = np.zeros(len(fevii_binned))

for i, sub in enumerate(fevii_binned):
    fevii_scores_binned[i] = sum(sub)

print(fevii_bins)
print(fevii_scores_binned)

"""# FeX
fex_scores = np.zeros((len(fex_pEQW), 2))
fex_scores[:,0] = fex_pEQW
fex_scores[:,1] = ecle_candidate

nans = np.argwhere(fex_scores[:,0] <= -999)
new_fex_scores = np.delete(fex_scores, nans, axis = 0)

fex_sort_ind = new_fex_scores[:,0].argsort()
fex_scores_sort = fex_scores[fex_sort_ind]


# FeXI
fexi_scores = np.zeros((len(fexi_pEQW), 2))
fexi_scores[:,0] = fexi_pEQW
fexi_scores[:,1] = ecle_candidate

nans = np.argwhere(fexi_scores[:,0] <= -999)
new_fexi_scores = np.delete(fexi_scores, nans, axis = 0)

fexi_sort_ind = new_fexi_scores[:,0].argsort()
fexi_scores_sort = new_fexi_scores[fexi_sort_ind]


# FeXIV
fexiv_scores = np.zeros((len(fexiv_pEQW), 2))
fexiv_scores[:,0] = fexiv_pEQW
fexiv_scores[:,1] = ecle_candidate

nans = np.argwhere(fexiv_scores[:,0] <= -999)
new_fexiv_scores = np.delete(fexiv_scores, nans, axis = 0)

fexiv_sort_ind = new_fexiv_scores[:,0].argsort()
fexiv_scores_sort = fexiv_scores[fexiv_sort_ind]
"""