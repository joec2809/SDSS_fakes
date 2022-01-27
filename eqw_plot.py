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

# Determine number of ECLE candidates

ecle_candidate = np.zeros(len(ecle_candidate_score))
for i, score in enumerate(ecle_candidate_score):
    if score >= 7 or (fevii_flag[i]+fex_flag[i]+fexi_flag[i]+fexiv_flag[i])> 0:
        ecle_candidate[i] = 1

# FeVII
fevii_pEQW = np.array(fevii_pEQW)
fevii_scores = np.array(ecle_candidate)

# Remove nan values
nans = np.argwhere(fevii_pEQW <= -999)
new_fevii_pEQW = np.delete(fevii_pEQW, nans)
new_fevii_scores = np.delete(fevii_scores, nans)

# Sort arrays by pEQW
fevii_sort_ind = new_fevii_pEQW.argsort()
fevii_pEQW_sort = fevii_pEQW[fevii_sort_ind]
fevii_scores_sort = fevii_scores[fevii_sort_ind]

# Sort pEQWs and scores into bins
fevii_info = pd.cut(fevii_pEQW_sort, 20, retbins=True)

fevii_cuts = np.cumsum(fevii_info[0].value_counts())[:-1]

fevii_bins = fevii_info[1]

fevii_scores_binned = np.split(fevii_scores_sort, fevii_cuts)

# Find total ECLE score per bin and determine efficiency by which they have been detected
total_fevii_scores = np.zeros(len(fevii_scores_binned))

for i, array in enumerate(fevii_scores_binned):
    total_fevii_scores[i] = sum(array)

bin_sizes = np.zeros(len(fevii_scores_binned))
for i, array in enumerate(fevii_scores_binned):
    bin_sizes[i] = len(array)

fevii_detection_efficiency = total_fevii_scores/bin_sizes

fevii_bin_centres = np.zeros(len(fevii_detection_efficiency))
for i in range(len(fevii_bin_centres)):
    fevii_bin_centres[i] = (fevii_bins[i] + fevii_bins[i+1])/2


# FeX, same process as above
fex_pEQW = np.array(fex_pEQW)
fex_scores = np.array(ecle_candidate)

nans = np.argwhere(fex_pEQW <= -999)
new_fex_pEQW = np.delete(fex_pEQW, nans)
new_fex_scores = np.delete(fex_scores, nans)

fex_sort_ind = new_fex_pEQW.argsort()
fex_pEQW_sort = fex_pEQW[fex_sort_ind]
fex_scores_sort = fex_scores[fex_sort_ind]

fex_info = pd.cut(fex_pEQW_sort, 20, retbins=True)

fex_cuts = np.cumsum(fex_info[0].value_counts())[:-1]

fex_bins = fex_info[1]

fex_scores_binned = np.split(fex_scores_sort, fex_cuts)

total_fex_scores = np.zeros(len(fex_scores_binned))

for i, array in enumerate(fex_scores_binned):
    total_fex_scores[i] = sum(array)

bin_sizes = np.zeros(len(fex_scores_binned))
for i, array in enumerate(fex_scores_binned):
    bin_sizes[i] = len(array)

fex_detection_efficiency = total_fex_scores/bin_sizes

fex_bin_centres = np.zeros(len(fex_detection_efficiency))
for i in range(len(fex_bin_centres)):
    fex_bin_centres[i] = (fex_bins[i] + fex_bins[i+1])/2


# FeXI, same process as above
fexi_pEQW = np.array(fexi_pEQW)
fexi_scores = np.array(ecle_candidate)

nans = np.argwhere(fexi_pEQW <= -999)
new_fexi_pEQW = np.delete(fexi_pEQW, nans)
new_fexi_scores = np.delete(fexi_scores, nans)

fexi_sort_ind = new_fexi_pEQW.argsort()
fexi_pEQW_sort = fexi_pEQW[fexi_sort_ind]
fexi_scores_sort = fexi_scores[fexi_sort_ind]

fexi_info = pd.cut(fexi_pEQW_sort, 20, retbins=True)

fexi_cuts = np.cumsum(fexi_info[0].value_counts())[:-1]

fexi_bins = fexi_info[1]

fexi_scores_binned = np.split(fexi_scores_sort, fexi_cuts)

total_fexi_scores = np.zeros(len(fexi_scores_binned))

for i, array in enumerate(fexi_scores_binned):
    total_fexi_scores[i] = sum(array)

bin_sizes = np.zeros(len(fexi_scores_binned))
for i, array in enumerate(fexi_scores_binned):
    bin_sizes[i] = len(array)

fexi_detection_efficiency = total_fexi_scores/bin_sizes

fexi_bin_centres = np.zeros(len(fexi_detection_efficiency))
for i in range(len(fexi_bin_centres)):
    fexi_bin_centres[i] = (fexi_bins[i] + fexi_bins[i+1])/2


# FeXIV, same process as above
fexiv_pEQW = np.array(fexiv_pEQW)
fexiv_scores = np.array(ecle_candidate)

nans = np.argwhere(fexiv_pEQW <= -999)
new_fexiv_pEQW = np.delete(fexiv_pEQW, nans)
new_fexiv_scores = np.delete(fexiv_scores, nans)

fexiv_sort_ind = np.argsort(new_fexiv_pEQW)
fexiv_pEQW_sort = fexiv_pEQW[fexiv_sort_ind]
fexiv_scores_sort = fexiv_scores[fexiv_sort_ind]

fexiv_info = pd.cut(fexiv_pEQW_sort, 20, retbins=True)

fexiv_cuts = np.cumsum(fexiv_info[0].value_counts())[:-1]

fexiv_bins = fexiv_info[1]

fexiv_scores_binned = np.split(fexiv_scores_sort, fexiv_cuts)

total_fexiv_scores = np.zeros(len(fexiv_scores_binned))

for i, array in enumerate(fexiv_scores_binned):
    total_fexiv_scores[i] = sum(array)

bin_sizes = np.zeros(len(fexiv_scores_binned))
for i, array in enumerate(fexiv_scores_binned):
    bin_sizes[i] = len(array)

fexiv_detection_efficiency = total_fexiv_scores/bin_sizes

fexiv_bin_centres = np.zeros(len(fexiv_detection_efficiency))
for i in range(len(fexiv_bin_centres)):
    fexiv_bin_centres[i] = (fexiv_bins[i] + fexiv_bins[i+1])/2

print(fevii_bin_centres)
print(f"fevii {fevii_detection_efficiency}")

print(fex_bin_centres)
print(f"fex {fex_detection_efficiency}")

print(fexi_bin_centres)
print(f"fexi {fexi_detection_efficiency}")

print(fexiv_bin_centres)
print(f"fexiv {fexiv_detection_efficiency}")


# Plot pEQw vs detection efficiency

fig, axs = plt.subplots(2,2, sharey=True, figsize = (10,10))

axs[0,0].plot(fevii_bin_centres, fevii_detection_efficiency)
axs[0,0].set_title("FeVII Lines")

axs[0,1].plot(fex_bin_centres, fex_detection_efficiency)
axs[0,1].set_title("FeX Lines")

axs[1,0].plot(fexi_bin_centres, fexi_detection_efficiency)
axs[1,0].set_title("FeXI Lines")

axs[1,1].plot(fexiv_bin_centres, fexiv_detection_efficiency)
axs[1,1].set_title("FeXIV Lines")

plt.setp(axs[0, :], xlabel = r'Equivalent Width, $\AA$')
plt.setp(axs[-1, :], xlabel = r'Equivalent Width, $\AA$')
plt.setp(axs[:, 0], ylabel = 'Detection Efficiency')

plt.show()