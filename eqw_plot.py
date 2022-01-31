import sys
import time
from timeit import default_timer
import warnings
import math
import os
import Hirogen_Functions
import mods_functions

import numpy as np
import scipy.constants as constants
import matplotlib.pyplot as plt
import pandas as pd

from scipy.interpolate import CubicSpline, interp1d
from astropy import units as u
from astropy.visualization import hist
from datetime import datetime, date

User = 'Joe'
User_Config = Hirogen_Functions.user_config(user=User)

Galaxy_TableID = "SDSS_Galaxy_Spectra"
Fakes_TableID = "SDSS_Fake_Spectra"

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
    f"FROM `{Database}`.`{Fakes_TableID}`"
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

    fake_Object_Name_List = [item[0] for item in Candidate_Data]
    fake_fevii_pEQW = [item[1] for item in Candidate_Data]
    fake_fex_pEQW = [item[2] for item in Candidate_Data]
    fake_fexi_pEQW = [item[3] for item in Candidate_Data]
    fake_fexiv_pEQW = [item[4] for item in Candidate_Data]
    fake_fevii_flag = [item[5] for item in Candidate_Data]
    fake_fex_flag = [item[6] for item in Candidate_Data]
    fake_fexi_flag = [item[7] for item in Candidate_Data]
    fake_fexiv_flag = [item[8] for item in Candidate_Data]
    fake_ecle_candidate_score = [item[9] for item in Candidate_Data]
        

else:
    print("Candidate_Data Length error: Check and try again.")

    sys.exit()

fake_ecle_candidate = np.zeros(len(fake_ecle_candidate_score))
for i, score in enumerate(fake_ecle_candidate_score):
    if score >= 7 or (fake_fevii_flag[i]+fake_fex_flag[i]+fake_fexi_flag[i]+fake_fexiv_flag[i])> 0:
        fake_ecle_candidate[i] = 1


cursor = Data.cursor()

cursor.execute(
    f"SELECT lin_con_pEQW_FeVII6008, lin_con_pEQW_FeX6376, lin_con_pEQW_FeXI7894, lin_con_pEQW_FeXIV5304, "
    f"Strong_FeVII_Flag, Strong_FeX_Flag, Strong_FeXI_Flag, Strong_FeXIV_Flag, ECLE_Candidate_Score "
    f"FROM `{Database}`.`{Galaxy_TableID}`"
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

    galaxy_fevii_pEQW = [item[0] for item in Candidate_Data]
    galaxy_fex_pEQW = [item[1] for item in Candidate_Data]
    galaxy_fexi_pEQW = [item[2] for item in Candidate_Data]
    galaxy_fexiv_pEQW = [item[3] for item in Candidate_Data]
    galaxy_fevii_flag = [item[4] for item in Candidate_Data]
    galaxy_fex_flag = [item[5] for item in Candidate_Data]
    galaxy_fexi_flag = [item[6] for item in Candidate_Data]
    galaxy_fexiv_flag = [item[7] for item in Candidate_Data]
    galaxy_ecle_candidate_score = [item[8] for item in Candidate_Data]

# Determine number of ECLE candidates

galaxy_ecle_candidate = np.zeros(len(galaxy_ecle_candidate_score))
for i, score in enumerate(galaxy_ecle_candidate_score):
    if score >= 7 or (galaxy_fevii_flag[i]+galaxy_fex_flag[i]+galaxy_fexi_flag[i]+galaxy_fexiv_flag[i])> 0:
        galaxy_ecle_candidate[i] = 1

fevii_pEQW = np.array(fake_fevii_pEQW)
nans = np.argwhere(fevii_pEQW <= -999)
fevii_pEQW[nans] = 0

fex_pEQW = np.array(fake_fex_pEQW)
nans = np.argwhere(fex_pEQW <= -999)
fex_pEQW[nans] = 0

fexi_pEQW = np.array(fake_fexi_pEQW)
nans = np.argwhere(fexi_pEQW <= -999)
fexi_pEQW[nans] = 0

fexiv_pEQW = np.array(fake_fexiv_pEQW)
nans = np.argwhere(fexiv_pEQW <= -999)
fexiv_pEQW[nans] = 0

total_pEQW = fevii_pEQW + fex_pEQW + fexi_pEQW + fexiv_pEQW

total_sort_ind = np.argsort(total_pEQW)
total_pEQW_sort = total_pEQW[total_sort_ind]
scores_sort = fake_ecle_candidate[total_sort_ind]

total_info = pd.cut(total_pEQW_sort, 20, retbins=True)

total_cuts = np.cumsum(total_info[0].value_counts())[:-1]

total_bins = total_info[1]

scores_binned = np.split(scores_sort, total_cuts)

total_scores = np.zeros(len(scores_binned))

for i, array in enumerate(scores_binned):
    total_scores[i] = sum(array)

bin_sizes = np.zeros(len(scores_binned))
for i, array in enumerate(scores_binned):
    bin_sizes[i] = len(array)

detection_efficiency = total_scores/bin_sizes

bin_centres = np.zeros(len(detection_efficiency))
for i in range(len(bin_centres)):
    bin_centres[i] = (total_bins[i] + total_bins[i+1])/2


#Errors

fevii_pEQW = np.array(galaxy_fevii_pEQW)
nans = np.argwhere(fevii_pEQW <= -999)
fevii_pEQW[nans] = 0

fex_pEQW = np.array(fake_fex_pEQW)
nans = np.argwhere(fex_pEQW <= -999)
fex_pEQW[nans] = 0

fexi_pEQW = np.array(fake_fexi_pEQW)
nans = np.argwhere(fexi_pEQW <= -999)
fexi_pEQW[nans] = 0

fexiv_pEQW = np.array(fake_fexiv_pEQW)
nans = np.argwhere(fexiv_pEQW <= -999)
fexiv_pEQW[nans] = 0

galaxy_pEQW = fevii_pEQW + fex_pEQW + fexi_pEQW + fexiv_pEQW

galaxy_sort_ind = np.argsort(galaxy_pEQW)
galaxy_pEQW_sort = galaxy_pEQW[galaxy_sort_ind]
galaxy_scores_sort = galaxy_ecle_candidate[galaxy_sort_ind]

bins_cut= pd.cut(galaxy_pEQW_sort, total_bins)

error_cuts = np.cumsum(bins_cut.value_counts())[:-1]

pEQW_for_error = np.split(galaxy_pEQW_sort, error_cuts)

error = np.zeros(len(pEQW_for_error))

for i, arr in enumerate(pEQW_for_error):
    error[i] = np.std(arr)





















not_nans = ~np.isnan(detection_efficiency)
fin_detection_eff = detection_efficiency[not_nans]
fin_bin_centres = bin_centres[not_nans]

cs = CubicSpline(fin_bin_centres, fin_detection_eff)
one_d = interp1d(fin_bin_centres, fin_detection_eff)
xs = np.arange(bin_centres[0], bin_centres[-1], 1)

fig, ax = plt.subplots()
ax.plot(xs, cs(xs), marker = ',')
ax.plot(xs, one_d(xs), marker = ',')
ax.errorbar(bin_centres, detection_efficiency, xerr = error, fmt = 'k.', ls = 'none', capsize = 5)

ax.set_xlabel(r'Equivalent Width, $\AA$')
ax.set_ylabel('Detection Efficiency')
ax.set_ylim(0,1.1)
ax.tick_params(top = True, right = True, direction = 'in')

plt.savefig("example_detection_efficiency_plot.png")

plt.show()