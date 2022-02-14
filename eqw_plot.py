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

from scipy.interpolate import CubicSpline, interp1d
from scipy.optimize import curve_fit

User = 'Joe'
User_Config = Hirogen_Functions.user_config(user=User)

Galaxy_TableID = "SDSS_Galaxy_Spectra"
Fakes_TableID = "SDSS_Non_FeVII_Fake_Spectra"

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


# Non-FeVII
total_pEQW = (fex_pEQW + fexi_pEQW + fexiv_pEQW)/3

#FeVII
#total_pEQW = (fevii_pEQW + fex_pEQW + fexi_pEQW + fexiv_pEQW)/4

total_sort_ind = np.argsort(total_pEQW)
total_pEQW_sort = total_pEQW[total_sort_ind]
scores_sort = fake_ecle_candidate[total_sort_ind]

# Main Plot

bin_width = 5

total_bins = mods_functions.create_bins(total_pEQW_sort[0], total_pEQW_sort[-1], bin_width)

total_info = pd.cut(total_pEQW_sort, total_bins)

total_cuts = np.cumsum(total_info.value_counts())[:-1]

scores_binned = np.split(scores_sort, total_cuts)

total_scores = np.zeros(len(scores_binned))

for i, array in enumerate(scores_binned):
    total_scores[i] = sum(array)

bin_sizes = np.zeros(len(scores_binned))
for i, array in enumerate(scores_binned):
    bin_sizes[i] = len(array)

det_eff = (total_scores)/(bin_sizes)

det_error = np.sqrt(((total_scores+1)*(total_scores+2))/((bin_sizes+2)*(bin_sizes+3)) - ((total_scores+1)**2)/((bin_sizes+2)**2))

det_err = np.zeros((2, len(det_error)))

for i, error in enumerate(det_error):
    det_err[0][i] = error
    det_err[1][i] = error

    if det_eff[i] - det_err[0][i] < 0:
        det_err[0][i] = det_eff[i]

    if det_eff[i] + det_err[1][i] > 1:
        det_err[1][i] = 1 - det_eff[i]

pEQW_err = np.sqrt(bin_width)

bin_centres = np.zeros(len(det_eff))
for i in range(len(bin_centres)):
    bin_centres[i] = (total_bins[i] + total_bins[i+1])/2

not_nans = ~np.isnan(det_eff)
fin_det_eff = det_eff[not_nans]
fin_bin_centres = bin_centres[not_nans]

# Inset

zoom_bin_width = 2

zoom_bins = mods_functions.create_bins(total_pEQW_sort[mods_functions.find_nearest(total_pEQW_sort, -60)], total_pEQW_sort[-1], zoom_bin_width)

zoom_info = pd.cut(total_pEQW_sort, zoom_bins)

zoom_cuts = np.cumsum(zoom_info.value_counts())[:-1]

zoom_scores_binned = np.split(scores_sort, zoom_cuts)

zoom_scores = np.zeros(len(zoom_scores_binned))

for i, array in enumerate(zoom_scores_binned):
    zoom_scores[i] = sum(array)

zoom_bin_sizes = np.zeros(len(zoom_scores_binned))
for i, array in enumerate(zoom_scores_binned):
    zoom_bin_sizes[i] = len(array)

zoom_det_eff = zoom_scores/zoom_bin_sizes

zoom_det_error = np.sqrt(((zoom_scores+1)*(zoom_scores+2))/((zoom_bin_sizes+2)*(zoom_bin_sizes+3)) - ((zoom_scores+1)**2)/((zoom_bin_sizes+2)**2))

zoom_det_err = np.zeros((2, len(zoom_det_error)))

for i, error in enumerate(zoom_det_error):
    zoom_det_err[0][i] = error
    zoom_det_err[1][i] = error

    if zoom_det_eff[i] - zoom_det_err[0][i] < 0:
        zoom_det_err[0][i] = zoom_det_eff[i]

    if zoom_det_eff[i] + zoom_det_err[1][i] > 1:
        zoom_det_err[1][i] = 1 - zoom_det_eff[i]

zoom_pEQW_err = np.sqrt(zoom_bin_width)

zoom_bin_centres = np.zeros(len(zoom_det_eff))
for i in range(len(zoom_bin_centres)):
    zoom_bin_centres[i] = (zoom_bins[i] + zoom_bins[i+1])/2

zoom_not_nans = ~np.isnan(zoom_det_eff)
zoom_fin_det_eff = zoom_det_eff[zoom_not_nans]
zoom_fin_bin_centres = zoom_bin_centres[zoom_not_nans]

# Curves

def func(x, w, s):
    return 1/(1 + np.exp((x-w)/s))

par_1, cov_1 = curve_fit(func, fin_bin_centres[fin_det_eff >= 0.5], fin_det_eff[fin_det_eff >= 0.5])
par_2, cov_2 = curve_fit(func, fin_bin_centres[fin_det_eff < 0.5], fin_det_eff[fin_det_eff < 0.5])

zoom_par_1, zoom_cov_1 = curve_fit(func, zoom_fin_bin_centres[zoom_fin_det_eff >= 0.5], zoom_fin_det_eff[zoom_fin_det_eff >= 0.5])
zoom_par_2, zoom_cov_2 = curve_fit(func, zoom_fin_bin_centres[zoom_fin_det_eff < 0.5], zoom_fin_det_eff[zoom_fin_det_eff < 0.5])

w_half = (par_1[0]+par_2[0])/2
zoom_w_half = (zoom_par_1[0]+zoom_par_2[0])/2

xs_1 = np.arange(-470, w_half, 0.1)
xs_2 = np.arange(w_half, 100, 0.1)
zoom_xs_1 = np.arange(-70, zoom_w_half, 0.1)
zoom_xs_2 = np.arange(zoom_w_half, 20, 0.1)

# Plotting

fig, ax = plt.subplots(figsize = (15,10))

ax.errorbar(bin_centres, det_eff, xerr = pEQW_err, yerr = det_err, fmt = 'k.', ls = 'none', capsize = 5)
ax.plot(xs_1, func(xs_1, w_half, par_1[1]), 'k')
ax.plot(xs_2, func(xs_2, w_half, par_2[1]), 'k')
ax.set_xlabel(r'Mean Equivalent Width, $\AA$')
ax.set_ylabel('Detection Efficiency')
ax.set_ylim(0,1.1)
ax.set_xlim(mods_functions.round_down(bin_centres[0], 50), 50)
ax.tick_params(top = True, right = True, direction = 'in')


left, bottom, width, height = [0.2, 0.2, 0.5, 0.4]
axins = fig.add_axes([left, bottom, width, height])
axins.errorbar(zoom_bin_centres, zoom_det_eff, xerr = zoom_pEQW_err, yerr = zoom_det_err, fmt = 'k.', ls = 'none', capsize = 5)
axins.plot(zoom_xs_1, func(zoom_xs_1, zoom_w_half, zoom_par_1[1]), 'k')
axins.plot(zoom_xs_2, func(zoom_xs_2, zoom_w_half, zoom_par_2[1]), 'k')
axins.tick_params(top = True, right = True, direction = 'in')
axins.set_xlim(mods_functions.round_down(zoom_bin_centres[0], 10), 10)
axins.set_ylim(0,1.1)

#plt.savefig("FeVII_detection_efficiency_plot.png")

plt.show()