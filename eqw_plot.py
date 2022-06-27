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
from scipy.optimize import curve_fit
from pynverse import inversefunc

import mpl_style
plt.style.use(mpl_style.scientific)

c = constants.c
H0 = Planck18.H(0)
H0_diff = H0.to(u.m/u.m/u.s)

User = 'Joe'
User_Config = Hirogen_Functions.user_config(user=User)

Database = User_Config[0]
Database_User = User_Config[1]
Database_Password = User_Config[2]
Main_Spectra_Path = User_Config[3]


# Select which fakes to use

Fakes_TableID = "SDSS_Sample_Spectra"

Data = Hirogen_Functions.database_connection(user=Database_User, password=Database_Password, database=Database)

print("\n")
cursor = Data.cursor()

cursor.execute(
    f"SELECT DR16_Spectroscopic_ID, lin_con_pEQW_FeVII6008, lin_con_pEQW_FeX6376, lin_con_pEQW_FeXI7894, lin_con_pEQW_FeXIV5304, "
    f"Strong_FeVII_Flag, Strong_FeX_Flag, Strong_FeXI_Flag, Strong_FeXIV_Flag, ECLE_Candidate_Score, Candidate_Flag "
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
    candidate_flag = [item[10] for item in Candidate_Data]
        

else:
    print("Candidate_Data Length error: Check and try again.")

    sys.exit()

# Determine if object is candidate
ecle_candidate = np.zeros(len(ecle_candidate_score))
for i, score in enumerate(ecle_candidate_score):
    if score >= 7 or (fevii_flag[i]+fex_flag[i]+fexi_flag[i]+fexiv_flag[i])> 0:
        ecle_candidate[i] = 1

idx = np.where(ecle_candidate == 1)

print(idx)
# Remove nans so won't affect plots
fevii_pEQW = np.array(fevii_pEQW)
nans = np.argwhere(fevii_pEQW <= -999)
fevii_pEQW[nans] = 0

fex_pEQW = np.array(fex_pEQW)
nans = np.argwhere(fex_pEQW <= -999)
fex_pEQW[nans] = 0

fexi_pEQW = np.array(fexi_pEQW)
nans = np.argwhere(fexi_pEQW <= -999)
fexi_pEQW[nans] = 0

fexiv_pEQW = np.array(fexiv_pEQW)
nans = np.argwhere(fexiv_pEQW <= -999)
fexiv_pEQW[nans] = 0

# Calculate average pEQW

total_pEQW = (fevii_pEQW + fex_pEQW + fexi_pEQW + fexiv_pEQW)/4


total_sort_ind = np.argsort(total_pEQW)
total_pEQW_sort = total_pEQW[total_sort_ind]
scores_sort = ecle_candidate[total_sort_ind]

# Main Plot
# Bin pEQWs and calculate detection efficiencies

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

fevii_det_eff = (total_scores)/(bin_sizes)

fevii_det_error = np.sqrt(((total_scores+1)*(total_scores+2))/((bin_sizes+2)*(bin_sizes+3)) - ((total_scores+1)**2)/((bin_sizes+2)**2))

fevii_det_err = np.zeros((2, len(fevii_det_error)))

for i, fevii_error in enumerate(fevii_det_error):
    fevii_det_err[0][i] = fevii_error
    fevii_det_err[1][i] = fevii_error

    if fevii_det_eff[i] - fevii_det_err[0][i] < 0:
        fevii_det_err[0][i] = fevii_det_eff[i]

    if fevii_det_eff[i] + fevii_det_err[1][i] > 1:
        fevii_det_err[1][i] = 1 - fevii_det_eff[i]

fevii_pEQW_err = np.full(len(fevii_det_error), np.sqrt(bin_width))

fevii_bin_centres = np.zeros(len(fevii_det_eff))
for i in range(len(fevii_bin_centres)):
    fevii_bin_centres[i] = (total_bins[i] + total_bins[i+1])/2

not_nans = ~np.isnan(fevii_det_eff)
fevii_fin_det_eff = fevii_det_eff[not_nans]
fevii_fin_det_err = fevii_det_error[not_nans]
fevii_fin_bin_centres = fevii_bin_centres[not_nans]


fevii_fin_det_err[fevii_fin_bin_centres < -200] = 10**-20

# Inset

zoom_bin_width = 1

zoom_bins = mods_functions.create_bins(total_pEQW_sort[0], total_pEQW_sort[-1], zoom_bin_width)

zoom_info = pd.cut(total_pEQW_sort, zoom_bins)

zoom_cuts = np.cumsum(zoom_info.value_counts())[:-1]

zoom_scores_binned = np.split(scores_sort, zoom_cuts)

zoom_scores = np.zeros(len(zoom_scores_binned))

for i, array in enumerate(zoom_scores_binned):
    zoom_scores[i] = sum(array)

zoom_bin_sizes = np.zeros(len(zoom_scores_binned))
for i, array in enumerate(zoom_scores_binned):
    zoom_bin_sizes[i] = len(array)

fevii_zoom_det_eff = zoom_scores/zoom_bin_sizes

fevii_zoom_det_error = np.sqrt(((zoom_scores+1)*(zoom_scores+2))/((zoom_bin_sizes+2)*(zoom_bin_sizes+3)) - ((zoom_scores+1)**2)/((zoom_bin_sizes+2)**2))

fevii_zoom_det_err = np.zeros((2, len(fevii_zoom_det_error)))

for i, error in enumerate(fevii_zoom_det_error):
    fevii_zoom_det_err[0][i] = error
    fevii_zoom_det_err[1][i] = error

    if fevii_zoom_det_eff[i] - fevii_zoom_det_err[0][i] < 0:
        fevii_zoom_det_err[0][i] = fevii_zoom_det_eff[i]

    if fevii_zoom_det_eff[i] + fevii_zoom_det_err[1][i] > 1:
        fevii_zoom_det_err[1][i] = 1 - fevii_zoom_det_eff[i]

fevii_zoom_pEQW_err = np.full(len(fevii_zoom_det_error), np.sqrt(zoom_bin_width))

fevii_zoom_bin_centres = np.zeros(len(fevii_zoom_det_eff))
for i in range(len(fevii_zoom_bin_centres)):
    fevii_zoom_bin_centres[i] = (zoom_bins[i] + zoom_bins[i+1])/2

# Curves

def sigmoid(x, A, K, B, v, Q):
    return A+((K-A)/((1+Q*np.exp(-B*x))**(1/v)))

# Fit curve parameters to data
par, cov = curve_fit(sigmoid, fevii_fin_bin_centres, fevii_fin_det_eff, sigma = fevii_fin_det_err, bounds = ((0, 0, -np.inf, -np.inf, -np.inf,), (1, 1, np.inf, np.inf, np.inf)), maxfev = 5000)

xs = np.arange(-500, 50, 0.1)
zoom_xs =np.arange(-50, 10, 0.1)

#np.savetxt("fevii_parameters.csv", fevii_par)

"""fitted_sig = lambda x: par[0]+((par[1]-par[0])/((1+par[4]*np.exp(-par[2]*x))**(1/par[3])))

invsig = inversefunc(fitted_sig)

half_det_eff = invsig(0.5)

print(half_det_eff)"""

# Plotting

fig, ax = plt.subplots(figsize = (20,12))

ax.errorbar(fevii_bin_centres, fevii_det_eff, yerr = fevii_det_err, xerr = fevii_pEQW_err, fmt = '.k', ls = 'none', capsize = 3)
ax.plot(xs, sigmoid(xs, par[0], par[1], par[2], par[3], par[4]), 'k')

ax.set_xlabel(r'Average Equivalent Width of Coronal Lines [$\AA$]')
ax.set_ylim(0,1.1)
ax.set_xlim(mods_functions.round_down(fevii_bin_centres[0], 50), 50)
ax.set_ylabel('Detection Efficiency')


left, bottom, width, height = [0.17, 0.2, 0.5, 0.4]
axins = fig.add_axes([left, bottom, width, height])

axins.errorbar(fevii_zoom_bin_centres, fevii_zoom_det_eff, xerr = fevii_zoom_pEQW_err, yerr = fevii_zoom_det_err, fmt = '.k', ls = 'none', capsize = 3)
axins.plot(zoom_xs, sigmoid(zoom_xs, par[0], par[1], par[2], par[3], par[4]), 'k')

axins.set_xlim(-50, 10)
axins.set_ylim(0,1.1)

#axins.vlines(half_det_eff, 0, 0.5, ls = '--', color = 'k')
#axins.hlines(0.5, -50, half_det_eff, ls = '--', color = 'k')


#plt.savefig("detection_efficiency.png", bbox_inches='tight')


plt.show()