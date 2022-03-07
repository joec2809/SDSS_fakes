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
from scipy.interpolate import CubicSpline, interp1d
from scipy.optimize import curve_fit

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

plot = "Non FeVII"
TDE_type_dist = True

if TDE_type_dist:
    if plot == "FeVII":
        Fakes_TableID = "SDSS_TDE_Dist_FeVII_Fake_Spectra"
    elif plot == "Non FeVII":
        Fakes_TableID = "SDSS_TDE_Dist_Non_FeVII_Fake_Spectra"
    else:
        print("Table not set correctly")
        sys.exit()
else:
    if plot == "FeVII":
        Fakes_TableID = "SDSS_FeVII_Fake_Spectra"
    elif plot == "Non FeVII":
        Fakes_TableID = "SDSS_Non_FeVII_Fake_Spectra"
    else:
        print("Table not set correctly")
        sys.exit()

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


if plot == "FeVII":
    total_pEQW = (fevii_pEQW + fex_pEQW + fexi_pEQW + fexiv_pEQW)/4
elif plot == "Non FeVII":
    total_pEQW = (fex_pEQW + fexi_pEQW + fexiv_pEQW)/3

total_sort_ind = np.argsort(total_pEQW)
total_pEQW_sort = total_pEQW[total_sort_ind]
scores_sort = ecle_candidate[total_sort_ind]

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
fin_det_err = det_error[not_nans]
fin_bin_centres = bin_centres[not_nans]

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
zoom_fin_det_err = zoom_det_error[zoom_not_nans]
zoom_fin_bin_centres = zoom_bin_centres[zoom_not_nans]

# Curves

def richards_curve(x, A, K, B, v, Q):
    return A+((K-A)/((1+Q*np.exp(B*x))**(1/v)))

par_rich, cov_rich = curve_fit(richards_curve, fin_bin_centres, fin_det_eff, sigma = fin_det_err)
zoom_par_rich, zoom_cov_rich = curve_fit(richards_curve, zoom_fin_bin_centres, zoom_fin_det_eff, sigma = zoom_fin_det_err)

xs = np.arange(-400, 50, 0.1)
zoom_xs =np.arange(-50, 10, 0.1)

# Plotting

fig, ax = plt.subplots(figsize = (20,10))

ax.errorbar(bin_centres, det_eff, xerr = pEQW_err, yerr = det_err, fmt = 'k.', ls = 'none', capsize = 5)

ax.plot(xs, richards_curve(xs, par_rich[0], par_rich[1], par_rich[2], par_rich[3], par_rich[4]), 'k')

if plot == "FeVII":
    ax.set_title('Detection efficiency for fakes generated using spectra with FeVII lines')
elif plot == "Non FeVII":
    ax.set_title('Detection efficiency for fakes generated using spectra without FeVII lines')
ax.set_xlabel(r'Average Equivalent Width of Coronal Lines, $\AA$')
ax.set_ylabel('Detection Efficiency')
ax.set_ylim(0,1.1)
ax.set_xlim(mods_functions.round_down(bin_centres[0], 50), 50)


left, bottom, width, height = [0.17, 0.2, 0.5, 0.4]
axins = fig.add_axes([left, bottom, width, height])

axins.errorbar(zoom_bin_centres, zoom_det_eff, xerr = zoom_pEQW_err, yerr = zoom_det_err, fmt = 'k.', ls = 'none', capsize = 5)

axins.plot(zoom_xs, richards_curve(zoom_xs, zoom_par_rich[0], zoom_par_rich[1], zoom_par_rich[2], zoom_par_rich[3], zoom_par_rich[4]), 'k')

axins.set_xlim(-50, 10)
axins.set_ylim(0,1.1)

if TDE_type_dist:
    if plot == "FeVII":
        plt.savefig("FeVII_detection_efficiency_TDE_dist.png")
    elif plot == "Non FeVII":
        plt.savefig("Non_FeVII_detection_efficiency_TDE_dist.png")
else:
    if plot == "FeVII":
        plt.savefig("FeVII_detection_efficiency.png")
    elif plot == "Non FeVII":
        plt.savefig("Non_FeVII_detection_efficiency.png")


plt.show()