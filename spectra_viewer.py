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

import mpl_style
plt.style.use(mpl_style.scientific)

User = 'Joe'
User_Config = Hirogen_Functions.user_config(user=User)

TableID = "SDSS_ECLE_Spectra"

# Search for 'QUERY' to find the actual database access location where parameters can be adjusted

Database = User_Config[0]
Database_User = User_Config[1]
Database_Password = User_Config[2]
Main_Spectra_Path = User_Config[3]

Data = Hirogen_Functions.database_connection(user=Database_User, password=Database_Password, database=Database)

all_lines = Hirogen_Functions.lines_for_analysis()
fe_lines = Hirogen_Functions.coronal_lines()
line_labels = Hirogen_Functions.presentation_zoom_line_labels()

short_name = 'SDSS J0952'
IDs = (0,2)

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
    f"WHERE SDSS_Very_ShortName = 'SDSS J0952' and Follow_Up_ID in {IDs} "
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

filepaths = Hirogen_Functions.sdss_spectra_file_path_generator(Main_Spectra_Path, Plate_List, MJD_List, FiberID_List, Survey_List, run2d_List, Path_Override_Flag_List, Path_Override_List)

new_spectra = Hirogen_Functions.ascii_spectrum_reader(filepaths[0], Redshift_List[0], Extinction_List[0], smoothing=False, sav_gol=True)

new_wave, new_flux = new_spectra[0],new_spectra[1]

og_spectra = Hirogen_Functions.sdss_spectrum_reader(filepaths[1], Redshift_List[1], Extinction_List[1], smoothing=False, sav_gol=True)

og_wave, og_flux = og_spectra[0],og_spectra[1]


try:
    og_Lower_Spec_Wave_Limit = Hirogen_Functions.round_down_to_nearest(np.min(np.array(og_wave)), 250)
except ValueError:
    og_Lower_Spec_Wave_Limit = 0

try:
    og_Upper_Spec_Wave_Limit = Hirogen_Functions.round_up_to_nearest(np.max(np.array(og_wave)), 250)
except ValueError:
    og_Upper_Spec_Wave_Limit = 10000

try:
    new_Lower_Spec_Wave_Limit = Hirogen_Functions.round_down_to_nearest(np.min(np.array(new_wave)), 250)
except ValueError:
    new_Lower_Spec_Wave_Limit = 0

try:
    new_Upper_Spec_Wave_Limit = Hirogen_Functions.round_up_to_nearest(np.max(np.array(new_wave)), 250)
except ValueError:
    new_Upper_Spec_Wave_Limit = 10000

Lower_Spec_Wave_Limit = np.min((og_Lower_Spec_Wave_Limit, new_Lower_Spec_Wave_Limit))
Upper_Spec_Wave_Limit = np.max((og_Upper_Spec_Wave_Limit, new_Upper_Spec_Wave_Limit))


fig, (ax1, ax2) = plt.subplots(2,1, sharex = True, sharey = True, figsize = (20,15))

line_locs = []

for i, line in enumerate(fe_lines):
    line_loc = all_lines[line][0]
    val, idx = Hirogen_Functions.find_nearest(og_wave, line_loc)
    ax1.axvline(val, linestyle = '--', color = 'r')
    ax2.axvline(val, linestyle = '--', color = 'r')
    line_locs.append(line_loc)
    

labels_ax = ax1.twiny()
labels_ax.set_xticks(line_locs)
labels_ax.set_xticklabels(line_labels, rotation = 15, fontsize = 18)
labels_ax.tick_params(top=False, bottom=False, left=False, right=False)
labels_ax.set_xlim(xmin=Lower_Spec_Wave_Limit, xmax=Upper_Spec_Wave_Limit)

ax1.plot(og_wave, og_flux, 'k')
ax1.set_xlim(xmin=Lower_Spec_Wave_Limit, xmax=Upper_Spec_Wave_Limit)

ax2.plot(new_wave, new_flux, 'k')
ax2.set_xlim(xmin=Lower_Spec_Wave_Limit, xmax=Upper_Spec_Wave_Limit)
ax2.set_xlabel(r'Rest Wavelength [$\AA$]')

fig.text(0.09, 0.5, r'Flux [10$^{-17}$ erg/cm$^{2}$/s/$\AA$]', ha='center', va='center', rotation='vertical', fontsize = 22)

plt.savefig("example_spectrum.png", bbox_inches='tight')

plt.show()