import sys
import time
import warnings
import math
import os
from xml.dom.expatbuilder import CDATA_SECTION_NODE

from click import Path
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

User = 'Joe'
User_Config = Hirogen_Functions.user_config(user=User)

Database = User_Config[0]
Database_User = User_Config[1]
Database_Password = User_Config[2]
Main_Spectra_Path = User_Config[3]

plot = "FeVII"

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
    f"Strong_FeVII_Flag, Strong_FeX_Flag, Strong_FeXI_Flag, Strong_FeXIV_Flag, ECLE_Candidate_Score, lin_con_LineFlux_Halpha, z_SDSS_spec, "
    f"spec_Plate, spec_MJD, spec_FiberID, generalised_extinction, survey, run2d, Standard_Inclusion, Path_Override_Flag, Path_Override "
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
    halpha_flux = [item[10] for item in Candidate_Data]
    redshift = [item[11] for item in Candidate_Data]
    Plate_List = [item[12] for item in Candidate_Data]
    MJD_List = [item[13] for item in Candidate_Data]
    FiberID_List = [item[14] for item in Candidate_Data]
    Extinction_List = [item[15] for item in Candidate_Data]
    Survey_List = [item[16] for item in Candidate_Data]
    run2d_List = [item[17] for item in Candidate_Data]
    Standard_Inclusion_List = [item[18] for item in Candidate_Data]
    Path_Override_Flag_List = [item[19] for item in Candidate_Data]
    Path_Override_List = [item[20] for item in Candidate_Data]


FilePaths = Hirogen_Functions.sdss_peaks_spectra_file_path_generator(
    Main_Spectra_Path, Plate_List, MJD_List, FiberID_List, Survey_List, run2d_List, 
    Path_Override_Flag_List, Path_Override_List, 'fakes'
)

h_alpha_flux = np.zeros(len(Candidate_Data))

for i, object in enumerate(Candidate_Data):
    spectral_data = Hirogen_Functions.sdss_spectrum_reader(
        FilePaths[i],
        redshift[i],
        Extinction_List[i]
        )
    
    wave = spectral_data[0]
    flux = spectral_data[1]

    idx = mods_functions.find_nearest(wave.value, 6563)
    h_alpha_flux[i] = flux.value[idx]

    print(i)


if plot == "FeVII":
    np.savetxt("FeVII_H_alpha_flux.dat", h_alpha_flux)
elif plot == "Non FeVII":
    np.savetxt("Non_FeVII_H_alpha_flux.dat", h_alpha_flux)