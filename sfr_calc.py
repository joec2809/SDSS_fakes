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

c = constants.c * u.m/u.s
H0 = Planck18.H(0)
H0_diff = H0.to(u.m/u.m/u.s)
Om = Planck18.Om0
Ol = 1 - Om

q0 = 0.5*Om - Ol

User = 'Joe'
User_Config = Hirogen_Functions.user_config(user=User)

Database = User_Config[0]
Database_User = User_Config[1]
Database_Password = User_Config[2]
Main_Spectra_Path = User_Config[3]

table = "FeVII"

if table == "FeVII":
    Fakes_TableID = "SDSS_FeVII_Fake_Spectra"
elif table == "Non FeVII":
    Fakes_TableID = "SDSS_Non_FeVII_Fake_Spectra"
else:
    print("Table not set correctly")
    sys.exit()

Data = Hirogen_Functions.database_connection(user=Database_User, password=Database_Password, database=Database)

print("\n")
cursor = Data.cursor()

cursor.execute(
    f"SELECT DR16_Spectroscopic_ID, z_SDSS_spec " 
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
    redshift = [item[1] for item in Candidate_Data]

redshift = np.array(redshift)

lum_dist = (((redshift*c)/H0_diff)*(1 + (redshift/2)*(1-q0))*100).value

if table == "FeVII":
    h_alpha_flux = np.loadtxt("Non_FeVII_H_alpha_flux.dat")*6563*10**-17
elif table == "Non FeVII":
    h_alpha_flux = np.loadtxt("FeVII_H_alpha_flux.dat")*6563*10**-17


luminosity = Hirogen_Functions.flux_to_luminosity(h_alpha_flux, lum_dist)

sfr = 4.6*10**-42 * luminosity

print(min(sfr))