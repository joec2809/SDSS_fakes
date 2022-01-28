import matplotlib.pyplot as plt

import Hirogen_Functions
import mods_functions

import sys
import os
import math

import numpy as np
import scipy.constants as constants
import mysql.connector
import decimal

from astropy.io import fits
from astropy import units as u
from mysql.connector import errorcode
from astropy.convolution import convolve, Box1DKernel
from random import uniform

number = -175
factor = 20

print(mods_functions.next_lowest_multiple(number, factor))