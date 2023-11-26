
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 22:53:17 2019
@author: Joe
Version PhD 1.0
"""

from __future__ import (absolute_import, division, print_function)

import pandas as pd
import numpy as np
import os
import shutil
import warnings
import logging
import datetime
import shlex, subprocess, os
from os.path import splitext
import glob
from scipy.optimize import curve_fit

import matplotlib  as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import quiver, cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import csv
from scipy.signal import savgol_filter

from .proc_force_sep import binscatter  # Importing the binscatter function from the proc_force_sep module
mpl.rcParams['lines.linewidth'] = 2  # Setting the default line width for plots

# Function to calculate DL (Double Layer) CR (Colloidal Rod) force
def F_DL_CR(y, lam_D=11.6, dp=1540.0, zeta=104, epsilon=50, p=0.5):
    f0 = np.pi * epsilon * (8.8542e-9) * (dp / lam_D) * np.power(zeta, 2)
    top = np.exp((-1) * y / lam_D) + (2 * p - 1) * np.exp((-2) * y / lam_D)
    bottom = 1 - np.power(2 * p - 1, 2) * np.exp((-2) * y / lam_D)
    return f0 * top / bottom

# Function to calculate van der Waals force
def F_vdw(y, AH=0.0025, dp=1540.0):
    return (-1) * AH * (dp / (12 * np.power(y, 2)))

# Function to calculate screening length
def screen_len(I, epsilon=50):  # I in mM
    return 0.304 * np.sqrt(epsilon / 68.0) * 1 / np.sqrt(I / 1000.0)

