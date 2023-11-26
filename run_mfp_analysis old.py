# -*- coding: utf-8 -*-
"""
Created on Feb 2019

@author: Thalpy
"""


import pandas as pd
import numpy as np
import os
import shutil
import warnings
import logging

import pyAFM_FC as afm

name = 'S3/550mMS1.0.5Hz.nD.12nN'
prefix = 'Run15/'
outdir = name+'_force_curves/'
txtdir = name+"_force_curves"+os.sep+prefix
k_c = 0.156
fitbin = 5
cfit_min = 1
cfit_max = 25
cthresh = 50
dfit_win = 40
dfit_off = 65
ext = '.jpg'
out = True

#ThalpyParams
"""
0.6mM = 9.942266946
1.6mM = 6.088369852
100mM = 0.770124686
25mM = 1.540249372
5mM = 3.444102299
10mM = 2.435347941
230mM = 0.507805149
550mM = 0.328382249
50mM = 1.089120709
"""
h = 9.942266946
hGraphs = False
squareRange = 20
hLow = 85
hHigh = 135
linHigh = 155

if not os.path.exists(outdir):
    res = afm.split_curves(name, k_c = k_c, ext = ext)


if os.path.isfile(txtdir+'approach_force_curves.csv'):
    temp = input("Path exists, proceed?")
    res_df = pd.read_csv(txtdir+'approach_force_curves.csv')
else:
    slope = deriv_curves = afm.comp_def_deriv(name, cfit_min = cfit_min, cfit_max = cfit_max,
        cthresh = cthresh, dfit_win = dfit_win, dfit_off = dfit_off, fitbin =fitbin, k_c = k_c,
        prefix = prefix, out = True, ext = ext, clear = True)

    res_df, param_bin = afm.proc_force_sep(name, cfit_min = cfit_min, cfit_max = cfit_max,
        cthresh = cthresh,dfit_win = dfit_win, dfit_off = dfit_off, fitbin =fitbin, k_c = k_c,
        prefix = prefix, out = True, ext = ext, clear = False, min_cont_pts=3, slope = slope, 
        h = h, hGraphs = hGraphs, squareRange = squareRange, hLow = hLow, hHigh = hHigh)


afm.plot_force_sep_res(name, res_df, prefix = prefix, binsize = 1.0, ext = '.jpg')
afm.get_contact_forces(name, res_df, prefix = prefix, binsize = 1.0, ext = '.jpg')

#Check for pygame and disable sound if not present. 
try:
    import pygame
    #for sfx:
    pygame.init()
    pygame.mixer.init()
    mysound = pygame.mixer.Sound("pyAFM_FC/bleep_high.ogg")
    mysound.play()
except:
    print("Completed")
