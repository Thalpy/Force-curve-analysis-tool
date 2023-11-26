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

name = 'S1/T6.2.550mM.Si.S1deflection'
prefix = 'Run3_Ret/'
outdir = name+'_force_curves/'
txtdir = name+"_force_curves"+os.sep+'ret_'+prefix
k_c = 0.156
fitbin = 20
cfit_min = 30
cfit_max = 50
cthresh = 50
dfit_win = 80
dfit_off = 85
ext = '.jpg'
out = True
approach = False

if not os.path.exists(outdir):
    res = afm.split_curves(name, k_c = k_c, ext = ext)


if os.path.isfile(txtdir+'retract_force_curves.csv'):
    res_df = pd.read_csv(txtdir+'retract_force_curves.csv')
else:
    slope = afm.comp_def_deriv(name, cfit_min = cfit_min, cfit_max = cfit_max,
        cthresh = cthresh, dfit_win = dfit_win, dfit_off = dfit_off, fitbin =fitbin, k_c = k_c,
        prefix = prefix, out = True, ext = ext, clear = True, approach = approach)
    res_df, param_bin = afm.proc_force_sep(name, cfit_min = cfit_min, cfit_max = cfit_max,
        cthresh = cthresh,dfit_win = dfit_win, dfit_off = dfit_off, fitbin =fitbin, k_c = k_c,
        prefix = prefix, out = True, ext = ext, clear = False, min_cont_pts=3, slope = slope, approach = approach)


#afm.plot_force_sep_res(name, res_df, prefix = prefix, binsize = 1.0, ext = '.pdf', approach = False)
#afm.fit_force_sep(name, res_df, prefix = prefix, binsize = 1.0, ext = '.pdf')
force_curves_c, force_curves_nc = afm.profile_attractive_forces(name, res_df, prefix = prefix, binsize = 1.0, ext = ext, approach = approach)
afm.plot_force_sep_c_nc(name, force_curves_c, force_curves_nc, prefix = prefix, binsize = 2.0, ext = ext, approach = approach)
