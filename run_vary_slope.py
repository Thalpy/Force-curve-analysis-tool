# -*- coding: utf-8 -*-
"""
Created on Feb 2019

@author: jrr1
"""


import pandas as pd
import numpy as np
import os
import shutil
import warnings
import logging

import pyAFM_FC as afm

import matplotlib  as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import quiver, cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import csv
mpl.rcParams['lines.linewidth'] = 2

name = '2018-05-23'
prefix = 'test/'
outdir = name+'_force_curves/'
txtdir = name+"_force_curves"+os.sep+prefix
k_c = 0.188
fitbin = 7
cfit_min = 60
cfit_max = 120
cthresh = 50
dfit_win = 60
dfit_off = 65
ext = '.pdf'
out = True

cols = ['slope', 'f_c']
slope_f_c = pd.DataFrame(columns=cols)

resdir = name+"_force_curves"+os.sep
slopes = [1.375, 1.377, 1.379, 1.381,1.383, 1.385, 1.387, 1.389]

for s in slopes:

    prefix = 'slope_'+str(s)+os.sep
    res_df, param_bin = afm.proc_force_sep(name, cfit_min = cfit_min, cfit_max = cfit_max,
        cthresh = cthresh,dfit_win = dfit_win, dfit_off = dfit_off, fitbin =fitbin, k_c = k_c,
        prefix = prefix, out = True, ext = ext, clear = False, min_cont_pts=3, slope = s)


    afm.plot_force_sep_res(name, res_df, prefix = prefix, binsize = 1.0, ext = '.pdf')
    #afm.fit_force_sep(name, res_df, prefix = prefix, binsize = 1.0, ext = '.pdf')
    fc = afm.get_contact_forces(name, res_df, prefix = prefix, binsize = 1.0, ext = '.pdf')

    dfs = pd.DataFrame({
        'slope': [s],
        'f_c': [fc],
    })
    slope_f_c=slope_f_c.append(dfs, sort=True)

slope_f_c.to_csv(resdir+'slope_f_c.csv')

fig, ax = plt.subplots(figsize=(10,4))
ax.plot(slope_f_c['slope'], slope_f_c['f_c'], '-ko' )
ax.set_ylabel(r'Contact Force (nN)')
ax.set_xlabel(r'set deflection slope')
plt.tight_layout()
fig.savefig(resdir+'slope_f_c'+ext)
plt.close(fig)
