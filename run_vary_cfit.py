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

cols = ['cmin', 'slope', 'f_c']
cm_f_c = pd.DataFrame(columns=cols)

resdir = name+"_force_curves"+os.sep
cmins = [50, 55, 60, 65]

for cm in cmins:

    prefix = 'cfit_min_'+str(cm)+os.sep

    slope = deriv_curves = afm.comp_def_deriv(name, cfit_min = cm, cfit_max = cm+60,
        cthresh = cthresh, dfit_win = dfit_win, dfit_off = dfit_off, fitbin =fitbin, k_c = k_c,
        prefix = prefix, out = True, ext = ext, clear = True)
    res_df, param_bin = afm.proc_force_sep(name, cfit_min = cm, cfit_max = cm+60,
        cthresh = cthresh,dfit_win = dfit_win, dfit_off = dfit_off, fitbin =fitbin, k_c = k_c,
        prefix = prefix, out = True, ext = ext, clear = False, min_cont_pts=3, slope = slope)


    afm.plot_force_sep_res(name, res_df, prefix = prefix, binsize = 1.0, ext = '.pdf')
    #afm.fit_force_sep(name, res_df, prefix = prefix, binsize = 1.0, ext = '.pdf')
    fc = afm.get_contact_forces(name, res_df, prefix = prefix, binsize = 1.0, ext = '.pdf')

    dfs = pd.DataFrame({
        'slope': [slope],
        'cmin': [cm],
        'f_c': [fc],
    })
    cm_f_c=cm_f_c.append(dfs, sort=True)

cm_f_c.to_csv(resdir+'cfit_min_f_c.csv')

fig, ax = plt.subplots(figsize=(10,4))
ax.plot(cm_f_c['cmin'], cm_f_c['f_c'], '-ko' )
ax.set_ylabel(r'Contact Force (nN)')
ax.set_xlabel(r'cfit_min')
plt.tight_layout()
fig.savefig(resdir+'cfit_min_f_c'+ext)
plt.close(fig)
