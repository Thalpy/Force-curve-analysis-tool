# -*- coding: utf-8 -*-
"""
Created on Feb 2019

@author: Thalpy + Thalpy
V.7.0 Edited 10.8.20
CURRENT MASTER - THIS IS THE FILE TO EDIT
"""

import pandas as pd
import numpy as np
import os
import shutil
import warnings
import logging
import datetime
import matplotlib
import matplotlib.pyplot as plt
import sys
import gc

matplotlib.use('Agg')

import pyAFM_FC as afm

# Assuming the parameters are passed in the following order:
# name, k_c, fitbin, cfit_min, cfit_max, cthresh, dfit_win, dfit_off, ext, approach
# Example command line call:
# python run_mfp_analysis.py S3/550mMS1.0.5Hz.nD.12nN 0.156 10 15 40 50 40 50 .jpg True

# Command line arguments
if len(sys.argv) < 13:
    print("Not enough arguments provided.")
    sys.exit(1)

# The first argument (sys.argv[0]) is the script name, so it is skipped.
name = sys.argv[1]
k_c = float(sys.argv[2])
fitbin = int(sys.argv[3])
cfit_min = int(sys.argv[4])
cfit_max = int(sys.argv[5])
cthresh = int(sys.argv[6])
dfit_win = int(sys.argv[7])
dfit_off = int(sys.argv[8])
ext = sys.argv[9]
out = sys.argv[10].lower() == 'true'
clear = sys.argv[11].lower() == 'true'
approach = sys.argv[12].lower() == 'true'


# Make sure we don't have .txt at the end of the name
if name[-4:] == '.txt':
    name = name[:-4]

# Autoprocessed stuff
outdir = name + '_force_curves/'

#Autoprocessed stuff
#Site = "S1" #unused
prefix = 'Run' #unused
outdir = name+'_force_curves/'
#txtdir = name+"_force_curves"+os.sep+prefix

#Check if Run folder exists

for i in range (1, 101):
    txtdir = name+"_force_curves"+os.sep+(prefix+str(i)+"/")
    if not os.path.exists(txtdir):
        prefix += str(i)+"/"
        break
    
    elif i == 100:
        print("Over 100 Runfolders found, overwriting run 100.")
    


#Extract processing
if not os.path.exists(outdir):
    res = afm.split_curves(name, k_c = k_c, ext = ext) # Works
    
#Initial run:
if os.path.isfile(txtdir+'approach_force_curves.csv'):
    res_df = pd.read_csv(txtdir+'approach_force_curves.csv')
else:
    slope = deriv_curves = afm.comp_def_deriv(name, cfit_min = cfit_min, cfit_max = cfit_max,
        dfit_win = dfit_win, dfit_off = dfit_off, fitbin =fitbin, k_c = k_c, prefix = prefix, out = True, ext = ext, clear = True)

    plt.close('all')
    gc.collect()
    matplotlib.use('Agg')

    res_df, param_bin = afm.proc_force_sep(name, cfit_min = cfit_min, cfit_max = cfit_max,
        dfit_win = dfit_win, dfit_off = dfit_off, fitbin =fitbin, k_c = k_c, prefix = prefix,
        out = True, ext = ext, clear = False, min_cont_pts=3, slope = slope, approach = approach)

plt.close('all')
gc.collect()

afm.plot_force_sep_res(name, res_df, prefix = prefix, binsize = 1.0, ext = ext)
afm.plot_force_h_scaling(name, res_df, prefix = prefix, binsize = 0.5, ext = ext)
Fc = afm.get_contact_forces(name, res_df, prefix = prefix, binsize = 1.0, ext = ext)

prefix = "cache/"

#Scroll along binfit
binfit_list = [0,0,Fc,0,0]
xbinfit_list = range(-2,3)
print(Fc)
_fitbin = 0
for i in range(-2,3):
    if i == 0:
        continue
    
    #Replace binfit with adjusted value
    _fitbin = fitbin + i
    
    res_df, param_bin = afm.proc_force_sep(name, cfit_min = cfit_min, cfit_max = cfit_max,
        dfit_win = dfit_win, dfit_off = dfit_off, fitbin =_fitbin, k_c = k_c, prefix = prefix,
        out = True, ext = ext, clear = False, min_cont_pts=3, slope = slope, approach = approach)
    
    
    binfit_list[i+2] = afm.get_contact_forces(name, res_df, prefix = prefix, binsize = 1.0, ext = ext)
    
    
#Plot graph
fig1, ax1 = plt.subplots(figsize=(10,4))
ax1.plot(xbinfit_list, binfit_list, '+-')
ax1.set_ylabel("Force at contact (N)")
ax1.set_xlabel("Binsize range (+/- 2)")
#ax1.set_xlim([-200,np.ceil(cfit_max+30)])
plt.tight_layout()
fig1.savefig(outdir+"_binsweep"+ext)
plt.close(fig1)

#Scroll along cfit_min
cmin_list = [0,0,0,0,0,Fc,0,0,0,0,0]
xc_list = range(-5,6)
print(Fc)
_cmin = 0
for i in range(-5,6):
    if i == 0:
        continue
    
    #Replace cfit with adjusted value
    _cmin = cfit_min + i
    
    res_df, param_bin = afm.proc_force_sep(name, cfit_min = _cmin, cfit_max = cfit_max,
        dfit_win = dfit_win, dfit_off = dfit_off, fitbin = fitbin, k_c = k_c, prefix = prefix,
        out = True, ext = ext, clear = False, min_cont_pts=3, slope = slope, approach = approach)
    
    
    cmin_list[i+5] = afm.get_contact_forces(name, res_df, prefix = prefix, binsize = 1.0, ext = ext)

#Plot graph
fig1, ax1 = plt.subplots(figsize=(10,4))
ax1.plot(xc_list, cmin_list, '+-')
ax1.set_ylabel("Force at contact (N)")
ax1.set_xlabel("Minimum contact parameter range (+/- 5)")
#ax1.set_xlim([-200,np.ceil(cfit_max+30)])
plt.tight_layout()
fig1.savefig(outdir+"_cminsweep"+ext)
plt.close(fig1)

#Scroll along cfit_max
cmax_list = [0,0,0,0,0,Fc,0,0,0,0,0]
print(Fc)
_cmax = 0
for i in range(-5,6):
    if i == 0:
        continue
    
    #Replace cfit with adjusted value
    _cmax = cfit_max + i
    
    res_df, param_bin = afm.proc_force_sep(name, cfit_min = cfit_min, cfit_max = _cmax,
        dfit_win = dfit_win, dfit_off = dfit_off, fitbin = fitbin, k_c = k_c, prefix = prefix,
        out = True, ext = ext, clear = False, min_cont_pts=3, slope = slope, approach = approach)
    
    
    cmax_list[i+5] = afm.get_contact_forces(name, res_df, prefix = prefix, binsize = 1.0, ext = ext)


#Plot graph
fig1, ax1 = plt.subplots(figsize=(10,4))
ax1.plot(xc_list, cmax_list, '+-')
ax1.set_ylabel("Force at contact (N)")
ax1.set_xlabel("Maximum contact parameter range (+/- 5)")
#ax1.set_xlim([-200,np.ceil(cfit_max+30)])
plt.tight_layout()
fig1.savefig(outdir+"_cmaxsweep"+ext)
plt.close(fig1)