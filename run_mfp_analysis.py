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

def remove_skipped_files(folder_path, skip_list):
    removed_dir = os.path.join(folder_path, 'Removed')
    if not os.path.exists(removed_dir):
        os.makedirs(removed_dir)
    
    files = os.listdir(folder_path)
    for file in files:
        if file.endswith('.txt'):
            # Extract the file number
            file_number = int(''.join(filter(str.isdigit, file)))
            # Check if the file number is in the skip list
            if file_number in skip_list:
                # Move the file to the 'Removed' directory
                shutil.move(os.path.join(folder_path, file), os.path.join(removed_dir, file))
                print(f"File {file} moved to 'Removed' directory.")


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
binsize = float(sys.argv[9])  # This was the 8th parameter in the GUI args list
ext = sys.argv[10]
out = sys.argv[11].lower() == 'true'
clear = sys.argv[12].lower() == 'true'
approach = sys.argv[13].lower() == 'true'
extra = sys.argv[14].lower() == 'true'
dwell = sys.argv[15].lower() == 'true'

# Make sure we don't have .txt at the end of the name
if name[-4:] == '.txt':
    name = name[:-4]

# Autoprocessed stuff
outdir = name + '_force_curves/'

#Autoprocessed stuff
#Site = "S1" #unused
if approach:
    prefix = 'Run'
else:
    prefix = 'ret_Run_Ret'
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
    
#Check if approach or retract folder is there
if not os.path.exists(outdir + "approach/") or not os.path.exists(outdir + "retract/"):
    # Use split_dwell_curves if dwell is True, else use split_curves
    if dwell:
        print("Dwell is enabled. Running split_dwell_curves.")
        res = afm.split_dwell_curves2(name, k_c=k_c, ext=ext)
    else:
        print("Dwell is disabled. Running split_curves.")
        res = afm.split_curves(name, k_c=k_c, ext=ext)

# Skip designated files
# List of file numbers to skip
skip_list = []  # Update this list as needed

# Assuming 'approach' variable is True if processing approach files, False for retract files
approach_or_retract_folder = "approach" if approach else "retract"
approach_or_retract_path = os.path.join(outdir, approach_or_retract_folder)

# Call the function to remove skipped files
remove_skipped_files(approach_or_retract_path, skip_list)

#Initial run:
if os.path.isfile(txtdir+'approach_force_curves.csv'):
    res_df = pd.read_csv(txtdir+'approach_force_curves.csv')
else:
    slope = deriv_curves = afm.comp_def_deriv(name, cfit_min = cfit_min, cfit_max = cfit_max, binsize = binsize,
        dfit_win = dfit_win, dfit_off = dfit_off, fitbin =fitbin, k_c = k_c, prefix = prefix, out = True, ext = ext, clear = True, approach = approach)

    plt.close('all')
    gc.collect()
    matplotlib.use('Agg')

    res_df, param_bin = afm.proc_force_sep(name, cfit_min = cfit_min, cfit_max = cfit_max, binsize = binsize,
        dfit_win = dfit_win, dfit_off = dfit_off, fitbin =fitbin, k_c = k_c, prefix = prefix,
        out = True, ext = ext, clear = False, min_cont_pts=3, slope = slope, approach = approach)

plt.close('all')
gc.collect()

afm.plot_force_sep_res(name, res_df, prefix = prefix, binsize = binsize, ext = ext, approach = approach)
afm.plot_force_h_scaling(name, res_df, prefix = prefix, binsize = binsize, ext = ext, approach = approach)
Fc, fc_stdev = afm.get_contact_forces(name, res_df, prefix = prefix, ext = ext, approach = approach)
force_curves_c, force_curves_nc = afm.profile_attractive_forces(name, res_df, prefix = prefix, binsize = binsize, ext = ext, approach = approach)  
afm.plot_force_sep_c_nc(name, force_curves_c, force_curves_nc, prefix = prefix, binsize = 2.0, ext = ext, approach = approach)

prefix = "cache/"

if extra: # If you can't get a good fit and want to run this 25 times in a row out of desparation.
    #Scroll along binfit
    binfit_list = [0,0,Fc,0,0]
    stdev_binfit_list = [0,0,fc_stdev,0,0]
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
        
        
        binfit_list[i+2], stdev_binfit_list[i+2] = afm.get_contact_forces(name, res_df, prefix = prefix, binsize = 1.0, ext = ext)
                
        
    #Plot graph
    fig1, ax1 = plt.subplots(figsize=(10,4))
    ax1.errorbar(xbinfit_list, binfit_list, yerr=stdev_binfit_list, fmt='+-', ecolor='red', capsize=5)
    ax1.set_ylabel("Force at contact (N)")
    ax1.set_xlabel("Binsize range (+/- 2)")
    plt.tight_layout()
    fig1.savefig(outdir+"_binsweep"+ext)
    plt.close(fig1)


    #Scroll along cfit_min
    cmin_list = [0,0,0,0,0,Fc,0,0,0,0,0]
    cmin_stdev_list = [0,0,0,0,0,fc_stdev,0,0,0,0,0]
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
        
        
        cmin_list[i+5], cmin_stdev_list[i+5] = afm.get_contact_forces(name, res_df, prefix = prefix, binsize = 1.0, ext = ext)

    #Plot graph
    fig1, ax1 = plt.subplots(figsize=(10,4))
    ax1.errorbar(xc_list, cmin_list, yerr=cmin_stdev_list, fmt='+-', ecolor='red', capsize=5)
    ax1.set_ylabel("Force at contact (N)")
    ax1.set_xlabel("Minimum contact parameter range (+/- 5)")
    plt.tight_layout()
    fig1.savefig(outdir+"_cminsweep"+ext)
    plt.close(fig1)


    #Scroll along cfit_max
    cmax_list = [0,0,0,0,0,Fc,0,0,0,0,0]
    cmax_stdev_list = [0,0,0,0,0,fc_stdev,0,0,0,0,0]
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
        
        
        cmax_list[i+5], cmax_stdev_list[i+5] = afm.get_contact_forces(name, res_df, prefix = prefix, binsize = 1.0, ext = ext)


    #Plot graph
    fig1, ax1 = plt.subplots(figsize=(10,4))
    ax1.errorbar(xc_list, cmax_list, yerr=cmax_stdev_list, fmt='+-', ecolor='red', capsize=5)
    ax1.set_ylabel("Force at contact (N)")
    ax1.set_xlabel("Maximum contact parameter range (+/- 5)")
    plt.tight_layout()
    fig1.savefig(outdir+"_cmaxsweep"+ext)
    plt.close(fig1)
