# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 2017

@author: Thalpy
Version PhD 1.0
"""

from __future__ import (absolute_import, division, print_function)

import pandas as pd
import numpy as np
import os
import warnings
import logging
from matplotlib.pyplot import cm

from matplotlib import pyplot as plt

import shlex, subprocess, os
from os.path import splitext
#from linkpost import *

import matplotlib  as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import quiver
from mpl_toolkits.axes_grid1 import make_axes_locatable
import csv
#from .proc_force_sep import *
#from .post_proc_mfp import *
mpl.rcParams['lines.linewidth'] = 2


def split_curves(name, k_c = 0.188, ext = '.pdf'):
    # Function to split up a bunch of force curves


    outdir = name+'_force_curves/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    forcedir = outdir+'Raw_force_curves/'
    if not os.path.exists(forcedir):
        os.makedirs(forcedir)

    atxtdir = outdir+'approach/'
    if not os.path.exists(atxtdir):
        os.makedirs(atxtdir)

    rtxtdir = outdir+'retract/'
    if not os.path.exists(rtxtdir):
        os.makedirs(rtxtdir)

    #dis_pos_full = np.loadtxt(name+'.txt', delimiter="\t", skiprows=1)
    dis_pos_full = np.genfromtxt(name+'.txt', delimiter="\t", skip_header=1)

    #removes pesky velocity curves in a hamfisted way:
    CheckType = np.genfromtxt(name+'.txt', delimiter="\t", max_rows=1, dtype=None)
    delcol=0
    while True:
        try:
            for i in range(len(dis_pos_full)-1):
                test = str(CheckType[i])
                if "Velo" in test:
                    dis_pos_full = np.delete(dis_pos_full, i-delcol, axis=1)
                    delcol = delcol+1
                    #print('Velo'+str(i))
                if "Position" in test:
                    #print(i)
                    dis_pos_full = np.delete(dis_pos_full, i-delcol, axis=1)
                    delcol = delcol+1
            break
        except:
            #print(len(dis_pos_full))
            break

    nrows, ncols = np.shape(dis_pos_full)

    #print(ncols)
    ncurves = int(ncols/2)

    for i in range(ncurves):
        try:
            #print(i)
            disp = dis_pos_full[:,2*i]
            disp = disp[~np.isnan(disp)]
            Zpiezo = dis_pos_full[:,2*i+1]
            Zpiezo = Zpiezo[~np.isnan(Zpiezo)]

            #print(len(disp), len(Zpiezo))
            #disp_offset = np.mean(disp[0:navg])
            disp_offset = 0
            force = k_c*(disp-disp_offset)*1.0e9

            Zpiezo = Zpiezo*1.0e6

            split = np.argmax(Zpiezo)
            # print(np.max(Zpiezo))
            #print(split)

            #plot stuff and save files

            fig, ax = plt.subplots(figsize=(10,4), dpi=100)
            ax.plot(Zpiezo[0:split-1],force[0:split-1] , 'r-')
            ax.plot(Zpiezo[split:-1],force[split:-1] , 'b-')
            ax.set_ylabel('Force (nN)')
            ax.set_xlabel('Zpiezo (um)')
            plt.tight_layout()
            fig.savefig(forcedir+str(i).zfill(3)+'_force_Zpiezo'+ext)
            plt.close(fig)

            np.savetxt(atxtdir+str(i).zfill(3)+'_approach.txt', np.vstack((Zpiezo[0:split-1], disp[0:split-1], force[0:split-1]) ).T, delimiter='\t', newline='\n', header='Zpiezo \t Defl \t Force')
            np.savetxt(rtxtdir+str(i).zfill(3)+'_retract.txt', np.vstack((Zpiezo[split:-1], disp[split:-1], force[split:-1]) ).T, delimiter='\t', newline='\n', header='Zpiezo \t Defl \t Force')

            # Zoom in on the transition region for the first curve as an example
            disp = dis_pos_full[:, 0]
            disp = disp[~np.isnan(disp)]
            Zpiezo = dis_pos_full[:, 1]
            Zpiezo = Zpiezo[~np.isnan(Zpiezo)]
            disp_offset = 0
            force = k_c * (disp - disp_offset) * 1e9
            Zpiezo = Zpiezo * 1e6  # Convert from meters to micrometers

            # Finding the transition point for zooming
            # Assuming the transition occurs at the point where the force starts increasing from the baseline
            baseline_noise = np.std(force[0:split-1])  # Calculate the standard deviation of the baseline noise
            baseline = np.mean(force[0:split-1])  # Calculate baseline force before the transition
            threshold = baseline + 5 * baseline_noise  # Set threshold as a multiple of noise level above the baseline

            transition_indices = np.where(force > threshold)[0]
            transition_index = transition_indices[0] if transition_indices.size > 0 else split
            zoom_window=100e-3  # Zoom window size in micrometers

            """ broken atm
            # Define the zoom region around the transition
            half_window = int(zoom_window / (Zpiezo[1] - Zpiezo[0]) / 2)
            start_index = max(transition_index - half_window, 0)
            end_index = min(transition_index + half_window, len(Zpiezo))

            # Create a plot for the zoomed region
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(Zpiezo[start_index:end_index], force[start_index:end_index], 'b-')
            ax.set_ylabel('Force (nN)')
            ax.set_xlabel('Zpiezo (um)')
            plt.tight_layout()
            #fig.savefig(forcedir+str(i).zfill(3)+'_force_Zpiezo_zoom'+ext) #Broken atm
            """
        except:
            print('!!! Error in curve '+str(i))

def validate_curve_shape(curve_data, force_min=-10, force_max=20, zpiezo_min_range=5e-7):
    """
    Validates the shape of the force curve.
    
    :param curve_data: DataFrame containing the curve data with 'Force' and 'Zpiezo' columns.
    :param force_min: The minimum acceptable force value in newtons.
    :param force_max: The maximum acceptable force value in newtons.
    :param zpiezo_min_range: The minimum acceptable range of Zpiezo values in meters.
    :return: True if the curve is valid, False otherwise.
    """
    #print("Force min:", curve_data['Force'].min(), "Force max:", curve_data['Force'].max())
    #print("Zpiezo min:", curve_data['Zpiezo'].min(), "Zpiezo max:", curve_data['Zpiezo'].max())


    # Check if the force range is within the specified window
    force_range_condition = curve_data['Force'].between(force_min, force_max).all()

    # Check if the Zpiezo range is at least 1um
    zpiezo_range_condition = (curve_data['Zpiezo'].max() - curve_data['Zpiezo'].min()) >= zpiezo_min_range

    # If both conditions are met, the curve is valid
    return force_range_condition and zpiezo_range_condition

def split_dwell_curves(name, k_c=0.188, ext='.pdf'):
    # Read the data file
    data = pd.read_csv(name + '.txt', sep="\t", header=0)

    # Correcting column names if they contain '#'
    data.columns = data.columns.str.strip('# ').str.strip()

    # Create directories to save the split curves
    outdir = name + '_force_curves/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    atxtdir = outdir + 'approach/'
    if not os.path.exists(atxtdir):
        os.makedirs(atxtdir)

    rtxtdir = outdir + 'retract/'
    if not os.path.exists(rtxtdir):
        os.makedirs(rtxtdir)

    forcedir = outdir+'Raw_force_curves/'
    if not os.path.exists(forcedir):
        os.makedirs(forcedir)

    # Iterate over each curve in the dataset
    nrows, ncols = np.shape(dis_pos_full)

        #removes pesky velocity curves in a hamfisted way:
    CheckType = np.genfromtxt(name+'.txt', delimiter="\t", max_rows=1, dtype=None)
    delcol=0
    while True:
        try:
            for i in range(len(dis_pos_full)-1):
                test = str(CheckType[i])
                if "Velo" in test:
                    dis_pos_full = np.delete(dis_pos_full, i-delcol, axis=1)
                    delcol = delcol+1
                    #print('Velo'+str(i))
                if "Position" in test:
                    #print(i)
                    dis_pos_full = np.delete(dis_pos_full, i-delcol, axis=1)
                    delcol = delcol+1
            break
        except:
            #print(len(dis_pos_full))
            break

    #print(ncols)
    ncurves = int(ncols/2)
    for i in range(ncurves):
        try:
            Zpiezo = data.iloc[:, 2 * i] * 1e6  # Convert from meters to micrometers
            defl = data.iloc[:, 2 * i + 1]  # Deflection in meters
            force = k_c * defl * 1e9  # Convert deflection to force (nN)

            # Combine approach and retract data for validation
            curve_data = pd.DataFrame({'Zpiezo': Zpiezo, 'Defl': defl, 'Force': force})

            if not validate_curve_shape(curve_data):
                print(f"Curve {i} discarded due to not meeting validation criteria.")
                continue

            # Identify the turning point
            turning_point_index = Zpiezo.idxmax()       

            # Format the data with the desired precision
            approach_data = curve_data.iloc[:turning_point_index].copy()
            retract_data = curve_data.iloc[turning_point_index:].copy()
            approach_data['Zpiezo'] = approach_data['Zpiezo'].map('{:.18e}'.format)
            approach_data['Defl'] = approach_data['Defl'].map('{:.18e}'.format)
            approach_data['Force'] = approach_data['Force'].map('{:.18e}'.format)
            retract_data['Zpiezo'] = retract_data['Zpiezo'].map('{:.18e}'.format)
            retract_data['Defl'] = retract_data['Defl'].map('{:.18e}'.format)
            retract_data['Force'] = retract_data['Force'].map('{:.18e}'.format)

            # Save the approach and retract data with headers
            header_string = '# Zpiezo \t Defl \t Force\n'
            with open(atxtdir + str(i).zfill(3) + '_approach.txt', 'w') as f_approach:
                f_approach.write(header_string)
                approach_data.to_csv(f_approach, sep='\t', index=False, header=False)

            with open(rtxtdir + str(i).zfill(3) + '_retract.txt', 'w') as f_retract:
                f_retract.write(header_string)
                retract_data.to_csv(f_retract, sep='\t', index=False, header=False)

            # Optional: Plotting and saving the figure Broke for some reason??
            plt.figure(figsize=(10, 4))
            plt.plot(approach_data['Zpiezo'], approach_data['Force'], label='Approach', color='blue')
            plt.plot(retract_data['Zpiezo'], retract_data['Force'], label='Retract', color='red')
            plt.xlabel('Zpiezo (nm)')
            plt.ylabel('Force (nN)')
            plt.title(f'Curve {i} - Approach and Retract')
            plt.legend()
            plt.savefig(forcedir + str(i).zfill(3) + '_force_curve' + ext)
            plt.close()

        except Exception as e:
            print(f'Error in processing curve {i}: {e}')


# Usage example
# split_curves('path_to_your_data_file_without_extension')
