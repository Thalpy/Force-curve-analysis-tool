#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 22:53:17 2019

@author: Thalpy
Version: PhD 1.0
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

from .proc_force_sep import *
mpl.rcParams['lines.linewidth'] = 2


def plot_force_sep_res(name, force_curves, prefix = '', binsize = 1.0, ext = '.pdf', approach = True):

    if approach:
        outbase = name+"_force_curves"+os.sep+prefix+"approach_"
    else:
        outbase = name+"_force_curves"+os.sep+prefix+"retract_"

    zp_arr = force_curves['zp'].to_numpy()
    force_arr = force_curves['force'].to_numpy()
    df_arr = force_curves['df'].to_numpy()
    sep_arr = force_curves['sep'].to_numpy()

    k_c = force_arr[0]/df_arr[0]

    zp_bin, force_bin, force_std = binscatter(zp_arr, force_arr, binsize, minn=3)
    zp_bin, sep_bin, sep_std = binscatter(zp_arr, sep_arr, binsize, minn=3)

    zp_bin = zp_bin[3:-3]
    force_bin = force_bin[3:-3]
    force_std = force_std[3:-3]
    sep_std = sep_std[3:-3]
    sep_bin = sep_bin[3:-3]
    df_bin = force_bin/k_c
    df_std = force_std/k_c

    wz = (zp_bin <= 0)
    wnz = (zp_bin > 0)

    back = 0*wz+ zp_bin*wnz

    fig, ax = plt.subplots(figsize=(10,4))
    ax.errorbar(zp_bin,force_bin, yerr = force_std , fmt='o')
    # ax.plot(zp_bin_co,df_bin_co, 'go')
    # ax.plot(zp_bin_nc,df_bin_nc, 'ro')
    ax.plot(zp_bin, k_c*back, 'k-',zorder = 100)
    ax.set_ylabel('Force (nN)')
    ax.set_xlabel(r'$z_{piezo}$ (nm)')
    ax.set_xlim([-100,np.floor(np.max(zp_bin))])
    plt.tight_layout()
    fig.savefig(outbase+'force_zp_bin'+ext)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.errorbar(zp_bin,df_bin, yerr = df_std , fmt='+')
    ax.plot(zp_bin, back, 'k-')
    ax.set_ylabel('def (nm)')
    ax.set_xlabel(r'$z_{piezo}$ (nm)')
    ax.set_xlim([-100,np.floor(np.max(zp_bin))])
    ax.set_ylim([-5,np.floor(np.max(df_bin))])
    plt.tight_layout()
    fig.savefig(outbase+'df_zp_bin'+ext)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.errorbar(sep_bin, force_bin, yerr = force_std , fmt='-o')
    ax.axvline(x=0.0, color='k', zorder = 100)
    ax.set_ylabel('Force (nN)')
    ax.set_xlabel(r'$z_{sep}$ (nm)')
    ax.set_xlim([-1,30])
    plt.tight_layout()
    fig.savefig(outbase+'force_sep'+ext)
    plt.close(fig)

    fmax = np.ceil(np.max(force_bin))

    fig, ax = plt.subplots(figsize=(8,5))
    ax.errorbar(sep_bin,force_bin, yerr = force_std , fmt='-o')
    ax.axvline(x=0.0, color='k', zorder = 100)
    ax.set_ylabel('Force (nN)')
    ax.set_xlabel(r'$z_{sep}$ (nm)')
    ax.set_xlim([-1,50])
    ax.set_ylim([.01,fmax])
    ax.set_yscale("log")
    plt.tight_layout()
    fig.savefig(outbase+'force_sep_log_lin'+ext)
    plt.close(fig)

def plot_force_sep_c_nc(name, force_curves_c, force_curves_nc, prefix = '', binsize = 1.0,
    ext = '.pdf', approach = False, binsize_nc = 0.3):

    if approach:
        outbase = name+"_force_curves"+os.sep+prefix+"approach_"
    else:
        outbase = name+"_force_curves"+os.sep+prefix+"retract_"

    zp_c_arr = force_curves_c['zp'].to_numpy()
    force_c_arr = force_curves_c['force'].to_numpy()
    df_c_arr = force_curves_c['df'].to_numpy()
    sep_c_arr = force_curves_c['sep'].to_numpy()

    zp_nc_arr = force_curves_nc['zp'].to_numpy()
    force_nc_arr = force_curves_nc['force'].to_numpy()
    df_nc_arr = force_curves_nc['df'].to_numpy()
    sep_nc_arr = force_curves_nc['sep'].to_numpy()

    k_c = force_c_arr[0]/df_c_arr[0]

    zp_c_bin, force_c_bin, force_c_std = binscatter(zp_c_arr, force_c_arr, binsize, minn=3)
    zp_c_bin, sep_c_bin, sep_c_std = binscatter(zp_c_arr, sep_c_arr, binsize, minn=3)

    zp_c_bin = zp_c_bin[3:-3]
    force_c_bin = force_c_bin[3:-3]
    force_c_std = force_c_std[3:-3]
    sep_c_std = sep_c_std[3:-3]
    sep_c_bin = sep_c_bin[3:-3]

    sep_nc_arr[sep_nc_arr<0]=0
    sep_nc_bin, force_nc_bin, force_nc_std = binscatter(sep_nc_arr, force_nc_arr, binsize_nc, minn=3)
    df_nc_bin = force_nc_bin/k_c
    zp_nc_bin = df_nc_bin-sep_nc_bin

    #blending curves where zp overlaps
    zp_nc_max = np.max(zp_nc_bin)
    zp_c_min = np.min(zp_c_bin)

    

    fig, ax = plt.subplots(figsize=(10,4))
    ax.errorbar(zp_c_bin,force_c_bin, yerr = force_c_std , fmt='o')
    ax.errorbar(zp_nc_bin,force_nc_bin, yerr = force_nc_std, fmt='o', color='green')
    # ax.plot(zp_bin_nc,df_bin_nc, 'ro')
    ax.set_ylabel('Force (nN)')
    ax.set_xlabel(r'$z_{piezo}$ (nm)')
    ax.set_xlim([-100,np.floor(np.max(zp_c_bin))])
    plt.tight_layout()
    fig.savefig(outbase+'force_zp_c_c_bin'+ext)
    plt.close(fig)


    fig, ax = plt.subplots(figsize=(10,4))
    ax.errorbar(sep_nc_bin, force_nc_bin, yerr = force_nc_std , fmt='-o')
    ax.errorbar(sep_c_bin, force_c_bin, yerr = force_c_std , fmt='-o', color = 'red')
    ax.axvline(x=0.0, color='k', zorder = 100)
    ax.set_ylabel('Force (nN)')
    ax.set_xlabel(r'$z_{sep}$ (nm)')
    ax.set_xlim([-1,30])
    plt.tight_layout()
    fig.savefig(outbase+'force_sep_nc'+ext)
    plt.close(fig)
    #
    # fmax = np.ceil(np.max(force_bin))
    #
    # fig, ax = plt.subplots(figsize=(8,5))
    # ax.errorbar(sep_bin,force_bin, yerr = force_std , fmt='-o')
    # ax.axvline(x=0.0, color='k', zorder = 100)
    # ax.set_ylabel('Force (nN)')
    # ax.set_xlabel(r'$z_{sep}$ (nm)')
    # ax.set_xlim([-1,50])
    # ax.set_ylim([.01,fmax])
    # ax.set_yscale("log")
    # plt.tight_layout()
    # fig.savefig(outbase+'force_sep_log_lin'+ext)
    # plt.close(fig)

def get_contact_forces(name, force_curves, prefix = '',
                       thresh = 0.0, binsize = 1.0, ext = '.pdf', out = True, 
                       smooth_win = 15, approach = True):

    # Adjust the base path for retract curves
    curve_type = "retract" if not approach else "approach"
    outbase = f"{name}_force_curves{os.sep}{prefix}{curve_type}_"
    outdir = outbase + 'force_curves' + os.sep
    contdir = outdir + 'contact' + os.sep

    # Read curve parameters
    curve_param = pd.read_csv(outbase + 'curve_param.csv')
    curve_param['F_cont'] = np.nan

    if out:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        if not os.path.exists(contdir):
            os.makedirs(contdir)

    # Process force curve data
    zp_arr = force_curves['zp'].to_numpy()
    force_arr = force_curves['force'].to_numpy()
    df_arr = force_curves['df'].to_numpy()
    sep_arr = force_curves['sep'].to_numpy()

    k_c = force_arr[0] / df_arr[0]

    # Binning the data
    zp_bin, force_bin, force_std = binscatter(zp_arr, force_arr, binsize, minn=3)
    zp_bin, sep_bin, sep_std = binscatter(zp_arr, sep_arr, binsize, minn=3)

    # Trimming the data
    zp_bin = zp_bin[3:-3]
    force_bin = force_bin[3:-3]
    force_std = force_std[3:-3]
    sep_std = sep_std[3:-3]
    sep_bin = sep_bin[3:-3]
    df_bin = force_bin / k_c
    df_std = force_std / k_c

    # Identifying the threshold crossing point
    cross = np.where(np.diff(np.sign(sep_bin - thresh)))[0]
    cross = np.min(cross)

    # Linear interpolation to find contact point
    x0 = sep_bin[cross + 1]
    x1 = sep_bin[cross]
    y0 = force_bin[cross + 1]
    y1 = force_bin[cross]
    f_c = (y0 * x1 - y1 * x0) / (x1 - x0)

    # Plotting the overall force-separation curve
    fmax = np.ceil(np.max(force_bin))
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(sep_bin, force_bin, yerr=force_std, fmt='-o')
    ax.plot(0, f_c, 'r^', zorder=100)
    ax.axvline(x=0.0, color='k', zorder=100)
    ax.set_ylabel('Force (nN)')
    ax.set_xlabel(r'$z_{sep}$ (nm)')
    ax.set_xlim([-1, 50])
    ax.set_ylim([0.01, fmax])
    ax.set_yscale("log")
    plt.tight_layout()
    fig.savefig(outbase + 'force_sep_log_lin_cont' + ext)
    plt.close(fig)

    # Processing individual force curves
    curve_param = curve_param.set_index('index')
    for i in force_curves['index'].unique():
        print('finding contact: ' + str(i))
        dfi = force_curves.groupby(['index']).get_group(i)
        force = dfi['force'].to_numpy()
        sep = dfi['sep'].to_numpy()

        # Applying smoothing filter
        sep = savgol_filter(sep, smooth_win, 3)
        force = savgol_filter(force, smooth_win, 3)

        # Finding contact point
        cross = np.where(np.diff(np.sign(sep - thresh)))[0]
        cross = np.min(cross)

        # Linear interpolation for contact force
        x0 = sep[cross + 1]
        x1 = sep[cross]
        y0 = force[cross + 1]
        y1 = force[cross]
        fci = (y0 * x1 - y1 * x0) / (x1 - x0)

        # Plotting individual force-separation curves
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(sep, force, '-o')
        ax.plot(0, fci, 'r^', zorder=100)
        ax.axvline(x=0.0, color='k', zorder=100)
        ax.set_ylabel('Force (nN)')
        ax.set_xlabel(r'$z_{sep}$ (nm)')
        ax.set_xlim([-1, 50])
        ax.set_ylim([0.01, fmax])
        ax.set_yscale("log")
        plt.tight_layout()
        fig.savefig(contdir + str(i).zfill(3) + '_force_sep_cont' + ext)
        plt.close(fig)

        curve_param.loc[i, 'F_cont'] = fci

    # Updating curve parameters
    curve_param.to_csv(outbase + 'curve_param.csv')

    # Plotting the histogram of contact forces
    fig, ax = plt.subplots(figsize=(10,4))
    f_c_arr = curve_param['F_cont'].to_numpy()
    f_c_arr = f_c_arr[~np.isnan(f_c_arr)]
    ax.hist(f_c_arr, bins=50)
    mean_f_c = np.mean(f_c_arr)
    f_c_std = np.std(f_c_arr)
    ax.set_ylabel('Counts')
    ax.set_xlabel(r'Contact Force (nN)')
    ax.set_title('Contact Force = %.3f +/- %.3f' % (mean_f_c, f_c_std))
    plt.tight_layout()
    fig.savefig(outbase + 'f_c_hist' + ext)
    plt.close(fig)

    return f_c, f_c_std

def plot_force_h_scaling(name, force_curves, prefix = '',
                         thresh = 0.0, binsize = 1.0, ext = '.pdf', out = True, approach = True):

    # Adjust the base path for retract curves
    curve_type = "retract" if not approach else "approach"
    outbase = f"{name}_force_curves{os.sep}{prefix}{curve_type}_"
    outdir = outbase + 'force_curves' + os.sep
    contdir = outdir + 'contact' + os.sep

    # Read curve parameters
    curve_param = pd.read_csv(outbase + 'curve_param.csv')
    if out:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        if not os.path.exists(contdir):
            os.makedirs(contdir)

    # Process force curve data
    zp_arr = force_curves['zp'].to_numpy()
    force_arr = force_curves['force'].to_numpy()
    df_arr = force_curves['df'].to_numpy()
    sep_arr = force_curves['sep'].to_numpy()

    k_c = force_arr[0] / df_arr[0]

    # Binning the data
    zp_bin, force_bin, force_std = binscatter(zp_arr, force_arr, binsize, minn=3)
    zp_bin, sep_bin, sep_std = binscatter(zp_arr, sep_arr, binsize, minn=3)

    # Trimming the data
    zp_bin = zp_bin[3:-3]
    force_bin = force_bin[3:-3]
    force_std = force_std[3:-3]
    sep_std = sep_std[3:-3]
    sep_bin = sep_bin[3:-3]
    df_bin = force_bin / k_c
    df_std = force_std / k_c

    # Identifying the threshold crossing point
    cross = np.where(np.diff(np.sign(sep_bin - thresh)))[0]
    cross = np.min(cross)

    # Taking only the relevant part of the force curve
    # For retract curves, this might be adjusted based on specific analysis requirements
    force_bin = force_bin[0:cross]
    force_std = force_std[0:cross]
    sep_std = sep_std[0:cross]
    sep_bin = sep_bin[0:cross]

    # Scaling the force
    force_h_scaled = force_bin * sep_bin
    force_h_std = np.sqrt((sep_bin**2) * (force_std**2) + (force_bin**2) * (sep_std**2))

    # Plotting
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(sep_bin, force_h_scaled, yerr=force_std, fmt='-o')
    ax.set_ylabel('Fh (nN*nm)')
    ax.set_xlabel(r'$h$ (nm)')
    ax.set_xlim([0.1, 100])
    ax.set_yscale("log")
    ax.set_xscale("log")
    plt.tight_layout()
    fig.savefig(outbase + 'force_h' + ext)
    plt.close(fig)


def profile_attractive_forces(name, force_curves, prefix = '', far_thresh = 60,
    thresh = 0.0, binsize = 1.0, ext = '.pdf', out = True, smooth_win=15, approach = False):

    if approach:
        outbase = name+"_force_curves"+os.sep+prefix+"approach_"
    else:
        outbase = name+"_force_curves"+os.sep+prefix+"retract_"

    outdir = outbase+'force_curves'+os.sep
    attrdir = outdir+'attractive'+os.sep

    curve_param = pd.read_csv(outbase+'curve_param.csv')
    curve_param['F_attr'] = np.nan
    curve_param['attr_sep'] = np.nan
    curve_param['F_cont'] = np.nan

    #create empty dataframes for contact and non-contact forces
    my_cols =['zp', 'df', 'sep', 'force', 'index']
    force_curves_c = pd.DataFrame(columns=my_cols)
    force_curves_nc = pd.DataFrame(columns=my_cols)

    if out:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        if not os.path.exists(attrdir):
            os.makedirs(attrdir)

    force_arr = force_curves['force'].to_numpy()
    df_arr = force_curves['df'].to_numpy()

    k_c = force_arr[0]/df_arr[0]

    curve_param = curve_param.set_index('index')
    #force_curves = force_curves.groupby(['index'])

    for i in force_curves['index'].unique() :
        try:
            print('finding contact: '+str(i))
            dfi = force_curves.groupby(['index']).get_group(i)
            force = dfi['force'].to_numpy()
            sep = dfi['sep'].to_numpy()
            zp = dfi['zp'].to_numpy()
            df = dfi['df'].to_numpy()


            sep_s = savgol_filter(sep, smooth_win, 3)
            force_s = savgol_filter(force, smooth_win, 3)
            zp_s = savgol_filter(zp, 2*smooth_win+1, 3)
            df_s = savgol_filter(df, 2*smooth_win+1, 3)

            deriv_df = np.gradient(df_s, zp_s)

            cross= np.where(np.diff(np.sign(sep_s-thresh)))[0]
            cross = np.min(cross)

            #linear interp to get Contact (sep =0)
            x0=sep_s[cross+1]
            x1=sep_s[cross]
            y0 = force_s[cross+1]
            y1 = force_s[cross]
            fci = (y0*x1 - y1*x0)/(x1-x0)


            far= np.where(np.diff(np.sign(sep-far_thresh)))[0]
            far = np.min(far)

            f_min = np.min(force[far:cross])
            min_loc = np.argmin(force[far:cross])+far
            sep_min = sep[min_loc] #separation at F_min
            if (sep_min > 30):
                f_min = 0  # if min is too far away from contact point, it's just noise, not attraction

            #find points right before/after any jump to/from contact
            #data here captures transient motion of cantelever, not ture force profile
            bad = (deriv_df <= -1.25)

            bad[min_loc::] = False
            if np.sum(bad) >2:
                badrange = np.where(bad)[0]
                bright = badrange[-1]
                if len(badrange[badrange>far])>0:
                    bleft = badrange[badrange>far][0]
                else:
                    bleft = bright
                #print(far, bleft, bright)
                bad[bleft:bright] = True

            sep_g = sep[~bad]
            force_g = force[~bad]
            df_g = df[~bad]
            zp_g = zp[~bad]


            fig, ax = plt.subplots(figsize=(8,5))
            ax.plot(sep_g,force_g, '-o')
            ax.plot(sep[cross:-1],force[cross:-1], 'o', color='purple')
            ax.plot(0,fci, 'r^', zorder =100)
            ax.plot(sep_min,f_min, 'gs', zorder =100)
            ax.axvline(x=0.0, color='k', zorder = 100)
            ax.set_ylabel('Force (nN)')
            ax.set_xlabel(r'$z_{sep}$ (nm)')
            ax.set_xlim([-1,50])

            plt.tight_layout()
            fig.savefig(attrdir+str(i).zfill(3)+'_force_sep_attr'+ext)
            plt.close(fig)

            curve_param.loc[i,'F_cont']=fci
            curve_param.loc[i,'F_attr']=f_min
            curve_param.loc[i,'attr_sep']=sep_min
            #print(cross, len(zp_g[cross:-1]), len(zp_g[0:cross]))

            dfi_c = pd.DataFrame({
                'zp': zp_g[cross:-1],
                'df': df_g[cross:-1],
                'sep': sep_g[cross:-1],
                'force': force_g[cross:-1],
                'index': i
            })

            dfi_nc = pd.DataFrame({
                'zp': zp_g[0:cross],
                'df': df_g[0:cross],
                'sep': sep_g[0:cross],
                'force': force_g[0:cross],
                'index': i
            })

            force_curves_c=force_curves_c.append(dfi_c, sort=True)
            force_curves_nc=force_curves_nc.append(dfi_nc, sort=True)
        except:
            print('error with curve: '+str(i)+" skipping! ")

    curve_param.to_csv(outbase+'curve_param.csv')
    force_curves_c.to_csv(outbase+'force_curves_c.csv')
    force_curves_nc.to_csv(outbase+'force_curves_nc.csv')


    fig, ax = plt.subplots(figsize=(10,4))
    f_c_arr = curve_param['F_cont'].to_numpy()
    f_c_arr = f_c_arr[~np.isnan(f_c_arr)]
    ax.hist(f_c_arr, bins=50)
    mean_f_c = np.mean(f_c_arr)
    f_c_std = np.std(f_c_arr)
    ax.set_ylabel('Counts')
    ax.set_xlabel(r'Contact Force (nN)')
    ax.set_title('Contact Force = %.3f +/- %.3f' % (mean_f_c, f_c_std ) )

    plt.tight_layout()
    fig.savefig(outbase+'f_c_hist'+ext)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10,4))
    f_a_arr = curve_param['F_attr'].to_numpy()
    f_a_arr = f_a_arr[~np.isnan(f_a_arr)]
    ax.hist(f_a_arr, bins=50)
    mean_f_a = np.mean(f_a_arr)
    f_a_std = np.std(f_a_arr)
    ax.set_ylabel('Counts')
    ax.set_xlabel(r'Pull Off Force (nN)')
    ax.set_title(r'$F_{attr}$ = %.3f +/- %.3f' % (mean_f_a, f_a_std ) )
    plt.tight_layout()
    fig.savefig(outbase+'f_a_hist'+ext)
    plt.close(fig)



    return force_curves_c, force_curves_nc
