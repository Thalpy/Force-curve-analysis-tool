# -*- coding: utf-8 -*-
"""
Created on Sun Jan  21 2018

@author: Thalpy
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
import gc


import matplotlib  as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import quiver, cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import csv
mpl.rcParams['lines.linewidth'] = 2
plt.ioff()

from .utils import *

def proc_force_sep(name, k_c = 0.188, binsize = 1.0, dfit_win = 50, dfit_off = 60,
    cfit_min = 30, cfit_max = 90, fitbin = 5, cthresh = 45, ext = '.pdf',
    clear = True, out=True, prefix = '', min_cont_pts=5, slope = 1.38, approach =True):

    #read in a bunch of force curves
    if approach:
        in_txtdir = name+"_force_curves"+os.sep+"approach"
    else:
        in_txtdir = name+"_force_curves"+os.sep+"retract"

    namelist = sorted(glob.glob(in_txtdir+os.sep+'*.txt'))

    txtdir = name+"_force_curves"+os.sep+prefix
    if not os.path.exists(txtdir):
        os.makedirs(txtdir)

    if approach:
        outbase = txtdir +"approach_"
    else:
        outbase = txtdir +"retract_"

    outdir = outbase+'force_curves'+os.sep
    derivdir = outdir+'deriv'+os.sep
    backdir = outdir+'back'+os.sep
    defldir = outdir+'defl'+os.sep

    if clear and os.path.exists(outdir):
        shutil.rmtree(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if out:
        if not os.path.exists(derivdir):
            os.makedirs(derivdir)
        if not os.path.exists(backdir):
            os.makedirs(backdir)
        if not os.path.exists(defldir):
            os.makedirs(defldir)

    zp_arr = []
    df_arr = []
    zp_bin_arr = []
    df_deriv_arr = []
    zp_deriv_used = []
    df_deriv_used = []

    #create empty dataframes, one for data, one for parameters
    my_cols =['zp', 'df', 'sep', 'force', 'index']
    force_curves = pd.DataFrame(columns=my_cols)

    para_cols = ['index','slope', 'k_c', 'F_max', 'F_min', 'minLoc', 'z_c', 'z_c1' ]
    curve_param = pd.DataFrame(columns=para_cols)

    slope_arr = np.array([])
    zc_arr = np.array([])

    plt.ioff()

    fig, ax = plt.subplots(figsize=(10,4))
    i=0


    for fullname in namelist:
        try:
            #the files in 'namelist' have names structured like '005_approach.txt'
            #this parses the name and pulls out '5' as the index
            #would need to change if files names are different
            numbCountJ = [int(s) for s in str(fullname) if s.isdigit()]
            numbCountJ = ("".join(map(str, numbCountJ[-3:])))

            da = np.loadtxt(fullname)
            zp = da[:,0]*1000.0 #convert microns to nm
            df = da[:,1]*1.0e9 #deflection in nm
            if not approach:
                zp = zp[::-1]
                df = df[::-1]

            avgdf_far = np.mean(df[0:(fitbin*10)])
            df = df - avgdf_far

            if np.max(df) < cthresh:
                print('skipping ' +str(numbCountJ).zfill(3))
                continue

            # CHANGED
            #guess contact point and work backwards
            #cthresh = 'first pass' guess deflection threshold for contact
            #returns de-drifted zp, df.  data used to de-drift runs from dfit_win:2*dfit_win
            zp, df, dfit_used = remove_farfield_drift(zp,df,fitbin,cfit_min,cthresh, dfit_off,
                dfit_win, out, backdir, numbCountJ, ext)

            f_max = k_c*np.max(df)
            #binning stuff - helps with fits to find contact point
            zp, df, zpb, dfb, dfb_std = bin_z_df(zp, df, fitbin)

            # fit region where deflection is between cfit_min and cfit_max to straight line
            # if there aren't enough points in here (>5 seems to work), skip this curve
            # note that this is sensitive to fitbin -less binning will naturally give
            # you more (noisy) points
            w = np.where( (dfb > cfit_min) & (dfb <cfit_max))
            if np.size(w) < min_cont_pts:
                print('skipping ' +str(numbCountJ).zfill(3))
                continue

            print('aligning '+str(numbCountJ).zfill(3))

            zp, df, sep, zpb, sepb, deriv_df, zc, zc1 = align_contact_reg(zpb, dfb, dfb_std, zp,
                df, w, dfit_used, cfit_max, out, defldir, derivdir, numbCountJ, ext, slope, cthresh)

            zp_arr = np.concatenate((zp_arr, zp), axis=0)
            df_arr = np.concatenate((df_arr, df), axis=0)
            zp_bin_arr = np.concatenate((zp_bin_arr, zpb), axis=0)
            df_deriv_arr = np.concatenate((df_deriv_arr, deriv_df), axis=0)
            zp_deriv_used = np.concatenate( (zp_deriv_used , zpb[w]), axis=0)
            df_deriv_used = np.concatenate( (df_deriv_used ,deriv_df[w]), axis=0)

            zc_arr = np.append(zc_arr, [zc], axis=0)

            dfi = pd.DataFrame({
                'zp': zp,
                'df': df,
                'sep': sep,
                'force': k_c*df,
                'index': int(numbCountJ)
            })

            f_min = np.min(k_c*dfb)
            min_loc = sepb[np.argmin(dfb)]
            f_noise = np.std(k_c*dfb[0:np.argmin(dfb)])

            dfp = pd.DataFrame({
                'index': [int(numbCountJ)],
                'slope': [slope],
                'k_c': [k_c],
                'F_max': [f_max],
                'z_c': [zc],
                'z_c1': [zc1],
                'F_min': [f_min],
                'minLoc': [min_loc],
                'F_noise': [f_noise]
            })

            force_curves=force_curves.append(dfi, sort=True)
            curve_param=curve_param.append(dfp, sort=True)
            #del dfi

            ax.plot(zp,df , '-')
            ax.set_ylabel('def (nm)')
            ax.set_xlabel('zp (nm)')
            plt.tight_layout()
            i=i+1
        except:
            print('error in curve ' +str(numbCountJ).zfill(3))
            print("failed to process curve!")
            pass

    print('%d curves processed out of %d total force curves' % (i, len(namelist)))

    sep_all = df_arr - zp_arr
    force_all = df_arr*k_c

    zp_bin, df_bin, df_std = binscatter(zp_arr, df_arr, binsize, minn=3)
    zpb_bin, df_deriv_bin, df_deriv_std = binscatter(zp_bin_arr, df_deriv_arr, binsize, minn=2)

    #ends of binned data are biased by how data is aligned (see curve overlays)
    #cut this off
    zp_bin = zp_bin[3:-3]
    df_bin = df_bin[3:-3]
    df_std = df_std[3:-3]

    sep = df_bin - zp_bin
    force = df_bin*k_c
    force_std = df_std*k_c


    wz = (zp_bin <= 0)
    wnz = (zp_bin > 0)

    back = 0*wz+ zp_bin*wnz

    ax.plot(zp_bin, back, 'k-')
    fig.savefig(outdir+'df_zp'+ext)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.errorbar(zpb_bin,df_deriv_bin, yerr = df_deriv_std , fmt='b+')
    ax.plot(zp_deriv_used, df_deriv_used, 'o', color = 'purple', alpha = 0.6  )
    ax.set_ylabel('defection slope')
    ax.set_xlabel(r'$z_{piezo}$ (nm)')
    ax.set_xlim([-100,np.ceil(cfit_max)])
    plt.tight_layout()
    fig.savefig(outdir+'df_deriv_bin'+ext)
    plt.close(fig)


    fig, (ax1, ax2) = plt.subplots(nrows =1, ncols =2, figsize=(8,5))
    ax1.hist(slope_arr, bins =20)
    ax1.set_xlabel('Contact slope (nm/nm)')
    ax2.hist(zc_arr, bins =20)
    ax2.set_xlabel(r'Contact point $z_c$ (nm)')
    fig.savefig(outdir+'cfit_hists'+ext)
    plt.close(fig)

    np.savetxt(outbase+'_force_sep.txt', np.c_[sep,force,force_std])
    np.savetxt(outdir+'deriv.txt', np.c_[zpb_bin,df_deriv_bin,df_deriv_std])

    f = open(outbase+'parameter_log.txt','w')
    f.write('Timestamp: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + '\n \n')

    f.write('name = \'%s\'  \n \n' % name )
    f.write('proc_force_sep(name, k_c = %.3f, binsize = %.1f, dfit_win = %d, dfit_off = %d, cfit_min = %d, cfit_max = %d, fitbin = %d, cthresh = %d, ext = \'%s\', clear = %s, out=%s) \n \n'
        % ( k_c, binsize, dfit_win, dfit_off, cfit_min, cfit_max, fitbin, cthresh, ext, str(clear), str(out)))
    f.write('%d curves processed out of %d total force curves \n' % (i, len(namelist)))
    #f.write('Spring const. kc = %.3f \n' % k_c )
    #f.write('binsize = %.2f \n' % binsize )
    #f.write('fitbin = %d \n' % fitbin )
    #f.write('dfit_win = %d , dfit_off = %d \n' % (dfit_win, dfit_off) )
    #f.write('cfit_min = %d , cfit_max = %d \n' % (cfit_min, cfit_max) )
    f.close()
    plt.ion()

    curve_param = curve_param.set_index('index')
    force_curves.to_csv(outbase+'force_curves.csv')
    curve_param.to_csv(outbase+'curve_param.csv')

    return force_curves, curve_param

def comp_def_deriv(name, k_c = 0.188, binsize = 1.0, dfit_win = 50, dfit_off = 60,
    cfit_min = 60, cfit_max = 120, fitbin = 5, cthresh = 45, ext = '.pdf',
    clear = True, out=True, prefix = '', approach = True):

    #read in a bunch of force curves
    if approach:
    	in_txtdir = name+"_force_curves"+os.sep+"approach"
    else:
        in_txtdir = name+"_force_curves"+os.sep+"retract"

    namelist = sorted(glob.glob(in_txtdir+os.sep+'*.txt'))

    txtdir = name+"_force_curves"+os.sep+prefix
    if not os.path.exists(txtdir):
        os.makedirs(txtdir)

    if approach:
        outbase = txtdir +"approach_"
    else:
        outbase = txtdir +"retract_"

    outdir = outbase+'force_curves'+os.sep
    derivdir = outdir+'deriv'+os.sep
    backdir = outdir+'back'+os.sep


    if clear and os.path.exists(outdir):
        shutil.rmtree(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if out:
        if not os.path.exists(derivdir):
            os.makedirs(derivdir)
        if not os.path.exists(backdir):
            os.makedirs(backdir)

    #create empty dataframe
    my_cols =['zp_bin', 'df_bin', 'deriv_df', 'index']
    deriv_curves = pd.DataFrame(columns=my_cols)

    slope_arr = np.array([])
    zc_arr = np.array([])

    plt.ioff()

    fig, ax = plt.subplots(figsize=(10,4))
    i=0
    for fullname in namelist:
        #the files in 'namelist' have names structured like '005_approach.txt'
        #this parses the name and pulls out '5' as the index
        #would need to change if files names are different
        numbCountJ = [int(s) for s in str(fullname) if s.isdigit()]
        numbCountJ = ("".join(map(str, numbCountJ[-3:])))

        da = np.loadtxt(fullname)
        zp = da[:,0]*1000.0 #convert microns to nm
        df = da[:,1]*1.0e9 #deflection in nm
        if not approach:
            zp = zp[::-1]
            df = df[::-1]

        avgdf_far = np.mean(df[0:(fitbin*10)])
        df = df - avgdf_far

        if np.max(df) <= cthresh:
            print('skipping: ' +str(numbCountJ).zfill(3))
            continue

        # CHANGED
        #guess contact point and work backwards
        #cthresh = 'first pass' guess deflection threshold for contact
        #returns de-drifted zp, df.  data used to de-drift runs from dfit_win:2*dfit_win
        zp, df, dfit_used = remove_farfield_drift(zp,df,fitbin,cfit_min,cthresh, dfit_off,
            dfit_win, out, backdir, numbCountJ, ext)


        #binning stuff - helps with computing derivatives
        zp, df, zpb, dfb, dfb_std = bin_z_df(zp, df, fitbin)

        print('comp. deriv: ' +str(numbCountJ).zfill(3))

        #deriv_df = np.gradient(dfb, zpb)
        deriv_df, zpb, dfb = diff_central(zpb, dfb)

        fig1, ax1 = plt.subplots(figsize=(10,4))
        ax1.plot(zpb,deriv_df , '+-')
        #ax1.plot(zpb[w],deriv_df[w] , 'o', color = 'purple')
        ax1.set_ylabel('deflection slope')
        ax1.set_xlabel(r'$z_p$ (nm)')
        #ax1.set_xlim([-200,np.ceil(cfit_max+30)])
        plt.tight_layout()
        fig1.savefig(derivdir+str(numbCountJ).zfill(3)+'_slope_zc'+ext)
        plt.close(fig1)

        dfi = pd.DataFrame({
            'zp_bin': zpb,
            'df_bin': dfb,
            'deriv_df': deriv_df,
            'index': numbCountJ
        })

        deriv_curves=deriv_curves.append(dfi, sort=True)
        #del dfi

        ax.plot(dfb,deriv_df , '-')
        ax.set_ylabel('deriv def')
        ax.set_xlabel('deflection (nm)')
        plt.tight_layout()
        i=i+1

    print('%d curves processed out of %d total force curves' % (i, len(namelist)))

    df_all = deriv_curves['df_bin'].to_numpy()
    deriv_df_all = deriv_curves['deriv_df'].to_numpy()

    df_bin, deriv_df_bin, deriv_df_std =binscatter(df_all, deriv_df_all, binsize, minn=2)

    #ends of binned data are biased by how data is aligned (see curve overlays)
    #cut this off
    df_bin = df_bin[2:-2]
    deriv_df_bin = deriv_df_bin[2:-2]
    deriv_df_std = deriv_df_std[2:-2]


    fig.savefig(outdir+'deriv_df_df'+ext)
    plt.close(fig)

    w = (df_all>cfit_min) & (df_all<cfit_max)
    fig, ax = plt.subplots(figsize=(10,4))
    ax.hist(deriv_df_all[w], bins=50)
    meanslope = np.mean(deriv_df_all[w])
    slope_std = np.std(deriv_df_all[w])/np.sqrt(len(w))  #stdev of mean
    ax.set_ylabel('Counts')
    ax.set_xlabel(r'deflection slope')
    ax.set_title('Deflection slope = %.3f +/- %.3f' % (meanslope, slope_std ) )
    #ax.set_xlim([-1,np.ceil(cfit_max)])
    #ax.set_ylim([-5,np.ceil(cfit_max)])
    plt.tight_layout()
    fig.savefig(outbase+'deriv_hist'+ext)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.errorbar(df_bin,deriv_df_bin, yerr = deriv_df_std , fmt='+')
    ax.axhline(meanslope, color='k')
    ax.axvline(cfit_min, color='k', linestyle =':')
    ax.axvline(cfit_max, color='k', linestyle =':')
    ax.set_ylabel('deriv def ')
    ax.set_xlabel(r'deflection (nm)')
    #ax.set_xlim([-1,np.ceil(cfit_max)])
    #ax.set_ylim([-5,np.ceil(cfit_max)])
    plt.tight_layout()
    fig.savefig(outdir+'deriv_df_df_bin'+ext)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.errorbar(df_bin,deriv_df_bin, yerr = deriv_df_std , fmt='+')
    ax.axhline(meanslope, color='k')
    ax.axvline(cfit_min, color='k', linestyle =':')
    ax.axvline(cfit_max, color='k', linestyle =':')
    ax.set_ylabel('deriv def ')
    ax.set_xlabel(r'deflection (nm)')
    ax.set_xlim([-1,np.ceil(2.0*cfit_max)])
    #ax.set_ylim([-5,np.ceil(cfit_max)])
    plt.tight_layout()
    fig.savefig(outbase+'deriv_df_df_bin_zoom'+ext)
    plt.close(fig)


    deriv_curves.to_csv(outbase+'deriv_curves.csv')

    return meanslope


def remove_farfield_drift(zp,df,fitbin,cfit_min,cthresh, dfit_off,
    dfit_win, out, backdir, numbCountJ , ext):

    #guess for contact point (array index)
    contguess = np.where(np.diff(np.sign(df-cthresh)))[0]
    #dfit_win and dfit_off in real units (nm)
    #convert to array index
    nm_per_pt = np.abs(zp[-1]-zp[0])/np.size(zp)
    dfit_win = int(dfit_win/nm_per_pt)
    dfit_off = int(dfit_off/nm_per_pt)
    dfit_max = contguess[0] - dfit_off
    dfit_min = dfit_max - dfit_win

    #remove background drift and offset from deflection, using points far away from contact
    pb = np.polyfit(zp[dfit_min:dfit_max], df[dfit_min:dfit_max], 1)
    backwidth = dfit_max-dfit_min
    blow = np.floor(backwidth/4).astype('i8')
    bhigh = np.floor(backwidth/2).astype('i8')
    zz = np.linspace(np.floor(zp[dfit_min-blow]),np.ceil(zp[dfit_max+bhigh]), num=100)
    if out:
        fig1, ax1 = plt.subplots(figsize=(10,4))
        ax1.plot(zp[dfit_min-blow:dfit_max+bhigh],df[dfit_min-blow:dfit_max+bhigh] , '+-')
        ax1.plot(zp[dfit_min:dfit_max],df[dfit_min:dfit_max] , 'o', color = 'orange')
        ax1.plot(zz,pb[1]+pb[0]*zz, 'k-')
        ax1.set_ylabel('raw deflection (nm)')
        ax1.set_xlabel(r'$z_p$ (nm)')
        plt.tight_layout()
        fig1.savefig(backdir+str(numbCountJ).zfill(3)+'_back'+ext)
        plt.close(fig1)

    df = df - (pb[1]+pb[0]*zp)
    df = df[dfit_min- dfit_win::]
    zp = zp[dfit_min- dfit_win::]

    return zp, df, dfit_win

def bin_z_df(zp, df, fitbin):
    npts = int(zp.shape[0]/fitbin)
    df = df[0:npts*fitbin]
    zp = zp[0:npts*fitbin]
    zpb = zp.reshape(npts, -1).mean(axis=1)
    dfb = df.reshape(npts, -1).mean(axis=1)
    dfb_std = df.reshape(npts, -1).std(axis=1)
    return zp, df, zpb, dfb, dfb_std


def align_contact_reg(zpb, dfb, dfb_std, zp, df, w, dfit_win, cfit_max, out,
    defldir, derivdir, numbCountJ, ext, slope, cthresh):

    zpcfit = zpb[w]
    dfcfit = dfb[w]

    intercept = np.mean(dfcfit-slope*zpcfit)

    # extrapolate back to where surfaces would contact (without finite ranged forces)
    zc = (-1)*intercept/slope

    #should have pa[0]=1 if deflection calibrated. adjustment here will correct for run to run drift
    df = df/slope
    dfb = dfb/slope
    dfcfit = dfcfit/slope

    deriv_df = np.gradient(dfb, zpb)
    #deriv_df, zpb, dfb = diff_central(zpb, dfb)
    zpb = zpb-zc

    #in case of jump to contact, take the jump point instead of the interpolated zc
    if np.min(deriv_df) <= -1:
        minloc = np.argmin(df)
        dfcut = df[0:minloc+1]
        jump = np.max(np.where(dfcut > (np.min(dfcut)+5)))
        zc1 = zp[jump]
        #print(zc)
    else:
        zc1 = zc

    # cut off high forces/deflections beyond cfit_max
    # this avoids the region where the photodiode is non-linear, giving spurious backwards bending
    keep = (df < cfit_max)
    zp = zp[keep]
    df = df[keep]

    # shift zp relative to extrapolated contact point
    zp = zp-zc
    zpcfit = zpcfit - zc

    zz = np.linspace(-200,np.ceil(cfit_max+20), num=20000)
    back = zz*(zz >0)-(zc-zc1)*(zz >0)
    zz = zz-(zc-zc1)

    if zc1 ==zc:
        sep = df-zp
        sepb = dfb-zpb
    else:
    # find where the undeflected tip is not in contact
        wz = (zp <= 0)
        wnz = (zp > 0)
        sep = df-zp
        #sep[wnz] = df[wnz]-zp[wnz]+(zc-zc1)
        wz = (zpb <= 0)
        wnz = (zpb > 0)
        sepb = dfb-zpb
        #sepb[wnz] = dfb[wnz]-zpb[wnz]+(zc-zc1)

    if out:
        fig1, ax1 = plt.subplots(figsize=(10,4))
        ax1.plot(zp, df , '+-')
        ax1.plot(zpcfit, dfcfit , 'o', color = 'purple')
        ax1.plot(zp[dfit_win:2*dfit_win], df[dfit_win:2*dfit_win] , 'o', color = 'orange')
        ax1.plot(zz, back, 'k-')

        # Find the z-position where the deflection first exceeds cthresh
        cthresh_crossing = np.argmax(df > cthresh)  # This assumes df is sorted in ascending z-position
        if cthresh_crossing.size > 0 and df[cthresh_crossing] > cthresh:
            zp_cthresh = zp[cthresh_crossing]
            ax1.plot(zp_cthresh, cthresh, 'rx', markersize=10)  # Plot a red 'X' at the threshold crossing

        ax1.set_ylabel('deflection (nm)')
        ax1.set_xlabel(r'$z_p-z_c$ (nm)')
        ax1.set_xlim([-200, np.ceil(cfit_max/slope+10)])
        plt.tight_layout()
        fig1.savefig(defldir+str(numbCountJ).zfill(3)+'_def_zc'+ext)
        plt.close(fig1)

    return zp, df, sep, zpb, sepb, deriv_df, zc, zc1
