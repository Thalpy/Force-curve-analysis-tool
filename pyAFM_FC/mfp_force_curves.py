# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 2017

@author: Joe
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
    #split up a bunch of force curves


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
