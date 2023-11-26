# -*- coding: utf-8 -*-
"""
Created on Thur 28 March 2019

Move some utility functions for binning data, computing derivatives etc into a
separate utility program

@author: Joe
Version PhD 1.0
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

import csv

def binscatter(x,y, binsize, minn = 3):
    xmin = np.min(x)
    x = x-xmin
    edges = np.arange(0, np.ceil(np.max(x)/binsize)+1)*binsize;
    n,b = np.histogram(x, bins=edges)
    sy, b = np.histogram(x, bins=edges, weights=y)
    sy2, b = np.histogram(x, bins=edges, weights=np.power(y,2))
    w = np.where(n>minn)
    mean = sy[w] / n[w]
    std = np.sqrt(sy2[w]/n[w] - mean*mean)
    xbin = (b[1:]+b[:-1])/2
    xbin = xbin[w] + xmin
    return xbin, mean, std
    #return xbin, mean

def diff_central(x, y):
    x0 = x[:-2]
    x1 = x[1:-1]
    x2 = x[2:]
    y0 = y[:-2]
    y1 = y[1:-1]
    y2 = y[2:]
    f = (x2 - x1)/(x2 - x0)

    return (1-f)*(y2 - y1)/(x2 - x1) + f*(y1 - y0)/(x1 - x0), x1, y1
