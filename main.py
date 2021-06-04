#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 14:20:15 2021

@author: arifeinberg
"""

import os
os.chdir('/Users/arifeinberg/target2/fs03/d0/arifein/python/pythonHgBenchmark')

from load_Hgmodel_data import open_Hg_spc
from TGMAndObs import SurfaceObsTGM

fn_old = '../../GEOS-Chem_runs/run0003/OutputDir/GEOSChem.SpeciesConc.alltime_m.nc4'
fn_new = '../../GEOS-Chem_runs/run0006/OutputDir/GEOSChem.SpeciesConc.alltime_m.nc4'

ds1, ds2 = open_Hg_spc(fn_old, fn_new)
#%%
SurfaceObsTGM(ds1, ds2, 2015)
