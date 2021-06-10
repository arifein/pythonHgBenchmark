#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 14:20:15 2021
Main file for running Hg benchmark in python and producing plots
@author: arifeinberg
"""

import os
os.chdir('/Users/arifeinberg/target2/fs03/d0/arifein/python/pythonHgBenchmark')

from load_Hgmodel_data import open_Hg_spc
from TGMAndObs import SurfaceObsTGM
from matplotlib.backends.backend_pdf import PdfPages

run_old = '0003'
run_new = '0007'
fn_old = '../../GEOS-Chem_runs/run' + run_old + '/OutputDir/GEOSChem.SpeciesConc.alltime_m.nc4'
fn_new = '../../GEOS-Chem_runs/run' + run_new + '/OutputDir/GEOSChem.SpeciesConc.alltime_m.nc4'

ds1, ds2 = open_Hg_spc(fn_old, fn_new)
#%%
plot1, plot2, plot3 = SurfaceObsTGM(ds1, ds2, 2015)

# Save all figures to one PDF file
pp  = PdfPages(('Figures/benchmark_' + run_old + '_' + run_new + '.pdf'))
pp.savefig(plot1)
pp.savefig(plot2)
pp.savefig(plot3)

pp.close()