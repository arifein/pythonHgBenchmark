#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 14:20:15 2021
Main file for running Hg benchmark in python and producing plots
@author: arifeinberg
"""

import os
os.chdir('/Users/arifeinberg/target2/fs03/d0/arifein/python/pythonHgBenchmark')

from helper_functions import open_Hg
from TGMAndObs import SurfaceObsTGM
from Hg2Plot import SurfaceHg2
from Latitudinal_Graphs import Seasonal_Lat_Regions, plot_gradient_TGM
from PlotSeasonSites import PlotSeasonSites
from wet_deposition import wet_dep_plots
from matplotlib.backends.backend_pdf import PdfPages

#%% Opening Hg species datasets
run_old = '0003'
run_new = '0007'
fn_old = '../../GEOS-Chem_runs/run' + run_old + '/OutputDir/GEOSChem.SpeciesConc.alltime_m.nc4'
fn_new = '../../GEOS-Chem_runs/run' + run_new + '/OutputDir/GEOSChem.SpeciesConc.alltime_m.nc4'

year_to_analyze = 2015 # year to analyze from simulations 
ds1, ds2 = open_Hg(fn_old, fn_new)
#%% Plots of surface TGM
plot1, plot2, plot3 = SurfaceObsTGM(ds1, ds2, year_to_analyze)

#%% Plot of surface Hg(II) + Hg(P)
plot4 = SurfaceHg2(ds1, ds2, year_to_analyze)

#%% Plot of seasonal cycle for different lat regions
plot5 = Seasonal_Lat_Regions(ds1, ds2, year_to_analyze)

#%% Plot seasonal cycles of TGM for all sites
plot6 = PlotSeasonSites(ds1, ds2, year_to_analyze)

#%% Plot latitudinal gradient of TGM vs. observations
plot7 = plot_gradient_TGM(ds1, ds2, year_to_analyze)
#%% Opening Hg wet deposition datasets
run_old_2 = '0005' # same as 0003, just with correct outputs

# total deposition, summed over all levels (see cdo_shell_scripts/ folder for postprocessing)
fn_old_wdep = '../../GEOS-Chem_runs/run' + run_old_2 + '/OutputDir/GEOSChem.WetLossTotal.alltime_m.nc4'
fn_new_wdep = '../../GEOS-Chem_runs/run' + run_new + '/OutputDir/GEOSChem.WetLossTotal.alltime_m.nc4'

ds1_wdep, ds2_wdep = open_Hg(fn_old_wdep, fn_new_wdep) # load deposition data
#%% Running wet deposition comparison plots
plot8, plot9, plot10, plot11, plot12, plot13 = wet_dep_plots(ds1_wdep, ds2_wdep, year_to_analyze)

#%% Save all figures to one PDF file
pp  = PdfPages(('Figures/benchmark_' + run_old + '_' + run_new + '.pdf'))
pp.savefig(plot1, bbox_inches = 'tight')
pp.savefig(plot2, bbox_inches = 'tight')
pp.savefig(plot3, bbox_inches = 'tight')
pp.savefig(plot4, bbox_inches = 'tight')
pp.savefig(plot5, bbox_inches = 'tight')
pp.savefig(plot6, bbox_inches = 'tight')
pp.savefig(plot7, bbox_inches = 'tight')
pp.savefig(plot8, bbox_inches = 'tight')
pp.savefig(plot9, bbox_inches = 'tight')
pp.savefig(plot10, bbox_inches = 'tight')
pp.savefig(plot11, bbox_inches = 'tight')
pp.savefig(plot12, bbox_inches = 'tight')
pp.savefig(plot13, bbox_inches = 'tight')

pp.close()