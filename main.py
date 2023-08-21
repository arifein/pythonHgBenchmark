#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 14:20:15 2021
Main file for running Hg benchmark in python and producing plots
@author: arifeinberg
"""

# import os
# os.chdir('/Users/arifeinberg/target2/fs03/d0/arifein/python/pythonHgBenchmark')
import matplotlib
matplotlib.use('Agg')

from helper_functions import open_Hg
from TGMAndObs import SurfaceObsTGM
from Hg2Plot import SurfaceHg2
from Latitudinal_Graphs import Seasonal_Lat_Regions, plot_gradient_TGM
from PlotSeasonSites import PlotSeasonSites
from wet_deposition import wet_dep_plots
from dry_deposition import dry_dep_plots
from Hg_chem_fluxes import chem_plots
from Hg_emis_fluxes import emis_plots
from Hg_ocean_fluxes import ocean_plots
from Hg_budget import budget_calc

from matplotlib.backends.backend_pdf import PdfPages

# create list of plots to add to
plotlist = []
#%% Opening Hg species datasets
run_old = '2000'
run_new = '2010'
fn_old = '../../GEOS-Chem_runs/run' + run_old + '/OutputDir/GEOSChem.SpeciesConc.alltime_m.nc4'
fn_new = '../../GEOS-Chem_runs/run' + run_new + '/OutputDir/GEOSChem.SpeciesConc.alltime_m.nc4'

year_1 = 2005 # year to analyze from old simulation 
year_2 = 2005 # year to analyze from new simulation 

ds1, ds2 = open_Hg(fn_old, fn_new)

#%% Plots of surface TGM
# add output plots to list of plots
plotlist.extend([SurfaceObsTGM(ds1, ds2, year_1, year_2)])

#%% Plot of surface Hg(II) + Hg(P)
plotlist.extend([SurfaceHg2(ds1, ds2, year_1, year_2)])

#%% Plot of seasonal cycle for different lat regions
plotlist.extend([Seasonal_Lat_Regions(ds1, ds2, year_1, year_2)])

#%% Plot seasonal cycles of TGM for all sites
plotlist.extend([PlotSeasonSites(ds1, ds2, year_1, year_2)])

#%% Plot latitudinal gradient of TGM vs. observations
plotlist.extend([plot_gradient_TGM(ds1, ds2, year_1, year_2)])
#%% Opening Hg wet deposition datasets
# total deposition, summed over all levels (see cdo_shell_scripts/ folder for postprocessing)
fn_old_wdep = '../../GEOS-Chem_runs/run' + run_old + '/OutputDir/GEOSChem.WetLossTotal.alltime_m.nc4'
fn_new_wdep = '../../GEOS-Chem_runs/run' + run_new + '/OutputDir/GEOSChem.WetLossTotal.alltime_m.nc4'

ds1_wdep, ds2_wdep = open_Hg(fn_old_wdep, fn_new_wdep) # load wet deposition data
#%% Running wet deposition comparison plots
plotlist.extend([wet_dep_plots(ds1_wdep, ds2_wdep, year_1, year_2)])
#%% Opening Hg dry deposition datasets

fn_old_ddep = '../../GEOS-Chem_runs/run' + run_old + '/OutputDir/GEOSChem.DryDep.alltime_m.nc4'
fn_new_ddep = '../../GEOS-Chem_runs/run' + run_new + '/OutputDir/GEOSChem.DryDep.alltime_m.nc4'

ds1_ddep, ds2_ddep = open_Hg(fn_old_ddep, fn_new_ddep) # load dry deposition data
#%% Running dry deposition comparison plots
plotlist.extend([dry_dep_plots(ds1_ddep, ds2_ddep, year_1, year_2)])

#%% Opening Hg chemistry datasets
fn_old_chem = '../../GEOS-Chem_runs/run' + run_old + '/OutputDir/GEOSChem.MercuryChem.alltime_m.nc4'
fn_new_chem = '../../GEOS-Chem_runs/run' + run_new + '/OutputDir/GEOSChem.MercuryChem.alltime_m.nc4'

ds1_chem, ds2_chem = open_Hg(fn_old_chem, fn_new_chem) # load chemistry data
#%% Running Hg chemistry plots
plotlist.extend([chem_plots(ds1_chem, ds2_chem, year_1, year_2)])

#%% Opening Hg emission datasets
fn_old_emis = '../../GEOS-Chem_runs/run' + run_old + '/OutputDir/GEOSChem.MercuryEmis.alltime_m.nc4'
fn_new_emis = '../../GEOS-Chem_runs/run' + run_new + '/OutputDir/GEOSChem.MercuryEmis.alltime_m.nc4'

ds1_emis, ds2_emis = open_Hg(fn_old_emis, fn_new_emis) # load emissions data
#%% Running Hg emissions plots
plotlist.extend([emis_plots(ds1_emis, ds2_emis, year_1, year_2)])

# #%% Opening Hg ocean datasets
# fn_old_ocean = '../../GEOS-Chem_runs/run' + run_old + '/OutputDir/GEOSChem.MercuryOcean.alltime_m.nc4'
# fn_new_ocean = '../../GEOS-Chem_runs/run' + run_new + '/OutputDir/GEOSChem.MercuryOcean.alltime_m.nc4'

# ds1_ocean, ds2_ocean = open_Hg(fn_old_ocean, fn_new_ocean) # load ocean data
# #%% Running Hg oceans plots
# plotlist.extend([ocean_plots(ds1_ocean, ds2_ocean, year_1, year_2)])

#%% Opening Hg budget datasets
fn_old_budget = '../../GEOS-Chem_runs/run' + run_old + '/OutputDir/GEOSChem.MercuryBudget_global_m.nc4'
fn_new_budget = '../../GEOS-Chem_runs/run' + run_new + '/OutputDir/GEOSChem.MercuryBudget_global_m.nc4'

ds1_budget, ds2_budget = open_Hg(fn_old_budget, fn_new_budget) # load ocean data
#%% Running Hg budget calc
plotlist.extend([budget_calc(ds1_budget, ds2_budget, year_1, year_2)])

#%% Flatten the list of plots into one list
flat_plotlist = []
for sublist in plotlist:
    if isinstance(sublist, list): 
        for item in sublist:
            flat_plotlist.append(item)
    else: # can't iterate over Figures
        flat_plotlist.append(sublist)
#%% Save all figures to one PDF file
if year_1 != year_2: # use different years for simulation
    fn_fig = 'Figures/benchmark_' + run_old + '_' + run_new + '_' + str(year_1) + '-' + str(year_2) + '.pdf'
else:
    fn_fig = 'Figures/benchmark_' + run_old + '_' + run_new + '.pdf'
pp  = PdfPages((fn_fig))

# Add all plots to PDF as a loop
for iplot in flat_plotlist:
    pp.savefig(iplot, bbox_inches = 'tight')
    
pp.close()
