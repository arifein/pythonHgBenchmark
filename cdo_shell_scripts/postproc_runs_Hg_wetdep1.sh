#!/bin/sh

# Postprocess simulation data for easy analysis in python
module load cdo/1.9.10_oel8

# Simulation name
RUN_NAME=0006
DIR_NAME=/net/fs03/d0/arifein/GEOS-Chem_runs/run${RUN_NAME}/OutputDir

cd ${DIR_NAME}

#merge into one timeseries
cdo mergetime GEOSChem.WetLossLS.*_0000z.nc4 GEOSChem.WetLossLS.alltime.nc4
cdo mergetime GEOSChem.WetLossConv.*_0000z.nc4 GEOSChem.WetLossConv.alltime.nc4

cdo monmean GEOSChem.WetLossLS.alltime.nc4 GEOSChem.WetLossLS.alltime_m.nc4
cdo monmean GEOSChem.WetLossConv.alltime.nc4 GEOSChem.WetLossConv.alltime_m.nc4
