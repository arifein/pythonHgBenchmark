#!/bin/sh

# Postprocess simulation data for easy analysis in python
module load cdo/1.9.10_oel8

# Simulation name
RUN_NAME=0007
DIR_NAME=/net/fs03/d0/arifein/GEOS-Chem_runs/run${RUN_NAME}/OutputDir

cd ${DIR_NAME}

# merge into one timeseries
cdo mergetime GEOSChem.WetLossLS.*_0000z.nc4 GEOSChem.WetLossLS.alltime.nc4
cdo mergetime GEOSChem.WetLossConv.*_0000z.nc4 GEOSChem.WetLossConv.alltime.nc4

# take monthly average of timeseries
cdo monmean GEOSChem.WetLossLS.alltime.nc4 GEOSChem.WetLossLS.alltime_m.nc4
cdo monmean GEOSChem.WetLossConv.alltime.nc4 GEOSChem.WetLossConv.alltime_m.nc4

# sum up individual components of wet deposition to get total wet deposition
cdo selvar,WetLossConv_HgP,WetLossConv_Hg2 GEOSChem.WetLossConv.alltime_m.nc4 temp1.nc4 # select Conv variables
cdo selvar,WetLossLS_HgP,WetLossLS_Hg2 GEOSChem.WetLossLS.alltime_m.nc4 temp2.nc4 # select LS variables
cdo merge temp1.nc4 temp2.nc4 temp_allwetdep.nc4 # merge datasets
cdo expr,'WetLossTot_Hg=WetLossConv_HgP+WetLossConv_Hg2+WetLossLS_HgP+WetLossLS_Hg2' temp_allwetdep.nc4 temp_allwetdep_sum.nc4 # take sum

# sum over all levels to get total deposition at surface
cdo vertsum temp_allwetdep_sum.nc4 GEOSChem.WetLossTotal.alltime_m.nc4

# edit metadata
ncatted -h -O -a  long_name,WetLossTot_Hg,o,c,"Total loss of Hg species by wet deposition" GEOSChem.WetLossTotal.alltime_m.nc4
ncatted -h -O -a  units,WetLossTot_Hg,o,c,"kg s-1" GEOSChem.WetLossTotal.alltime_m.nc4

# remove temporary files
rm temp*.nc4
