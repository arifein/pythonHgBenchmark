#!/bin/sh

# Postprocess simulation data for easy analysis in python
module load cdo/1.9.10_oel8

# Simulation name
RUN_NAME=0006
DIR_NAME=/net/fs03/d0/arifein/GEOS-Chem_runs/run${RUN_NAME}/OutputDir

cd ${DIR_NAME}

# # merge into one timeseries
# cdo mergetime GEOSChem.MercuryChem.*_0000z.nc4 GEOSChem.MercuryChem.alltime.nc4
# 
# # take monthly average of timeseries
# cdo monmean GEOSChem.MercuryChem.alltime.nc4 GEOSChem.MercuryChem.alltime_m.nc4
# 
# # sum up individual components of gross oxidation 
# cdo selvar,ProdHg2fromBr,ProdHg2fromBrY,ProdHg2fromClY,ProdHg2fromHgBrPlusBr2,ProdHg2fromHgBrPlusBrBrO,ProdHg2fromHgBrPlusBrClO,ProdHg2fromHgBrPlusBrHO2,ProdHg2fromHgBrPlusBrNO2,ProdHg2fromHgBrPlusBrOH,ProdHg2fromO3,ProdHg2fromOH GEOSChem.MercuryChem.alltime_m.nc4 temp1.nc4
# cdo expr,'Gross_Hg_Ox=ProdHg2fromBr+ProdHg2fromBrY+ProdHg2fromClY+ProdHg2fromHgBrPlusBr2+ProdHg2fromHgBrPlusBrBrO+ProdHg2fromHgBrPlusBrClO+ProdHg2fromHgBrPlusBrHO2+ProdHg2fromHgBrPlusBrNO2+ProdHg2fromHgBrPlusBrOH+ProdHg2fromO3+ProdHg2fromOH' temp1.nc4 temp_gross_ox.nc4 # take sum for gross oxidation
# 
# # move gross oxidation variable into GEOSChem.MercuryChem.alltime_m.nc4 file
# ncks -A -v Gross_Hg_Ox temp_gross_ox.nc4 GEOSChem.MercuryChem.alltime_m.nc4
# 
# # edit metadata of gross oxidation variable
# ncatted -h -O -a  long_name,Gross_Hg_Ox,o,c,"Gross oxidation of Hg0 to Hg2+" GEOSChem.MercuryChem.alltime_m.nc4
# ncatted -h -O -a  units,Gross_Hg_Ox,o,c,"kg s-1" GEOSChem.MercuryChem.alltime_m.nc4 

# take vertical sum of sea salt uptake over layers
cdo selvar,LossHg2bySeaSalt GEOSChem.MercuryChem.alltime_m.nc4 temp2.nc4
cdo vertsum temp2.nc4 temp3.nc4

# rename variable and move into GEOSChem.MercuryChem file
ncrename -v LossHg2bySeaSalt,LossHg2bySeaSalt_v temp3.nc4
ncks -A -v LossHg2bySeaSalt_v temp3.nc4 GEOSChem.MercuryChem.alltime_m.nc4

# remove temporary files
rm temp*.nc4
