#!/bin/sh

# Postprocess simulation data for easy analysis in python benchmark
module load cdo/1.9.10_oel8

# Simulation name and directory
RUN_NAME=2000
DIR_NAME=/net/fs03/d0/arifein/GEOS-Chem_runs/run${RUN_NAME}/OutputDir

cd ${DIR_NAME}

# --- Monthly average postprocessing ---
# merge into one timeseries
cdo mergetime GEOSChem.SpeciesConc.*_0000z.nc4 GEOSChem.SpeciesConc.alltime.nc4
cdo mergetime GEOSChem.MercuryEmis.*_0000z.nc4 GEOSChem.MercuryEmis.alltime.nc4
cdo mergetime GEOSChem.DryDep.*_0000z.nc4 GEOSChem.DryDep.alltime.nc4
cdo mergetime GEOSChem.MercuryOcean.*_0000z.nc4 GEOSChem.MercuryOcean.alltime.nc4
cdo mergetime GEOSChem.StateMet.*_0000z.nc4 GEOSChem.StateMet.alltime.nc4
cdo mergetime GEOSChem.WetLossLS.*_0000z.nc4 GEOSChem.WetLossLS.alltime.nc4
cdo mergetime GEOSChem.WetLossConv.*_0000z.nc4 GEOSChem.WetLossConv.alltime.nc4
cdo mergetime GEOSChem.MercuryChem.*_0000z.nc4 GEOSChem.MercuryChem.alltime.nc4

# take monthly average
cdo monmean GEOSChem.SpeciesConc.alltime.nc4 GEOSChem.SpeciesConc.alltime_m.nc4
cdo monmean GEOSChem.MercuryEmis.alltime.nc4 GEOSChem.MercuryEmis.alltime_m.nc4
cdo monmean GEOSChem.DryDep.alltime.nc4 GEOSChem.DryDep.alltime_m.nc4
cdo monmean GEOSChem.MercuryOcean.alltime.nc4 GEOSChem.MercuryOcean.alltime_m.nc4
cdo monmean GEOSChem.StateMet.alltime.nc4 GEOSChem.StateMet.alltime_m.nc4
cdo monmean GEOSChem.WetLossLS.alltime.nc4 GEOSChem.WetLossLS.alltime_m.nc4
cdo monmean GEOSChem.WetLossConv.alltime.nc4 GEOSChem.WetLossConv.alltime_m.nc4
cdo monmean GEOSChem.MercuryChem.alltime.nc4 GEOSChem.MercuryChem.alltime_m.nc4

# --- Wet Deposition postprocessing ---
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

# --- Chemistry postprocessing ---
# sum up individual components of gross oxidation 
cdo selvar,ProdHg2fromBr,ProdHg2fromBrY,ProdHg2fromClY,ProdHg2fromO3,ProdHg2fromOH GEOSChem.MercuryChem.alltime_m.nc4 temp1.nc4
cdo expr,'Gross_Hg_Ox=ProdHg2fromBr+ProdHg2fromBrY+ProdHg2fromClY+ProdHg2fromO3+ProdHg2fromOH' temp1.nc4 temp_gross_ox.nc4 # take sum for gross oxidation

# move gross oxidation variable into GEOSChem.MercuryChem.alltime_m.nc4 file
ncks -A -v Gross_Hg_Ox temp_gross_ox.nc4 GEOSChem.MercuryChem.alltime_m.nc4

# edit metadata of gross oxidation variable
ncatted -h -O -a  long_name,Gross_Hg_Ox,o,c,"Gross oxidation of Hg0 to Hg2+" GEOSChem.MercuryChem.alltime_m.nc4
ncatted -h -O -a  units,Gross_Hg_Ox,o,c,"kg s-1" GEOSChem.MercuryChem.alltime_m.nc4

# take vertical sum of sea salt uptake over layers
cdo selvar,LossHg2bySeaSalt GEOSChem.MercuryChem.alltime_m.nc4 temp2.nc4
cdo vertsum temp2.nc4 temp3.nc4

# rename variable and move into GEOSChem.MercuryChem file
ncrename -v LossHg2bySeaSalt,LossHg2bySeaSalt_v temp3.nc4
ncks -A -v LossHg2bySeaSalt_v temp3.nc4 GEOSChem.MercuryChem.alltime_m.nc4

# remove temporary files
rm temp*.nc4

# --- Mercury burden postprocessing ---

# molar mass Hg = 200.59, molar mass dry air = 28.965, conversion factor to mass mixing ratio = 200.59/28.965 = 6.925
# multiply by grid box air mass, and convert units
cdo -L mulc,6.925 -mul -selvar,Met_AD ../../run0013/OutputDir/GEOSChem.StateMet.alltime_m.nc4 -selvar,SpeciesConc_Hg0,SpeciesConc_Hg2,SpeciesConc_HgP GEOSChem.SpeciesConc.alltime_m.nc4 temp_mul.nc

# calculate total burden by summing vertically and horizontally
cdo -O -L -fldsum -vertsum temp_mul.nc GEOSChem.MercuryBurden_global_m.nc
# create mask for tropospheric values
#cdo -f nc -expr,'mask=lev<TR-PAUSE__TP-LEVEL' temp_dims2.nc imask.nc

# change units
ncatted -h -O -a  units,SpeciesConc_Hg0,o,c,"kg" GEOSChem.MercuryBurden_global_m.nc
ncatted -h -O -a  units,SpeciesConc_HgP,o,c,"kg" GEOSChem.MercuryBurden_global_m.nc
ncatted -h -O -a  units,SpeciesConc_Hg2,o,c,"kg" GEOSChem.MercuryBurden_global_m.nc

# remove temporary files
rm temp*.nc

# --- Mercury budget postprocessing ---
# calculate global total emissions from different sources (kg/s)
cdo -O -L -fldsum -selvar,EmisHg0anthro,EmisHg2HgPanthro,EmisHg0geogenic,EmisHg0soil,EmisHg0biomass,EmisHg0land,EmisHg0snow GEOSChem.MercuryEmis.alltime_m.nc4 temp_emiss.nc4

# calculate global total fluxes in air-sea exchange (kg/s)
cdo -O -L -fldsum -selvar,FluxHg0fromOceanToAir,FluxHg0fromAirToOcean GEOSChem.MercuryOcean.alltime_m.nc4 temp_ocean.nc4

# calculate global total wet deposition (kg/s)
# select variables for wet deposition
cdo selvar,WetLossConv_Hg2,WetLossConv_HgP GEOSChem.WetLossConv.alltime_m.nc4 temp1.nc4
cdo selvar,WetLossLS_Hg2,WetLossLS_HgP GEOSChem.WetLossLS.alltime_m.nc4 temp2.nc4
cdo merge temp1.nc4 temp2.nc4 temp3.nc4

# add variables for each Hg species total wetdep
cdo expr,'WetLossTot_Hg2=WetLossLS_Hg2+WetLossConv_Hg2' temp3.nc4 temp4.nc4
cdo expr,'WetLossTot_HgP=WetLossLS_HgP+WetLossConv_HgP' temp3.nc4 temp5.nc4
cdo merge temp4.nc4 temp5.nc4 temp6.nc4

cdo -O -L -fldsum -vertsum -selvar,WetLossTot_Hg2,WetLossTot_HgP temp6.nc4 temp_wetdep.nc4

# calculate global total dry deposition (kg/s)
# select variables for dry deposition
cdo selvar,DryDep_HgP,DryDep_Hg2,DryDep_Hg0 GEOSChem.DryDep.alltime_m.nc4 temp7.nc4

# convert from molec/cm2/s to kg/m2/s
cdo mulc,3.332059800664452e-21 temp7.nc4 temp8.nc4

# weight with grid area and sum to get kg/s globally
cdo -O -L fldsum -mul temp8.nc4 -gridarea temp8.nc4 temp_drydep.nc4


# calculate global total chemistry fluxes (kg/s)
cdo -O -L -fldsum -vertsum -selvar,LossHg2bySeaSalt,Gross_Hg_Ox,ProdHg2fromHg0 GEOSChem.MercuryChem.alltime_m.nc4 temp_chem.nc4

# merge files, clean up metadata, remove temporary files
# merge intermediate files
cdo merge temp_emiss.nc4 temp_wetdep.nc4 temp_drydep.nc4 temp_chem.nc4 temp_ocean.nc4 GEOSChem.MercuryBudget_global_m.nc4

# change units
ncatted -h -O -a  units,WetLossTot_Hg2,o,c,"kg s-1" GEOSChem.MercuryBudget_global_m.nc4
ncatted -h -O -a  units,WetLossTot_HgP,o,c,"kg s-1" GEOSChem.MercuryBudget_global_m.nc4
ncatted -h -O -a  units,DryDep_HgP,o,c,"kg s-1" GEOSChem.MercuryBudget_global_m.nc4
ncatted -h -O -a  units,DryDep_Hg2,o,c,"kg s-1" GEOSChem.MercuryBudget_global_m.nc4
ncatted -h -O -a  units,DryDep_Hg0,o,c,"kg s-1" GEOSChem.MercuryBudget_global_m.nc4

# remove temporary files
rm temp*.nc4


