#!/bin/sh

# Postprocess simulation data for easy analysis in python
module load cdo/1.9.10_oel8

# Simulation name
RUN_NAME=0007
DIR_NAME=/net/fs03/d0/arifein/GEOS-Chem_runs/run${RUN_NAME}/OutputDir

cd ${DIR_NAME}

#---calculate global total emissions from different sources (kg/s)---
cdo -O -L -fldsum -selvar,EmisHg0anthro,EmisHg2HgPanthro,EmisHg0geogenic,EmisHg0soil,EmisHg0biomass,EmisHg0land,EmisHg0snow GEOSChem.MercuryEmis.alltime_m.nc4 temp_emiss.nc4
#-------------------------------------------------------------------

#---calculate global total fluxes in air-sea exchange (kg/s)-------
cdo -O -L -fldsum -selvar,FluxHg0fromOceanToAir,FluxHg0fromAirToOcean GEOSChem.MercuryOcean.alltime_m.nc4 temp_ocean.nc4
#-------------------------------------------------------------------

#---calculate global total wet deposition (kg/s)--------------------
# select variables for wet deposition
cdo selvar,WetLossConv_Hg2,WetLossConv_HgP GEOSChem.WetLossConv.alltime_m.nc4 temp1.nc4
cdo selvar,WetLossLS_Hg2,WetLossLS_HgP GEOSChem.WetLossLS.alltime_m.nc4 temp2.nc4
cdo merge temp1.nc4 temp2.nc4 temp3.nc4

# add variables for each Hg species total wetdep
cdo expr,'WetLossTot_Hg2=WetLossLS_Hg2+WetLossConv_Hg2' temp3.nc4 temp4.nc4
cdo expr,'WetLossTot_HgP=WetLossLS_HgP+WetLossConv_HgP' temp3.nc4 temp5.nc4
cdo merge temp4.nc4 temp5.nc4 temp6.nc4

cdo -O -L -fldsum -vertsum -selvar,WetLossTot_Hg2,WetLossTot_HgP temp6.nc4 temp_wetdep.nc4
#-------------------------------------------------------------------


#---calculate global total dry deposition (kg/s)--------------------
# select variables for dry deposition
cdo selvar,DryDep_HgP,DryDep_Hg2,DryDep_Hg0 GEOSChem.DryDep.alltime_m.nc4 temp7.nc4

# convert from molec/cm2/s to kg/m2/s
cdo mulc,3.332059800664452e-21 temp7.nc4 temp8.nc4

# weight with grid area and sum to get kg/s globally
cdo -O -L fldsum -mul temp8.nc4 -gridarea temp8.nc4 temp_drydep.nc4 
#------------------------------------------------------------------


#---calculate global total chemistry fluxes (kg/s)------------------
cdo -O -L -fldsum -vertsum -selvar,LossHg2bySeaSalt,Gross_Hg_Ox,ProdHg2fromHg0 GEOSChem.MercuryChem.alltime_m.nc4 temp_chem.nc4
#------------------------------------------------------------------


#---merge files, clean up metadata, remove temporary files---------
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
#------------------------------------------------------------------
