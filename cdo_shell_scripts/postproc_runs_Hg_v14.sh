#!/bin/sh

# Postprocess v14 GEOS-Chem Hg run so that files are ready for benchmarking scripts 
module load cdo/1.9.10_oel8

# Simulation name
RUN_NAME=0208
DIR_NAME=/net/fs03/d0/arifein/GEOS-Chem_runs/run${RUN_NAME}/OutputDir

cd ${DIR_NAME}

#---------------------------------------------
# Edits to SpeciesConc monthly mean file
#---------------------------------------------
cdo mergetime GEOSChem.SpeciesConc.*_0000z.nc4 temp_GEOSChem.SpeciesConc.alltime.nc4
cdo monmean temp_GEOSChem.SpeciesConc.alltime.nc4 GEOSChem.SpeciesConc.alltime_m.nc4

mv GEOSChem.SpeciesConc.alltime_m.nc4 temp1.nc4
cdo selvar,AREA,SpeciesConc_Hg0 temp1.nc4 temp_Hg0.nc4 # select Hg0 variable
cdo expr,'SpeciesConc_Hg2=SpeciesConc_HgCl2+SpeciesConc_HgOHOH+SpeciesConc_HgOHBrO+SpeciesConc_HgOHClO+SpeciesConc_HgOHHO2+SpeciesConc_HgOHNO2+SpeciesConc_HgOH+SpeciesConc_HgClOH+SpeciesConc_HgClBr+SpeciesConc_HgClBrO+SpeciesConc_HgClClO+SpeciesConc_HgClHO2+SpeciesConc_HgClNO2+SpeciesConc_HgCl+SpeciesConc_HgBr2+SpeciesConc_HgBrOH+SpeciesConc_HgBrClO+SpeciesConc_HgBrBrO+SpeciesConc_HgBrHO2+SpeciesConc_HgBrNO2+SpeciesConc_HgBr' temp1.nc4 temp_Hg2.nc4 # select Hg2 variable
cdo expr,'SpeciesConc_HgP=SpeciesConc_Hg2STRP+SpeciesConc_Hg2ORGP+SpeciesConc_Hg2ClP' temp1.nc4 temp_HgP.nc4 # select HgP variable

cdo merge temp_Hg0.nc4 temp_Hg2.nc4 temp_HgP.nc4 GEOSChem.SpeciesConc.alltime_m.nc4
#---------------------------------------------
# Edits to Wet deposition monthly mean file
#---------------------------------------------
# merge into one timeseries
cdo mergetime GEOSChem.WetLossLS.*_0000z.nc4 temp_LS_alltime.nc4
cdo mergetime GEOSChem.WetLossConv.*_0000z.nc4 temp_CV_alltime.nc4

# take monthly average of timeseries
cdo monmean temp_LS_alltime.nc4 temp_LS_alltime_m.nc4 
cdo monmean temp_CV_alltime.nc4 temp_CV_alltime_m.nc4
 
# sum up individual files
cdo aexpr,'LS_TOT=WetLossLS_Hg2ORGP+WetLossLS_Hg2ClP+WetLossLS_HgCl2+WetLossLS_HgOHOH+WetLossLS_HgOHBrO+WetLossLS_HgOHClO+WetLossLS_HgOHHO2+WetLossLS_HgOHNO2+WetLossLS_HgClOH+WetLossLS_HgClBr+WetLossLS_HgClBrO+WetLossLS_HgClClO+WetLossLS_HgClHO2+WetLossLS_HgClNO2+WetLossLS_HgBr2+WetLossLS_HgBrOH+WetLossLS_HgBrClO+WetLossLS_HgBrBrO+WetLossLS_HgBrHO2+WetLossLS_HgBrNO2' temp_LS_alltime_m.nc4 temp_LS_total_m.nc4 # total LS wet deposition 
cdo aexpr,'CV_TOT=WetLossConv_Hg2ORGP+WetLossConv_Hg2ClP+WetLossConv_HgCl2+WetLossConv_HgOHOH+WetLossConv_HgOHBrO+WetLossConv_HgOHClO+WetLossConv_HgOHHO2+WetLossConv_HgOHNO2+WetLossConv_HgClOH+WetLossConv_HgClBr+WetLossConv_HgClBrO+WetLossConv_HgClClO+WetLossConv_HgClHO2+WetLossConv_HgClNO2+WetLossConv_HgBr2+WetLossConv_HgBrOH+WetLossConv_HgBrClO+WetLossConv_HgBrBrO+WetLossConv_HgBrHO2+WetLossConv_HgBrNO2' temp_CV_alltime_m.nc4 temp_CV_total_m.nc4 # total Conv wet deposition 
cdo merge temp_LS_total_m.nc4 temp_CV_total_m.nc4 temp_both_total_m.nc4
cdo expr,'WetLossTot_Hg=LS_TOT+CV_TOT' temp_both_total_m.nc4 temp_sum_m1.nc4
cdo expr,'WetLossTot_HgP=WetLossLS_Hg2ORGP+WetLossLS_Hg2ClP+WetLossConv_Hg2ORGP+WetLossConv_Hg2ClP' temp_both_total_m.nc4 temp_sum_m_P.nc4
cdo expr,'WetLossTot_Hg2=WetLossLS_HgCl2+WetLossLS_HgOHOH+WetLossLS_HgOHBrO+WetLossLS_HgOHClO+WetLossLS_HgOHHO2+WetLossLS_HgOHNO2+WetLossLS_HgClOH+WetLossLS_HgClBr+WetLossLS_HgClBrO+WetLossLS_HgClClO+WetLossLS_HgClHO2+WetLossLS_HgClNO2+WetLossLS_HgBr2+WetLossLS_HgBrOH+WetLossLS_HgBrClO+WetLossLS_HgBrBrO+WetLossLS_HgBrHO2+WetLossLS_HgBrNO2+WetLossConv_HgCl2+WetLossConv_HgOHOH+WetLossConv_HgOHBrO+WetLossConv_HgOHClO+WetLossConv_HgOHHO2+WetLossConv_HgOHNO2+WetLossConv_HgClOH+WetLossConv_HgClBr+WetLossConv_HgClBrO+WetLossConv_HgClClO+WetLossConv_HgClHO2+WetLossConv_HgClNO2+WetLossConv_HgBr2+WetLossConv_HgBrOH+WetLossConv_HgBrClO+WetLossConv_HgBrBrO+WetLossConv_HgBrHO2+WetLossConv_HgBrNO2' temp_both_total_m.nc4 temp_sum_m_2.nc4

cdo merge temp_sum_m1.nc4 temp_sum_m_P.nc4 temp_sum_m_2.nc4 temp_sum_m.nc4 

# sum over all levels to get total deposition at surface
cdo vertsum temp_sum_m.nc4 GEOSChem.WetLossTotal.alltime_m.nc4
 
# edit metadata
ncatted -h -O -a  long_name,WetLossTot_Hg,o,c,"Total loss of Hg species by wet deposition" GEOSChem.WetLossTotal.alltime_m.nc4
ncatted -h -O -a  units,WetLossTot_Hg,o,c,"kg s-1" GEOSChem.WetLossTotal.alltime_m.nc4
ncatted -h -O -a  long_name,WetLossTot_HgP,o,c,"Total loss of particle HgP by wet deposition" GEOSChem.WetLossTotal.alltime_m.nc4
ncatted -h -O -a  units,WetLossTot_HgP,o,c,"kg s-1" GEOSChem.WetLossTotal.alltime_m.nc4
ncatted -h -O -a  long_name,WetLossTot_Hg2,o,c,"Total loss of gas Hg2 by wet deposition" GEOSChem.WetLossTotal.alltime_m.nc4
ncatted -h -O -a  units,WetLossTot_Hg2,o,c,"kg s-1" GEOSChem.WetLossTotal.alltime_m.nc4

#---------------------------------------------
# Edits to Dry deposition monthly mean file
#---------------------------------------------
cdo mergetime GEOSChem.DryDep.*_0000z.nc4 temp_GEOSChem.DryDep.alltime.nc4
cdo monmean temp_GEOSChem.DryDep.alltime.nc4 GEOSChem.DryDep.alltime_m.nc4

mv GEOSChem.DryDep.alltime_m.nc4 temp1.nc4
cdo selvar,AREA,DryDep_Hg0 temp1.nc4 temp_Hg0_dd.nc4 # select Hg0 variable
cdo expr,'DryDep_Hg2=DryDep_HgCl2+DryDep_HgOHOH+DryDep_HgOHBrO+DryDep_HgOHClO+DryDep_HgOHHO2+DryDep_HgOHNO2+DryDep_HgClOH+DryDep_HgClBr+DryDep_HgClBrO+DryDep_HgClClO+DryDep_HgClHO2+DryDep_HgClNO2+DryDep_HgBr2+DryDep_HgBrOH+DryDep_HgBrClO+DryDep_HgBrBrO+DryDep_HgBrHO2+DryDep_HgBrNO2' temp1.nc4 temp_Hg2_dd.nc4 # select Hg2 variable
cdo expr,'DryDep_HgP=DryDep_Hg2ORGP+DryDep_Hg2ClP' temp1.nc4 temp_HgP_dd.nc4 # select HgP variable

cdo merge temp_Hg0_dd.nc4 temp_Hg2_dd.nc4 temp_HgP_dd.nc4 GEOSChem.DryDep.alltime_m.nc4
#--------------------------------------------
# Monthly averages of mercury emissions/ocean files 
#---------------------------------------------
cdo mergetime GEOSChem.MercuryEmis.*_0000z.nc4 temp_GEOSChem.MercuryEmis.alltime.nc4
cdo mergetime GEOSChem.MercuryOcean.*_0000z.nc4 temp_GEOSChem.MercuryOcean.alltime.nc4

cdo monmean temp_GEOSChem.MercuryEmis.alltime.nc4 temp_GEOSChem.MercuryEmis.alltime_m.nc4
cdo monmean temp_GEOSChem.MercuryOcean.alltime.nc4 GEOSChem.MercuryOcean.alltime_m.nc4

cdo mergetime HEMCO_diagnostics*.nc HEMCO_diagnostics_alltime.nc
# merge some of the HEMCO emissions in with the MercuryEmis files
cdo -O -L -mul -selvar,EmisHg2ClP_Anthro,EmisHgCl2_Anthro,EmisHg0_Natural,EmisHg0_BioBurn,EmisHg0_Artisanal,EmisHg0_Anthro HEMCO_diagnostics_alltime.nc -gridarea -selgrid,2 HEMCO_diagnostics_alltime.nc temp_hemco_1.nc4
cdo expr,'EmisHg2HgPanthro=EmisHg2ClP_Anthro+EmisHgCl2_Anthro' temp_hemco_1.nc4 temp_hemco_2.nc4 # select variable
cdo expr,'EmisHg0anthro=EmisHg0_Anthro+EmisHg0_Artisanal' temp_hemco_1.nc4 temp_hemco_3.nc4 # select variable
cdo selvar,EmisHg0_Natural,EmisHg0_BioBurn temp_hemco_1.nc4 temp_hemco_4.nc4 
ncrename -v EmisHg0_Natural,EmisHg0geogenic temp_hemco_4.nc4
ncrename -v EmisHg0_BioBurn,EmisHg0biomass temp_hemco_4.nc4
cdo merge temp_hemco_2.nc4 temp_hemco_3.nc4 temp_hemco_4.nc4 temp_GEOSChem.MercuryEmis.alltime_m.nc4 GEOSChem.MercuryEmis.alltime_m.nc4
ncatted -h -O -a  units,EmisHg2HgPanthro,o,c,"kg s-1" GEOSChem.MercuryEmis.alltime_m.nc4
ncatted -h -O -a  units,EmisHg0anthro,o,c,"kg s-1" GEOSChem.MercuryEmis.alltime_m.nc4
ncatted -h -O -a  units,EmisHg0geogenic,o,c,"kg s-1" GEOSChem.MercuryEmis.alltime_m.nc4
ncatted -h -O -a  units,EmisHg0biomass,o,c,"kg s-1" GEOSChem.MercuryEmis.alltime_m.nc4


#---------------------------------------------
# Edits to ProdLoss files 
#---------------------------------------------
cdo mergetime GEOSChem.ProdLoss.*_0000z.nc4 temp_GEOSChem.ProdLoss.alltime.nc4
cdo monmean temp_GEOSChem.ProdLoss.alltime.nc4 temp1_pl.nc4
# 
#cdo aexpr,'Gross_Hg_Ox=Prod_Hg2Br+Prod_Hg2OH+Prod_Hg2Cl' temp1_pl.nc4 temp2_pl.nc4 # Gross GOM oxidation # 0205
cdo aexpr,'Gross_Hg_Ox=Prod_Hg2' temp1_pl.nc4 temp2_pl.nc4 # Gross GOM oxidation 
cdo aexpr,'NetOx=Gross_Hg_Ox-Prod_Hg0' temp2_pl.nc4 GEOSChem.ProdLoss.alltime_m.nc4 # Net GOM oxidation

#---------------------------------------------
# Editing mercury chemistry files
#---------------------------------------------
cdo mergetime GEOSChem.MercuryChem.*_0000z.nc4 temp_GEOSChem.MercuryChem.alltime.nc4
cdo monmean temp_GEOSChem.MercuryChem.alltime.nc4 temp_GEOSChem.MercuryChem.alltime_m.nc4
cdo selvar,Hg2GasToSSA temp_GEOSChem.MercuryChem.alltime_m.nc4 temp_seasalt.nc4

# convert units from molec/cm3/s to kg/s
#cdo -O -L selvar,Met_AIRVOL -selyear,2014 /net/fs03/d0/arifein/GEOS-Chem_runs/run0009/OutputDir/GEOSChem.StateMet.alltime_m.nc4 airvol.nc4
cdo -O -L selvar,Met_AIRVOL /net/fs03/d0/arifein/GEOS-Chem_runs/run0009/OutputDir/GEOSChem.StateMet.alltime_m.nc4 airvol.nc4
ncks -A -v Met_AIRVOL airvol.nc4 temp_seasalt.nc4
cdo expr,'LossHg2bySeaSalt=Hg2GasToSSA*Met_AIRVOL' temp_seasalt.nc4 temp_seasalt2.nc4
# unit conv  = 200.59 (g mol^-1) /1000 (g kg^-1) /6.02*10^23 (molec mol^-1) * 10^6 (m^3/cm^3)= 3.332e-19
cdo mulc,3.332059800664452e-19 temp_seasalt2.nc4 temp_seasalt3.nc4

cdo vertsum temp_seasalt3.nc4 temp_seasalt4.nc4

ncrename -v LossHg2bySeaSalt,LossHg2bySeaSalt_v temp_seasalt4.nc4
ncks -A -v LossHg2bySeaSalt_v temp_seasalt4.nc4 GEOSChem.MercuryChem.alltime_m.nc4
ncks -A -v LossHg2bySeaSalt temp_seasalt3.nc4 GEOSChem.MercuryChem.alltime_m.nc4

#add chemistry fluxes
ncks -A -v Gross_Hg_Ox GEOSChem.ProdLoss.alltime_m.nc4 temp_chem.nc4
ncks -A -v NetOx GEOSChem.ProdLoss.alltime_m.nc4 temp_chem.nc4
ncks -A -v Met_AIRVOL airvol.nc4 temp_chem.nc4

cdo expr,'Gross_Hg_Ox=Gross_Hg_Ox*Met_AIRVOL' temp_chem.nc4 temp_chem1.nc4
cdo mulc,3.332059800664452e-19 temp_chem1.nc4 temp_chem1u.nc4

cdo expr,'ProdHg2fromHg0=NetOx*Met_AIRVOL' temp_chem.nc4 temp_chem2.nc4
cdo mulc,3.332059800664452e-19 temp_chem2.nc4 temp_chem2u.nc4

ncks -A -v Gross_Hg_Ox temp_chem1u.nc4 GEOSChem.MercuryChem.alltime_m.nc4
ncks -A -v ProdHg2fromHg0 temp_chem2u.nc4 GEOSChem.MercuryChem.alltime_m.nc4
#---------------------------------------------
# Create mercury budget files
#---------------------------------------------
#---calculate global total emissions from different sources (kg/s)---
cdo -O -L -fldsum -selvar,EmisHg0anthro,EmisHg2HgPanthro,EmisHg0geogenic,EmisHg0biomass,EmisHg0soil,EmisHg0land,EmisHg0snow,EmisHg0ocean GEOSChem.MercuryEmis.alltime_m.nc4 temp_emiss.nc4

#-------------------------------------------------------------------

#---calculate global total fluxes in air-sea exchange (kg/s)-------
cdo -O -L -fldsum -selvar,FluxHg0fromOceanToAir,FluxHg0fromAirToOcean GEOSChem.MercuryOcean.alltime_m.nc4 temp_ocean.nc4
#-------------------------------------------------------------------

#---calculate global total wet deposition (kg/s)--------------------
cdo -O -L -fldsum -vertsum -selvar,WetLossTot_Hg2,WetLossTot_HgP GEOSChem.WetLossTotal.alltime_m.nc4 temp_wetdep.nc4
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
cdo -O -L -fldsum -vertsum -selvar,Gross_Hg_Ox,ProdHg2fromHg0 GEOSChem.MercuryChem.alltime_m.nc4 temp_chem1.nc4
cdo -O -L -fldsum -selvar,LossHg2bySeaSalt_v GEOSChem.MercuryChem.alltime_m.nc4 temp_chem2.nc4
ncrename -v LossHg2bySeaSalt_v,LossHg2bySeaSalt temp_chem2.nc4
cdo merge temp_chem1.nc4 temp_chem2.nc4 temp_chem3.nc4
#------------------------------------------------------------------


#---merge files, clean up metadata, remove temporary files---------
# merge intermediate files
cdo merge temp_emiss.nc4 temp_wetdep.nc4 temp_drydep.nc4 temp_chem3.nc4 temp_ocean.nc4 GEOSChem.MercuryBudget_global_m.nc4

# change units
ncatted -h -O -a  units,WetLossTot_Hg2,o,c,"kg s-1" GEOSChem.MercuryBudget_global_m.nc4
ncatted -h -O -a  units,WetLossTot_HgP,o,c,"kg s-1" GEOSChem.MercuryBudget_global_m.nc4
ncatted -h -O -a  units,DryDep_HgP,o,c,"kg s-1" GEOSChem.MercuryBudget_global_m.nc4
ncatted -h -O -a  units,DryDep_Hg2,o,c,"kg s-1" GEOSChem.MercuryBudget_global_m.nc4
ncatted -h -O -a  units,DryDep_Hg0,o,c,"kg s-1" GEOSChem.MercuryBudget_global_m.nc4

#---------------------------------------------
# remove temporary files
#---------------------------------------------
rm temp*.nc4
