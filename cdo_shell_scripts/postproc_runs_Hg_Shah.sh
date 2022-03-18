#!/bin/sh

# Postprocess Shah et al. 2021 Hg chemistry run so that similar outputs to previous version GEOS-Chem
module load cdo/1.9.10_oel8

# Simulation name
RUN_NAME=0105
DIR_NAME=/net/fs03/d0/arifein/GEOS-Chem_runs/run${RUN_NAME}/OutputDir

cd ${DIR_NAME}

#---------------------------------------------
# Edits to SpeciesConc monthly mean file
#---------------------------------------------
cdo mergetime GEOSChem.SpeciesConc.*_0000z.nc4 temp_GEOSChem.SpeciesConc.alltime.nc4
cdo monmean temp_GEOSChem.SpeciesConc.alltime.nc4 GEOSChem.SpeciesConc.alltime_m.nc4

mv GEOSChem.SpeciesConc.alltime_m.nc4 temp1.nc4
cdo selvar,AREA,SpeciesConc_Hg0 temp1.nc4 temp_Hg0.nc4 # select Hg0 variable
cdo expr,'SpeciesConc_Hg2=SpeciesConc_HGCL2+SpeciesConc_HGOHOH+SpeciesConc_HGOHBRO+SpeciesConc_HGOHCLO+SpeciesConc_HGOHHO2+SpeciesConc_HGOHNO2+SpeciesConc_HGOH+SpeciesConc_HGCLOH+SpeciesConc_HGCLBR+SpeciesConc_HGCLBRO+SpeciesConc_HGCLCLO+SpeciesConc_HGCLHO2+SpeciesConc_HGCLNO2+SpeciesConc_HGCL+SpeciesConc_HGBR2+SpeciesConc_HGBROH+SpeciesConc_HGBRCLO+SpeciesConc_HGBRBRO+SpeciesConc_HGBRHO2+SpeciesConc_HGBRNO2+SpeciesConc_HGBR' temp1.nc4 temp_Hg2.nc4 # select Hg2 variable
cdo expr,'SpeciesConc_HgP=SpeciesConc_HG2STRP+SpeciesConc_HG2ORGP+SpeciesConc_HG2CLP' temp1.nc4 temp_HgP.nc4 # select HgP variable

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
cdo aexpr,'LS_TOT=WetLossLS_HG2ORGP+WetLossLS_HG2CLP+WetLossLS_HGCL2+WetLossLS_HGOHOH+WetLossLS_HGOHBRO+WetLossLS_HGOHCLO+WetLossLS_HGOHHO2+WetLossLS_HGOHNO2+WetLossLS_HGCLOH+WetLossLS_HGCLBR+WetLossLS_HGCLBRO+WetLossLS_HGCLCLO+WetLossLS_HGCLHO2+WetLossLS_HGCLNO2+WetLossLS_HGBR2+WetLossLS_HGBROH+WetLossLS_HGBRCLO+WetLossLS_HGBRBRO+WetLossLS_HGBRHO2+WetLossLS_HGBRNO2' temp_LS_alltime_m.nc4 temp_LS_total_m.nc4 # total LS wet deposition 
cdo aexpr,'CV_TOT=WetLossConv_HG2ORGP+WetLossConv_HG2CLP+WetLossConv_HGCL2+WetLossConv_HGOHOH+WetLossConv_HGOHBRO+WetLossConv_HGOHCLO+WetLossConv_HGOHHO2+WetLossConv_HGOHNO2+WetLossConv_HGCLOH+WetLossConv_HGCLBR+WetLossConv_HGCLBRO+WetLossConv_HGCLCLO+WetLossConv_HGCLHO2+WetLossConv_HGCLNO2+WetLossConv_HGBR2+WetLossConv_HGBROH+WetLossConv_HGBRCLO+WetLossConv_HGBRBRO+WetLossConv_HGBRHO2+WetLossConv_HGBRNO2' temp_CV_alltime_m.nc4 temp_CV_total_m.nc4 # total Conv wet deposition 
cdo merge temp_LS_total_m.nc4 temp_CV_total_m.nc4 temp_both_total_m.nc4
cdo expr,'WetLossTot_Hg=LS_TOT+CV_TOT' temp_both_total_m.nc4 temp_sum_m1.nc4
cdo expr,'WetLossTot_HgP=WetLossLS_HG2ORGP+WetLossLS_HG2CLP+WetLossConv_HG2ORGP+WetLossConv_HG2CLP' temp_both_total_m.nc4 temp_sum_m_P.nc4
cdo expr,'WetLossTot_Hg2=WetLossLS_HGCL2+WetLossLS_HGOHOH+WetLossLS_HGOHBRO+WetLossLS_HGOHCLO+WetLossLS_HGOHHO2+WetLossLS_HGOHNO2+WetLossLS_HGCLOH+WetLossLS_HGCLBR+WetLossLS_HGCLBRO+WetLossLS_HGCLCLO+WetLossLS_HGCLHO2+WetLossLS_HGCLNO2+WetLossLS_HGBR2+WetLossLS_HGBROH+WetLossLS_HGBRCLO+WetLossLS_HGBRBRO+WetLossLS_HGBRHO2+WetLossLS_HGBRNO2+WetLossConv_HGCL2+WetLossConv_HGOHOH+WetLossConv_HGOHBRO+WetLossConv_HGOHCLO+WetLossConv_HGOHHO2+WetLossConv_HGOHNO2+WetLossConv_HGCLOH+WetLossConv_HGCLBR+WetLossConv_HGCLBRO+WetLossConv_HGCLCLO+WetLossConv_HGCLHO2+WetLossConv_HGCLNO2+WetLossConv_HGBR2+WetLossConv_HGBROH+WetLossConv_HGBRCLO+WetLossConv_HGBRBRO+WetLossConv_HGBRHO2+WetLossConv_HGBRNO2' temp_both_total_m.nc4 temp_sum_m_2.nc4

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
cdo expr,'DryDep_Hg2=DryDep_HGCL2+DryDep_HGOHOH+DryDep_HGOHBRO+DryDep_HGOHCLO+DryDep_HGOHHO2+DryDep_HGOHNO2+DryDep_HGCLOH+DryDep_HGCLBR+DryDep_HGCLBRO+DryDep_HGCLCLO+DryDep_HGCLHO2+DryDep_HGCLNO2+DryDep_HGBR2+DryDep_HGBROH+DryDep_HGBRCLO+DryDep_HGBRBRO+DryDep_HGBRHO2+DryDep_HGBRNO2' temp1.nc4 temp_Hg2_dd.nc4 # select Hg2 variable
cdo expr,'DryDep_HgP=DryDep_HG2ORGP+DryDep_HG2CLP' temp1.nc4 temp_HgP_dd.nc4 # select HgP variable

cdo merge temp_Hg0_dd.nc4 temp_Hg2_dd.nc4 temp_HgP_dd.nc4 GEOSChem.DryDep.alltime_m.nc4
#--------------------------------------------
# Monthly averages of mercury emissions/ocean files 
#---------------------------------------------
cdo mergetime GEOSChem.MercuryEmis.*_0000z.nc4 temp_GEOSChem.MercuryEmis.alltime.nc4
cdo mergetime GEOSChem.MercuryOcean.*_0000z.nc4 temp_GEOSChem.MercuryOcean.alltime.nc4

cdo monmean temp_GEOSChem.MercuryEmis.alltime.nc4 GEOSChem.MercuryEmis.alltime_m.nc4
cdo monmean temp_GEOSChem.MercuryOcean.alltime.nc4 GEOSChem.MercuryOcean.alltime_m.nc4
#---------------------------------------------
# Edits to ProdLoss files 
#---------------------------------------------
cdo mergetime GEOSChem.ProdLoss.*_0000z.nc4 temp_GEOSChem.ProdLoss.alltime.nc4

mv temp_GEOSChem.ProdLoss.alltime.nc4 temp1.nc4
cdo expr,'GroHg0_HgI=Prod_M0pBtM1+Prod_M0pOtM1' temp1.nc4 temp_Gro_Hg0_HgI.nc4 # select variable
cdo expr,'GroHgI_Hg0=Prod_MBtM0+Prod_MBpBtM0+Prod_MBpHVtM0+Prod_MOtM0+Prod_MOpHVtM0' temp1.nc4 temp_Gro_HgI_Hg0.nc4 # select variable
cdo expr,'GroHgI_HgII=Prod_MBpO3tM2+Prod_MOpO3tM2' temp1.nc4 temp_Gro_HgI_HgII.nc4 # select variable
cdo expr,'GroHgII_Hg0=Prod_MOOHpHVtM0+Prod_MBOHpHVtM0+Prod_MORGpHVtM0' temp1.nc4 temp_Gro_HgII_Hg0.nc4 # select variable
cdo expr,'GroHgIIp_Hg0=Prod_MORGpHVtM0' temp1.nc4 temp_Gro_HgIIp_Hg0.nc4 # select variable
cdo expr,'GroHgII_HgI=Prod_MBOpCOtMB+Prod_MBOHpHVtMB+Prod_MBOHpHVtMO+Prod_MOOpCOtMO' temp1.nc4 temp_Gro_HgII_HgI.nc4 # select variable

cdo merge temp_Gro_Hg*.nc4 temp_GEOSChem.ProdLoss.alltime_m.nc4
#  mv GEOSChem.ProdLoss.alltime_m.nc4 temp_GEOSChem.ProdLoss.alltime_m.nc4
cdo aexpr,'NetOx=GroHg0_HgI-GroHgI_Hg0-GroHgII_Hg0' temp_GEOSChem.ProdLoss.alltime_m.nc4 GEOSChem.ProdLoss.alltime_m.nc4
#---------------------------------------------
# Editing mercury chemistry files
#---------------------------------------------
cdo mergetime GEOSChem.MercuryChem.*_0000z.nc4 temp_GEOSChem.MercuryChem.alltime.nc4
cdo monmean temp_GEOSChem.MercuryChem.alltime.nc4 temp_GEOSChem.MercuryChem.alltime_m.nc4
cdo selvar,Hg2GasToSSA temp_GEOSChem.MercuryChem.alltime_m.nc4 temp_seasalt.nc4

# convert units from molec/cm3/s to kg/s
cdo selvar,Met_AIRVOL /net/fs03/d0/arifein/GEOS-Chem_runs/run0009/OutputDir/GEOSChem.StateMet.alltime_m.nc4 airvol.nc4
ncks -A -v Met_AIRVOL airvol.nc4 temp_seasalt.nc4
cdo expr,'LossHg2bySeaSalt=Hg2GasToSSA*Met_AIRVOL' temp_seasalt.nc4 temp_seasalt2.nc4
# unit conv  = 200.59 (g mol^-1) /1000 (g kg^-1) /6.02*10^23 (molec mol^-1) * 10^6 (m^3/cm^3)= 3.332e-19
cdo mulc,3.332059800664452e-19 temp_seasalt2.nc4 temp_seasalt3.nc4

cdo vertsum temp_seasalt3.nc4 temp_seasalt4.nc4

ncrename -v LossHg2bySeaSalt,LossHg2bySeaSalt_v temp_seasalt4.nc4
ncks -A -v LossHg2bySeaSalt_v temp_seasalt4.nc4 GEOSChem.MercuryChem.alltime_m.nc4
ncks -A -v LossHg2bySeaSalt temp_seasalt3.nc4 GEOSChem.MercuryChem.alltime_m.nc4

#add chemistry fluxes
ncks -A -v GroHg0_HgI GEOSChem.ProdLoss.alltime_m.nc4 temp_chem.nc4
ncks -A -v NetOx GEOSChem.ProdLoss.alltime_m.nc4 temp_chem.nc4
ncks -A -v Met_AIRVOL airvol.nc4 temp_chem.nc4

cdo expr,'Gross_Hg_Ox=GroHg0_HgI*Met_AIRVOL' temp_chem.nc4 temp_chem1.nc4
cdo mulc,3.332059800664452e-19 temp_chem1.nc4 temp_chem1u.nc4

cdo expr,'ProdHg2fromHg0=NetOx*Met_AIRVOL' temp_chem.nc4 temp_chem2.nc4
cdo mulc,3.332059800664452e-19 temp_chem2.nc4 temp_chem2u.nc4

ncks -A -v Gross_Hg_Ox temp_chem1u.nc4 GEOSChem.MercuryChem.alltime_m.nc4
ncks -A -v ProdHg2fromHg0 temp_chem2u.nc4 GEOSChem.MercuryChem.alltime_m.nc4
#---------------------------------------------
# Create mercury budget files
#---------------------------------------------
#---calculate global total emissions from different sources (kg/s)---
cdo -O -L -fldsum -selvar,EmisHg0anthro,EmisHg2HgPanthro,EmisHg0geogenic,EmisHg0soil,EmisHg0biomass,EmisHg0land,EmisHg0snow GEOSChem.MercuryEmis.alltime_m.nc4 temp_emiss.nc4
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
