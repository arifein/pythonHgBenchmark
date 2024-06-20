#!/bin/sh

# Postprocess CESM2 WACCM before benchmark simulation 
module load nco
module load cdo

# Simulation name
SIM_NAME=FW.f19_f19_mg17.L70.cesm2.2-asdbranch_vsl03_CSIC_Hg_v2.Testing_HG_full.FWnudged_slh
DIR_NAME=/glade/derecho/scratch/afeinberg/Setup_Derecho/CESM2-SLH_Madrid24/Testing_CESM2-SLH/archive/${SIM_NAME}/atm/hist/

cd ${DIR_NAME}

#---------------------------------------------
# Select variables that will need for the benchmark 
#---------------------------------------------
cdo select,name=AREA,SFHG,SFHGCL2,SFHGP,HGBR2,HGBRO,HGBRNO2,HGBROH,HGBROOH,HGCL2,HGO_S3,HGO_S1,HGBR,HGCL,HGOH,HG,HG3P1,HG3P0,HGP,HGBRCL,HGCLO,HGOHOH,HGCLOH,HGOHONO,HGOHOOH,HGOHO,HGCLONO,HGCLOOH,DF_HGBR2,DF_HGBRNO2,DF_HGBROH,DF_HGBROOH,DF_HGCL2,DF_HGP,DF_HGBRCL,DF_HGCLO,DF_HGOHOH,DF_HGCLOH,DF_HGOHONO,DF_HGOHOOH,DF_HGOHO,DF_HGCLONO,DF_HGCLOOH,DF_HG,DV_HGP,DV_HGBR2,DV_HGCL2,DV_HGBROH,DV_HGOHOH,DV_HG,WD_HGBR2,WD_HGBRNO2,WD_HGBROH,WD_HGBROOH,WD_HGCL2,WD_HGO_S3,WD_HGO_S1,WD_HGP,WD_HGBRCL,WD_HGCLO,WD_HGOHOH,WD_HGCLOH,WD_HGOHONO,WD_HGOHOOH,WD_HGOHO,WD_HGCLONO,WD_HGCLOOH,DHGBR2CHM,CT_HGBR2,DHGBROCHM,CT_HGBRO,DHGBRNO2CHM,CT_HGBRNO2,DHGBROHCHM,CT_HGBROH,DHGBROOHCHM,CT_HGBROOH,DHGCL2CHM,CT_HGCL2,DHGO_S3CHM,CT_HGO_S3,DHGO_S1CHM,CT_HGO_S1,DHGBRCHM,CT_HGBR,DHGCLCHM,CT_HGCL,DHGOHCHM,CT_HGOH,DHGCHM,CT_HG,DHG3P1CHM,CT_HG3P1,DHG3P0CHM,CT_HG3P0,DHGPCHM,CT_HGP,DHGBRCLCHM,CT_HGBRCL,DHGCLOCHM,CT_HGCLO,DHGOHOHCHM,CT_HGOHOH,DHGCLOHCHM,CT_HGCLOH,DHGOHONOCHM,CT_HGOHONO,DHGOHOOHCHM,CT_HGOHOOH,DHGOHOCHM,CT_HGOHO,DHGCLONOCHM,CT_HGCLONO,DHGCLOOHCHM,CT_HGCLOOH ${SIM_NAME}.cam.h0.????-??.nc  temp_${SIM_NAME}.cam.h0_alltime_m_sel.nc

#---------------------------------------------
# Summarized concentrations for comparison
#---------------------------------------------
# Calculate total oxidized mercury concentration
cdo expr,"vmrhg2=HGBR2+HGBRO+HGBRNO2+HGBROH+HGBROOH+HGCL2+HGO_S3+HGO_S1+HGBR+HGCL+HGOH+HG3P1+HG3P0+HGBRCL+HGCLO+HGOHOH+HGCLOH+HGOHONO+HGOHOOH+HGOHO+HGCLONO+HGCLOOH" temp_${SIM_NAME}.cam.h0_alltime_m_sel.nc temp_${SIM_NAME}.SpeciesConc_alltime_m_hg2.nc
ncatted -h -O -a  long_name,vmrhg2,o,c,"HgII(g) Volume Mixing Ratio" temp_${SIM_NAME}.SpeciesConc_alltime_m_hg2.nc
ncatted -h -O -a  units,vmrhg2,o,c,"mol/mol" temp_${SIM_NAME}.SpeciesConc_alltime_m_hg2.nc

cdo expr,"vmrhg0=HG" temp_${SIM_NAME}.cam.h0_alltime_m_sel.nc temp_${SIM_NAME}.SpeciesConc_alltime_m_hg0.nc
cdo expr,"vmrhgp=HGP" temp_${SIM_NAME}.cam.h0_alltime_m_sel.nc temp_${SIM_NAME}.SpeciesConc_alltime_m_hgp.nc

cdo merge temp_${SIM_NAME}.SpeciesConc_alltime_m_hg2.nc temp_${SIM_NAME}.SpeciesConc_alltime_m_hg0.nc temp_${SIM_NAME}.SpeciesConc_alltime_m_hgp.nc temp_${SIM_NAME}.SpeciesConc_alltime_m.nc
nccopy -d1 temp_${SIM_NAME}.SpeciesConc_alltime_m.nc ${SIM_NAME}.SpeciesConc_alltime_m.nc

#---------------------------------------------
# Summarized wet deposition comparison 
#---------------------------------------------
# sum up individual files
# need to account for different molar mass of Hg2 species
mw_HG=200.59
mw_HGBR2=360.398
mw_HGBRNO2=326.4995
mw_HGBROH=297.50134
mw_HGBROOH=313.50074
mw_HGCL2=271.496
mw_HGO_S3=216.5894
mw_HGO_S1=216.5894
mw_HGP=200.59
mw_HGBRCL=315.947
mw_HGCLO=252.0424
mw_HGOHOH=234.60468
mw_HGCLOH=253.05034
mw_HGOHONO=263.60284
mw_HGOHOOH=250.60408
mw_HGOHO=233.59674
mw_HGCLONO=282.0485
mw_HGCLOOH=269.0497
# account for difference in molar mass, account for the fact that have negative values 
cdo expr,"wethg2=-(WD_HGBR2*${mw_HG}/${mw_HGBR2}+WD_HGBRNO2*${mw_HG}/${mw_HGBRNO2}+WD_HGBROH*${mw_HG}/${mw_HGBROH}+WD_HGBROOH*${mw_HG}/${mw_HGBROOH}+WD_HGCL2*${mw_HG}/${mw_HGCL2}+WD_HGO_S3*${mw_HG}/${mw_HGO_S3}+WD_HGO_S1*${mw_HG}/${mw_HGO_S1}+WD_HGBRCL*${mw_HG}/${mw_HGBRCL}+WD_HGCLO*${mw_HG}/${mw_HGCLO}+WD_HGOHOH*${mw_HG}/${mw_HGOHOH}+WD_HGCLOH*${mw_HG}/${mw_HGCLOH}+WD_HGOHONO*${mw_HG}/${mw_HGOHONO}+WD_HGOHOOH*${mw_HG}/${mw_HGOHOOH}+WD_HGOHO*${mw_HG}/${mw_HGOHO}+WD_HGCLONO*${mw_HG}/${mw_HGCLONO}+WD_HGCLOOH*${mw_HG}/${mw_HGCLOOH})" temp_${SIM_NAME}.cam.h0_alltime_m_sel.nc temp_${SIM_NAME}.WetDep_alltime_m_hg2.nc
cdo expr,'wethgp=-WD_HGP' temp_${SIM_NAME}.cam.h0_alltime_m_sel.nc temp_${SIM_NAME}.WetDep_alltime_m_hgp.nc

cdo merge temp_${SIM_NAME}.WetDep_alltime_m_hg2.nc temp_${SIM_NAME}.WetDep_alltime_m_hgp.nc temp_${SIM_NAME}.WetDep_alltime_m_0.nc 

# select only positive values to avoid issues (minor part of wet deposition flux globally)
cdo setrtoc,-1.e99,0,0 temp_${SIM_NAME}.WetDep_alltime_m_0.nc temp_${SIM_NAME}.WetDep_alltime_m_noneg.nc 

cdo aexpr,'wethg=wethgp+wethg2' temp_${SIM_NAME}.WetDep_alltime_m_noneg.nc temp_${SIM_NAME}.WetDep_alltime_m.nc 
 
# edit metadata
ncatted -h -O -a  long_name,wethg,o,c,"Wet deposition of Hg" temp_${SIM_NAME}.WetDep_alltime_m.nc
ncatted -h -O -a  units,wethg,o,c,"kg Hg m-2 s-1" temp_${SIM_NAME}.WetDep_alltime_m.nc
ncatted -h -O -a  long_name,wethg2,o,c,"Wet deposition of Hg2" temp_${SIM_NAME}.WetDep_alltime_m.nc
ncatted -h -O -a  units,wethg2,o,c,"kg Hg m-2 s-1" temp_${SIM_NAME}.WetDep_alltime_m.nc
ncatted -h -O -a  long_name,wethgp,o,c,"Wet deposition of HgP" temp_${SIM_NAME}.WetDep_alltime_m.nc
ncatted -h -O -a  units,wethgp,o,c,"kg Hg m-2 s-1" temp_${SIM_NAME}.WetDep_alltime_m.nc

nccopy -d1 temp_${SIM_NAME}.WetDep_alltime_m.nc ${SIM_NAME}.WetDep_alltime_m.nc

#---------------------------------------------
# Summarized dry deposition comparison 
#---------------------------------------------
# account for difference in molar mass, account for the fact that have negative values 
cdo selvar,DV_HG,DV_HGCL2,DV_HGP temp_${SIM_NAME}.cam.h0_alltime_m_sel.nc temp_${SIM_NAME}.DryDep.h0_alltime_m_dv.nc
cdo expr,'dryhg0=DF_HG' temp_${SIM_NAME}.cam.h0_alltime_m_sel.nc temp_${SIM_NAME}.DryDep_alltime_m_hg0.nc
cdo expr,'dryhgp=DF_HGP' temp_${SIM_NAME}.cam.h0_alltime_m_sel.nc temp_${SIM_NAME}.DryDep_alltime_m_hgp.nc
cdo expr,"dryhg2=DF_HGBR2*${mw_HG}/${mw_HGBR2}+DF_HGBRNO2*${mw_HG}/${mw_HGBRNO2}+DF_HGBROH*${mw_HG}/${mw_HGBROH}+DF_HGBROOH*${mw_HG}/${mw_HGBROOH}+DF_HGCL2*${mw_HG}/${mw_HGCL2}+DF_HGBRCL*${mw_HG}/${mw_HGBRCL}+DF_HGCLO*${mw_HG}/${mw_HGCLO}+DF_HGOHOH*${mw_HG}/${mw_HGOHOH}+DF_HGCLOH*${mw_HG}/${mw_HGCLOH}+DF_HGOHONO*${mw_HG}/${mw_HGOHONO}+DF_HGOHOOH*${mw_HG}/${mw_HGOHOOH}+DF_HGOHO*${mw_HG}/${mw_HGOHO}+DF_HGCLONO*${mw_HG}/${mw_HGCLONO}+DF_HGCLOOH*${mw_HG}/${mw_HGCLOOH}" temp_${SIM_NAME}.cam.h0_alltime_m_sel.nc temp_${SIM_NAME}.DryDep_alltime_m_hg2.nc

cdo merge temp_${SIM_NAME}.DryDep.h0_alltime_m_dv.nc temp_${SIM_NAME}.DryDep_alltime_m_hg0.nc temp_${SIM_NAME}.DryDep_alltime_m_hg2.nc temp_${SIM_NAME}.DryDep_alltime_m_hgp.nc temp_${SIM_NAME}.DryDep_alltime_m.nc 
 
# edit metadata
ncatted -h -O -a  long_name,dryhg0,o,c,"Dry deposition of Hg0(g) to land" temp_${SIM_NAME}.DryDep_alltime_m.nc
ncatted -h -O -a  units,dryhg0,o,c,"kg Hg m-2 s-1" temp_${SIM_NAME}.DryDep_alltime_m.nc
ncatted -h -O -a  long_name,dryhg2,o,c,"Dry deposition of HgII(g)" temp_${SIM_NAME}.DryDep_alltime_m.nc
ncatted -h -O -a  units,dryhg2,o,c,"kg Hg m-2 s-1" temp_${SIM_NAME}.DryDep_alltime_m.nc
ncatted -h -O -a  long_name,dryhgp,o,c,"Dry deposition of particulate mercury" temp_${SIM_NAME}.DryDep_alltime_m.nc
ncatted -h -O -a  units,dryhgp,o,c,"kg Hg m-2 s-1" temp_${SIM_NAME}.DryDep_alltime_m.nc

nccopy -d1 temp_${SIM_NAME}.DryDep_alltime_m.nc ${SIM_NAME}.DryDep_alltime_m.nc

#--------------------------------------------
# Monthly averages of mercury emissions/ocean files 
#---------------------------------------------
cdo selvar,SFHG,SFHGCL2,SFHGP temp_${SIM_NAME}.cam.h0_alltime_m_sel.nc temp_${SIM_NAME}.MercuryEmis_alltime_m_SF.nc
cdo select,name=OCE_FLUX_HG  ${SIM_NAME}.cam.h5.????-??.nc temp_${SIM_NAME}.MercuryEmis_alltime_m_oce.nc 
# adjust units
cdo aexpr,"SFHGCL2=SFHGCL2*${mw_HG}/${mw_HGCL2}" temp_${SIM_NAME}.MercuryEmis_alltime_m_SF.nc temp_${SIM_NAME}.MercuryEmis_alltime_m_SF_u.nc
# edit units
ncatted -h -O -a  units,SFHGCL2,o,c,"kg Hg m-2 s-1" temp_${SIM_NAME}.MercuryEmis_alltime_m_SF_u.nc
# merge files
cdo merge temp_${SIM_NAME}.MercuryEmis_alltime_m_SF_u.nc temp_${SIM_NAME}.MercuryEmis_alltime_m_oce.nc temp_${SIM_NAME}.MercuryEmis_alltime_m.nc 
# compress files
nccopy -d1 temp_${SIM_NAME}.MercuryEmis_alltime_m.nc ${SIM_NAME}.MercuryEmis_alltime_m.nc


# #---------------------------------------------
# # Edits to ProdLoss files 
# #---------------------------------------------
# # 
# #cdo aexpr,'Gross_Hg_Ox=Prod_Hg2Br+Prod_Hg2OH+Prod_Hg2Cl' temp1_pl.nc temp2_pl.nc # Gross GOM oxidation # 0205
# cdo aexpr,'Gross_Hg_Ox=Prod_Hg2' temp1_pl.nc temp2_pl.nc # Gross GOM oxidation 
# cdo aexpr,'NetOx=Gross_Hg_Ox-Prod_Hg0' temp2_pl.nc GEOSChem.ProdLoss.alltime_m.nc # Net GOM oxidation

#---------------------------------------------
# Create mercury budget files
#---------------------------------------------
#---calculate global total emissions from different sources (kg Hg/s)---
cdo -O -L -fldsum -mul -selvar,SFHG,SFHGCL2,SFHGP,OCE_FLUX_HG ${SIM_NAME}.MercuryEmis_alltime_m.nc -gridarea ${SIM_NAME}.MercuryEmis_alltime_m.nc temp_emiss.nc


#---calculate global total wet deposition (kg/s)--------------------
cdo -O -L -fldsum -mul -selvar,wethgp,wethg2 ${SIM_NAME}.WetDep_alltime_m.nc -gridarea ${SIM_NAME}.WetDep_alltime_m.nc temp_wetdep.nc

#---calculate global total dry deposition (kg/s)--------------------
cdo -O -L -fldsum -mul -selvar,dryhg0,dryhgp,dryhg2 ${SIM_NAME}.DryDep_alltime_m.nc -gridarea ${SIM_NAME}.DryDep_alltime_m.nc temp_drydep.nc

# #---calculate global total chemistry fluxes (kg/s)------------------
# cdo -O -L -fldsum -vertsum -selvar,Gross_Hg_Ox,ProdHg2fromHg0 GEOSChem.MercuryChem.alltime_m.nc temp_chem1.nc
# cdo -O -L -fldsum -selvar,LossHg2bySeaSalt_v GEOSChem.MercuryChem.alltime_m.nc temp_chem2.nc
# ncrename -v LossHg2bySeaSalt_v,LossHg2bySeaSalt temp_chem2.nc
# cdo merge temp_chem1.nc temp_chem2.nc temp_chem3.nc
#------------------------------------------------------------------


#---merge files, clean up metadata, remove temporary files---------
# merge intermediate files
cdo merge temp_emiss.nc temp_wetdep.nc temp_drydep.nc ${SIM_NAME}.MercuryBudget_global_m.nc

# change units
ncatted -h -O -a  units,wethgp,o,c,"kg s-1" ${SIM_NAME}.MercuryBudget_global_m.nc
ncatted -h -O -a  units,wethg2,o,c,"kg s-1" ${SIM_NAME}.MercuryBudget_global_m.nc
ncatted -h -O -a  units,dryhg0,o,c,"kg s-1" ${SIM_NAME}.MercuryBudget_global_m.nc
ncatted -h -O -a  units,dryhgp,o,c,"kg s-1" ${SIM_NAME}.MercuryBudget_global_m.nc
ncatted -h -O -a  units,dryhg2,o,c,"kg s-1" ${SIM_NAME}.MercuryBudget_global_m.nc
ncatted -h -O -a  units,SFHG,o,c,"kg s-1" ${SIM_NAME}.MercuryBudget_global_m.nc
ncatted -h -O -a  units,SFHGCL2,o,c,"kg Hg s-1" ${SIM_NAME}.MercuryBudget_global_m.nc
ncatted -h -O -a  units,SFHGP,o,c,"kg s-1" ${SIM_NAME}.MercuryBudget_global_m.nc
ncatted -h -O -a  units,OCE_FLUX_HG,o,c,"kg s-1" ${SIM_NAME}.MercuryBudget_global_m.nc

#---------------------------------------------
# remove temporary files
#---------------------------------------------
rm temp_*.nc*
