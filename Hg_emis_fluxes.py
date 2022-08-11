# -*- coding: utf-8 -*-
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from helper_functions import ds_sel_yr, annual_avg
from diff_plots_Hg import diff_plots

def emis_plots(Dataset_OLD, Dataset_NEW, Year1 = None, Year2 = None):
    """Main script for calling different routines that produce Hg emission plots
    
    Parameters
    ----------
    Dataset_OLD : xarray dataset
        Reference Model dataset (Hg emission fluxes)
                
    Dataset_NEW : xarray dataset
        New Model dataset (Hg emission fluxes)
            
    Year1 : int or list of int, optional
        Optional parameter to only select subset of years for old sim
    Year2 : int or list of int, optional
        Optional parameter to only select subset of years for new sim
    
    """
    #---Anthropogenic emissions of Hg0---
    
    # Allow subsetting for years, if inputted into the function
    OLD_Hg0_emis_yr = ds_sel_yr(Dataset_OLD, 'EmisHg0anthro', Year1) 
    NEW_Hg0_emis_yr = ds_sel_yr(Dataset_NEW, 'EmisHg0anthro', Year2)
    
    
    # calculate annual average
    OLD_Hg0_emis = annual_avg(OLD_Hg0_emis_yr)
    NEW_Hg0_emis = annual_avg(NEW_Hg0_emis_yr)
    
    # Convert model data from kg/s to kg/yr for annual average       
    s_in_yr = 365.2425 * 24 * 3600 # s in one year
    
    unit_conv = s_in_yr 
            
    OLD_Hg0_emis = OLD_Hg0_emis * unit_conv # kg/yr
    NEW_Hg0_emis = NEW_Hg0_emis * unit_conv # kg/yr
        
    # Plot Hg0 emissions difference plot 
    plot1 = diff_plots(OLD_Hg0_emis, NEW_Hg0_emis, 
                       Units="kg yr$^{-1}$ ",
                       Title="Anthro Emissions - Hg(0)")

    #---Anthropogenic emissions of HgII + HgP---
    
    # Allow subsetting for years, if inputted into the function
    OLD_Hg2_emis_yr = ds_sel_yr(Dataset_OLD, 'EmisHg2HgPanthro', Year1) 
    NEW_Hg2_emis_yr = ds_sel_yr(Dataset_NEW, 'EmisHg2HgPanthro', Year2)
    
    
    # calculate annual average
    OLD_Hg2_emis = annual_avg(OLD_Hg2_emis_yr)
    NEW_Hg2_emis = annual_avg(NEW_Hg2_emis_yr)
    
    # Convert model data from kg/s to kg/yr for annual average       
            
    OLD_Hg2_emis = OLD_Hg2_emis * unit_conv # kg/yr
    NEW_Hg2_emis = NEW_Hg2_emis * unit_conv # kg/yr
        
    # Plot Hg2 emissions difference plot 
    plot2 = diff_plots(OLD_Hg2_emis, NEW_Hg2_emis, 
                       Units="kg yr$^{-1}$ ",
                       Title="Anthro Emissions - Hg(II)+Hg(P)")

    #---Terrestrial emissions - GEO, BB, & Soil---
    
    # Allow subsetting for years, if inputted into the function
    # geogenic
    OLD_geo_emis_yr = ds_sel_yr(Dataset_OLD, 'EmisHg0geogenic', Year1) 
    NEW_geo_emis_yr = ds_sel_yr(Dataset_NEW, 'EmisHg0geogenic', Year2)
    
    # biomass burning
    OLD_bb_emis_yr = ds_sel_yr(Dataset_OLD, 'EmisHg0biomass', Year1) 
    NEW_bb_emis_yr = ds_sel_yr(Dataset_NEW, 'EmisHg0biomass', Year2)

    # soil
    OLD_soil_emis_yr = ds_sel_yr(Dataset_OLD, 'EmisHg0soil', Year1) 
    NEW_soil_emis_yr = ds_sel_yr(Dataset_NEW, 'EmisHg0soil', Year2)
    
    # calculate annual averages
    OLD_geo_emis = annual_avg(OLD_geo_emis_yr)
    NEW_geo_emis = annual_avg(NEW_geo_emis_yr)

    OLD_bb_emis = annual_avg(OLD_bb_emis_yr)
    NEW_bb_emis = annual_avg(NEW_bb_emis_yr)

    OLD_soil_emis = annual_avg(OLD_soil_emis_yr)
    NEW_soil_emis = annual_avg(NEW_soil_emis_yr)
    
    # sum of emissions
    OLD_terr_emis = OLD_geo_emis + OLD_bb_emis + OLD_soil_emis
    NEW_terr_emis = NEW_geo_emis + NEW_bb_emis + NEW_soil_emis
    
    # Convert model data from kg/s to kg/yr for annual average       
            
    OLD_terr_emis = OLD_terr_emis * unit_conv # kg/yr
    NEW_terr_emis = NEW_terr_emis * unit_conv # kg/yr
        
    # Plot Hg0 terrestrial emissions difference plot 
    plot3 = diff_plots(OLD_terr_emis, NEW_terr_emis, 
                       Units="kg yr$^{-1}$ ",
                       Title="Direct Terrestrial - Geo, BB, & Soil")

    #---Re-emission - Land and Snow---
    
    # Allow subsetting for years, if inputted into the function
    # land
    OLD_land_emis_yr = ds_sel_yr(Dataset_OLD, 'EmisHg0land', Year1) 
    NEW_land_emis_yr = ds_sel_yr(Dataset_NEW, 'EmisHg0land', Year2)
    
    # snow
    OLD_snow_emis_yr = ds_sel_yr(Dataset_OLD, 'EmisHg0snow', Year1) 
    NEW_snow_emis_yr = ds_sel_yr(Dataset_NEW, 'EmisHg0snow', Year2)
    
    # calculate annual averages
    OLD_land_emis = annual_avg(OLD_land_emis_yr)
    NEW_land_emis = annual_avg(NEW_land_emis_yr)

    OLD_snow_emis = annual_avg(OLD_snow_emis_yr)
    NEW_snow_emis = annual_avg(NEW_snow_emis_yr)
    
    # sum of emissions
    OLD_re_emis = OLD_land_emis + OLD_snow_emis
    NEW_re_emis = NEW_land_emis + NEW_snow_emis
    
    # Convert model data from kg/s to kg/yr for annual average       
            
    OLD_re_emis = OLD_re_emis * unit_conv # kg/yr
    NEW_re_emis = NEW_re_emis * unit_conv # kg/yr
        
    # Plot prompt re-emissions difference plot 
    plot4 = diff_plots(OLD_re_emis, NEW_re_emis, 
                       Units="kg yr$^{-1}$ ",
                       Title="Prompt Re-emission - Land & Snow")
    
    #---Re-emission - Ocean---
    
    # Allow subsetting for years, if inputted into the function
    # ocean
    OLD_oc_emis_yr = ds_sel_yr(Dataset_OLD, 'EmisHg0ocean', Year1) 
    NEW_oc_emis_yr = ds_sel_yr(Dataset_NEW, 'EmisHg0ocean', Year2)
    
    
    # calculate annual averages
    OLD_oc_emis = annual_avg(OLD_oc_emis_yr)
    NEW_oc_emis = annual_avg(NEW_oc_emis_yr)
    
    # Convert model data from kg/s to kg/yr for annual average       
            
    OLD_oc_emis = OLD_oc_emis * unit_conv # kg/yr
    NEW_oc_emis = NEW_oc_emis * unit_conv # kg/yr
        
    # Plot prompt re-emissions difference plot 
    plot5 = diff_plots(OLD_oc_emis, NEW_oc_emis, 
                       Units="kg yr$^{-1}$ ",
                       Title="Net Ocean Evasion")

    plotlist = [plot1, plot2, plot3, plot4, plot5]
    return plotlist