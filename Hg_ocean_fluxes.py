# -*- coding: utf-8 -*-
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from helper_functions import ds_sel_yr, annual_avg
from diff_plots_Hg import diff_plots

def ocean_plots(Dataset_OLD, Dataset_NEW, Year = None):
    """Main script for calling different routines that produce Hg ocean plots
    
    Parameters
    ----------
    Dataset_OLD : xarray dataset
        Reference Model dataset (Hg ocean fluxes)
                
    Dataset_NEW : xarray dataset
        New Model dataset (Hg ocean fluxes)
            
    Year : int or list of int, optional
        Optional parameter to only select subset of years    
    
    """
    #---Gross Ocean Evasion---
    
    # Allow subsetting for years, if inputted into the function
    OLD_gross_evas_yr = ds_sel_yr(Dataset_OLD, 'FluxHg0fromOceanToAir', Year) 
    NEW_gross_evas_yr = ds_sel_yr(Dataset_NEW, 'FluxHg0fromOceanToAir', Year)
    
    
    # calculate annual average
    OLD_gross_evas = annual_avg(OLD_gross_evas_yr)
    NEW_gross_evas = annual_avg(NEW_gross_evas_yr)
    
    # Convert model data from kg/s to kg/yr for annual average       
    s_in_yr = 365.2425 * 24 * 3600 # s in one year
    
    unit_conv = s_in_yr 
            
    OLD_gross_evas = OLD_gross_evas * unit_conv # kg/yr
    NEW_gross_evas = NEW_gross_evas * unit_conv # kg/yr
        
    # Plot Hg0 emissions difference plot 
    plot1 = diff_plots(OLD_gross_evas, NEW_gross_evas, 
                       Units="kg yr$^{-1}$ ",
                       Title="Gross Ocean Evasion")

    #---Gross Ocean Hg(0) Uptake---
    
    # Allow subsetting for years, if inputted into the function
    OLD_gross_uptake_yr = ds_sel_yr(Dataset_OLD, 'FluxHg0fromAirToOcean', Year) 
    NEW_gross_uptake_yr = ds_sel_yr(Dataset_NEW, 'FluxHg0fromAirToOcean', Year)
    
    
    # calculate annual average
    OLD_gross_uptake = annual_avg(OLD_gross_uptake_yr)
    NEW_gross_uptake = annual_avg(NEW_gross_uptake_yr)
    
    # Convert model data from kg/s to kg/yr for annual average       
    s_in_yr = 365.2425 * 24 * 3600 # s in one year
    
    unit_conv = s_in_yr 
            
    OLD_gross_uptake = OLD_gross_uptake * unit_conv # kg/yr
    NEW_gross_uptake = NEW_gross_uptake * unit_conv # kg/yr
        
    # Plot Hg0 emissions difference plot 
    plot2 = diff_plots(OLD_gross_uptake, NEW_gross_uptake, 
                       Units="kg yr$^{-1}$ ",
                       Title="Gross Ocean Hg(0) Uptake")
    
    #---Net Ocean Evasion---
    
    OLD_net_evas = OLD_gross_evas - OLD_gross_uptake # kg/yr
    NEW_net_evas = NEW_gross_evas - NEW_gross_uptake # kg/yr
        
    # Plot Hg0 emissions difference plot 
    plot3 = diff_plots(OLD_net_evas, NEW_net_evas, 
                       Units="kg yr$^{-1}$ ",
                       Title="Net Ocean Evasion")

    plotlist = [plot1, plot2, plot3]
    return plotlist