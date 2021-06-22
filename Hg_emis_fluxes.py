# -*- coding: utf-8 -*-
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from helper_functions import ds_sel_yr, annual_avg
from diff_plots_Hg import diff_plots

def emis_plots(Dataset_OLD, Dataset_NEW, Year = None):
    """Main script for calling different routines that produce Hg emission plots
    
    Parameters
    ----------
    Dataset_OLD : xarray dataset
        Reference Model dataset (Hg emission fluxes)
                
    Dataset_NEW : xarray dataset
        New Model dataset (Hg emission fluxes)
            
    Year : int or list of int, optional
        Optional parameter to only select subset of years    
    
    """
    # Anthropogenic emissions of Hg0
    
    # Allow subsetting for years, if inputted into the function
    OLD_Hg0_emis_yr = ds_sel_yr(Dataset_OLD, 'EmisHg0anthro', Year) 
    NEW_Hg0_emis_yr = ds_sel_yr(Dataset_NEW, 'EmisHg0anthro', Year)
    
    
    # calculate annual average
    OLD_Hg0_emis = annual_avg(OLD_Hg0_emis_yr)
    NEW_Hg0_emis = annual_avg(NEW_Hg0_emis_yr)
    
    # Convert model data from kg/s to kg/yr for annual average       
    s_in_yr = 365.2425 * 24 * 3600 # s in one year
    
    unit_conv = s_in_yr 
            
    OLD_Hg0_emis = OLD_Hg0_emis * unit_conv # kg/yr
    NEW_Hg0_emis = NEW_Hg0_emis * unit_conv # kg/yr
        
    # Plot Hg0 dry deposition difference plot 
    plot1 = diff_plots(OLD_Hg0_emis, NEW_Hg0_emis, 
                       Units="kg yr$^{-1}$ ",
                       Title="Anthro Emissions - Hg(0)")
    
    return plot1