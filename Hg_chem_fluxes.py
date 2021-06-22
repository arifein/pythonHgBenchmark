# -*- coding: utf-8 -*-
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from helper_functions import ds_sel_yr, annual_avg
from diff_plots_Hg import diff_plots

def chem_plots(Dataset_OLD, Dataset_NEW, Year = None):
    """Main script for calling different routines that produce Hg chemistry plots
    
    Parameters
    ----------
    Dataset_OLD : xarray dataset
        Reference Model dataset (Hg chemical fluxes)
                
    Dataset_NEW : xarray dataset
        New Model dataset (Hg chemical fluxes)
            
    Year : int or list of int, optional
        Optional parameter to only select subset of years    
    
    """
    # Sea salt uptake of Hg(II)
    
    # Allow subsetting for years, if inputted into the function
    OLD_Hg2_salt_yr = ds_sel_yr(Dataset_OLD, 'LossHg2bySeaSalt_v', Year) 
    NEW_Hg2_salt_yr = ds_sel_yr(Dataset_NEW, 'LossHg2bySeaSalt_v', Year)
    
    
    # calculate annual average
    OLD_Hg2_salt = annual_avg(OLD_Hg2_salt_yr)
    NEW_Hg2_salt = annual_avg(NEW_Hg2_salt_yr)
    
    # Convert model data from kg/s to ug/m^2/yr for annual average   
    # Load grid cell area for unit conversion of model
    fn_gbox = 'data/GEOSChem_2x25_gboxarea.nc'
    ds_gbox = xr.open_dataset(fn_gbox)
    gbox_GC = ds_gbox.cell_area
    
    s_in_yr = 365.2425 * 24 * 3600 # s in one year
    kg_ug = 1e9 # kg in ug 
    
    unit_conv = s_in_yr * kg_ug / gbox_GC
            
    OLD_Hg2_salt = OLD_Hg2_salt * unit_conv # ug/m^2/yr
    NEW_Hg2_salt = NEW_Hg2_salt * unit_conv # ug/m^2/yr
        
    # Plot Hg0 dry deposition difference plot 
    plot1 = diff_plots(OLD_Hg2_salt, NEW_Hg2_salt, 
                       Units="\u03BCg m$^{-2}$ yr$^{-1}$ ",
                       Title="Sea Salt Uptake")
    
    return plot1