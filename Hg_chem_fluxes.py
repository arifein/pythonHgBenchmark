# -*- coding: utf-8 -*-
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from helper_functions import ds_sel_yr, annual_avg
from diff_plots_Hg import diff_plots, diff_profiles

def chem_plots(Dataset_OLD, Dataset_NEW, Year1 = None, Year2 = None):
    """Main script for calling different routines that produce Hg chemistry plots
    
    Parameters
    ----------
    Dataset_OLD : xarray dataset
        Reference Model dataset (Hg chemical fluxes)
                
    Dataset_NEW : xarray dataset
        New Model dataset (Hg chemical fluxes)
            
    Year1 : int or list of int, optional
        Optional parameter to only select subset of years for old sim
    Year2 : int or list of int, optional
        Optional parameter to only select subset of years for new sim
    
    """
    #---Sea salt uptake of Hg(II)---
    
    # Allow subsetting for years, if inputted into the function
    OLD_Hg2_salt_yr = ds_sel_yr(Dataset_OLD, 'LossHg2bySeaSalt_v', Year1) 
    NEW_Hg2_salt_yr = ds_sel_yr(Dataset_NEW, 'LossHg2bySeaSalt_v', Year2)
    
    
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
        
    # Plot Hg2+ sea salt uptake difference plot 
    plot1 = diff_plots(OLD_Hg2_salt, NEW_Hg2_salt, 
                       Units="\u03BCg m$^{-2}$ yr$^{-1}$ ",
                       Title="Sea Salt Uptake")
    
    #---Zonal Gross Oxidation---
    
    # Allow subsetting for years, if inputted into the function
    OLD_gross_oxid_yr = ds_sel_yr(Dataset_OLD, 'Gross_Hg_Ox', Year1) 
    NEW_gross_oxid_yr = ds_sel_yr(Dataset_NEW, 'Gross_Hg_Ox', Year2)

    # calculate annual and zonal average
    OLD_gross_oxid = annual_avg(OLD_gross_oxid_yr.mean('lon'))
    NEW_gross_oxid = annual_avg(NEW_gross_oxid_yr.mean('lon'))

    # Convert model data from kg/s to kg/yr for annual average       
    s_in_yr = 365.2425 * 24 * 3600 # s in one year
    
    unit_conv = s_in_yr
    
    OLD_gross_oxid = OLD_gross_oxid * unit_conv # kg/yr
    NEW_gross_oxid = NEW_gross_oxid * unit_conv # kg/yr

    # Plot Hg gross oxidation difference plot 
    plot2 = diff_profiles(OLD_gross_oxid, NEW_gross_oxid, 
                       Units="kg yr$^{-1}$ ",
                       Title="Zonal Gross Oxidation")

    #---Zonal Net Oxidation---
    
    # Allow subsetting for years, if inputted into the function
    OLD_net_oxid_yr = ds_sel_yr(Dataset_OLD, 'ProdHg2fromHg0', Year1) 
    NEW_net_oxid_yr = ds_sel_yr(Dataset_NEW, 'ProdHg2fromHg0', Year2)

    # calculate annual and zonal average
    OLD_net_oxid = annual_avg(OLD_net_oxid_yr.mean('lon'))
    NEW_net_oxid = annual_avg(NEW_net_oxid_yr.mean('lon'))

    # Convert model data from kg/s to kg/yr for annual average       
    s_in_yr = 365.2425 * 24 * 3600 # s in one year
    
    unit_conv = s_in_yr
    
    OLD_net_oxid = OLD_net_oxid * unit_conv # kg/yr
    NEW_net_oxid = NEW_net_oxid * unit_conv # kg/yr

    # Plot Hg net oxidation difference plot 
    plot3 = diff_profiles(OLD_net_oxid, NEW_net_oxid, 
                       Units="kg yr$^{-1}$ ",
                       Title="Zonal Net Oxidation")

    #---Zonal Gross Reduction---
    # calculate gross reduction as gross minus net oxidation 
    OLD_gross_red = OLD_gross_oxid - OLD_net_oxid
    NEW_gross_red = NEW_gross_oxid - NEW_net_oxid
    
    # Plot Hg gross reduction difference plot 
    plot4 = diff_profiles(OLD_gross_red, NEW_gross_red, 
                       Units="kg yr$^{-1}$ ",
                       Title="Zonal Gross Reduction")
    
    plotlist = [plot1, plot2, plot3, plot4]
    return plotlist