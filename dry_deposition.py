# -*- coding: utf-8 -*-
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from helper_functions import ds_sel_yr, annual_avg
from diff_plots_Hg import diff_plots

def dry_dep_plots(Dataset_OLD, Dataset_NEW, Year = None):
    """Main script for calling different routines that produce dry deposition plots
    
    Parameters
    ----------
    Dataset_OLD : xarray dataset
        Reference Model dataset (total wet deposition)
                
    Dataset_NEW : xarray dataset
        New Model dataset (total wet deposition)
            
    Year : int or list of int, optional
        Optional parameter to only select subset of years    
    
    """
    # Hg0 dry deposition
    
    # Allow subsetting for years, if inputted into the function
    # temporarily set to 2014 since only have data from this year from my reference run
    OLD_Hg0_ddep_yr = ds_sel_yr(Dataset_OLD, 'DryDep_Hg0', 2014) # AF - must change

    NEW_Hg0_ddep_yr = ds_sel_yr(Dataset_NEW, 'DryDep_Hg0', Year)
    
    
    # calculate annual average
    OLD_Hg0_ddep = annual_avg(OLD_Hg0_ddep_yr)
    NEW_Hg0_ddep = annual_avg(NEW_Hg0_ddep_yr)
    
    # Convert model data from molec/cm^2/s to ug/m^2/yr for annual average   
    s_in_yr = 3.154e7 # seconds in year
    ug_g = 1e6 # ug in g
    cm2_m2 = 1e4 # cm^2 in m^2
    MW_Hg = 200.59 # g mol^-1
    avo = 6.02e23 # avogadro number molec mol^-1
    
    unit_conv = MW_Hg / avo * ug_g * cm2_m2 * s_in_yr # constant to convert units
        
    OLD_Hg0_ddep = OLD_Hg0_ddep * unit_conv # ug/m^2/yr
    NEW_Hg0_ddep = NEW_Hg0_ddep * unit_conv # ug/m^2/yr
        
    # Plot Hg0 dry deposition difference plot 
    plot1 = diff_plots(OLD_Hg0_ddep, NEW_Hg0_ddep, 
                       Units="\u03BCg m$^{-2}$ yr$^{-1}$ ",
                       Title="Hg(0) Dry Dep")

    # HgII + HgP dry deposition
    
    # Allow subsetting for years, if inputted into the function
    # temporarily set to 2014 since only have data from this year from my reference run
    OLD_Hg2_ddep_yr = ds_sel_yr(Dataset_OLD, 'DryDep_Hg2', 2014) # AF - must change
    OLD_HgP_ddep_yr = ds_sel_yr(Dataset_OLD, 'DryDep_HgP', 2014) # AF - must change

    NEW_Hg2_ddep_yr = ds_sel_yr(Dataset_NEW, 'DryDep_Hg2', Year)
    NEW_HgP_ddep_yr = ds_sel_yr(Dataset_NEW, 'DryDep_HgP', Year)
    
    
    # calculate annual average
    OLD_Hg2_ddep = annual_avg(OLD_Hg2_ddep_yr)
    OLD_HgP_ddep = annual_avg(OLD_HgP_ddep_yr)
    NEW_Hg2_ddep = annual_avg(NEW_Hg2_ddep_yr)
    NEW_HgP_ddep = annual_avg(NEW_HgP_ddep_yr)
    
    # Sum gas-phase Hg(II) with particulate Hg(P)
    OLD_Hg2_tot_ddep = OLD_Hg2_ddep + OLD_HgP_ddep
    NEW_Hg2_tot_ddep = NEW_Hg2_ddep + NEW_HgP_ddep
    
    # Converting units to ug/m^2/yr
    OLD_Hg2_tot_ddep = OLD_Hg2_tot_ddep * unit_conv # ug/m^2/yr
    NEW_Hg2_tot_ddep = NEW_Hg2_tot_ddep * unit_conv # ug/m^2/yr
    
    # Plot Hg(II) + Hg(P) dry deposition difference plot 
    plot2 = diff_plots(OLD_Hg2_tot_ddep, NEW_Hg2_tot_ddep, 
                       Units="\u03BCg m$^{-2}$ yr$^{-1}$ ",
                       Title="Hg(II)+Hg(P) Dry Dep")
    
    return plot1, plot2