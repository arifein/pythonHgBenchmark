import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from load_Hgmodel_data import ds_sel_yr, annual_avg

def wet_dep_plots(Dataset_OLD, Dataset_NEW, Year = None):
    """Main script for calling different routines that produce wet deposition map plots
    
    Parameters
    ----------
    Dataset_OLD : xarray dataset
        Reference Model dataset (total wet deposition)
                
    Dataset_NEW : xarray dataset
        New Model dataset (total wet deposition)
            
    Year : int or list of int, optional
        Optional parameter to only select subset of years    
    
    """
    # Allow subsetting for years, if inputted into the function
    
    # temporarily set to 2014 since only have data from this year from my reference run
    OLD_Hg_totwdep_yr = ds_sel_yr(Dataset_OLD, 'WetLossTot_Hg', 2014) # AF - must change

    NEW_Hg_totwdep_yr = ds_sel_yr(Dataset_NEW, 'WetLossTot_Hg', Year)
    
    
    # calculate annual average
    OLD_Hg_totwdep = annual_avg(OLD_Hg_totwdep_yr)
    NEW_Hg_totwdep = annual_avg(NEW_Hg_totwdep_yr)
    return OLD_Hg_totwdep, NEW_Hg_totwdep

    #for total wet deposition, sum variables    
    OLD_Hg_totwdep_yr = OLD_HgP_ls_yr + OLD_Hg2_ls_yr + OLD_HgP_cv_yr + OLD_Hg2_cv_yr
    NEW_Hg_totwdep_yr = NEW_HgP_ls_yr + NEW_Hg2_ls_yr + NEW_HgP_cv_yr + NEW_Hg2_cv_yr
        
    # sum all data over levels and calculate annual average
    OLD_Hg_totwdep = annual_avg(OLD_Hg_totwdep_yr.sum('lev'))
    NEW_Hg_totwdep = annual_avg(NEW_Hg_totwdep_yr.sum('lev'))
    
    # OLD_HgP_ls_t = annual_avg(OLD_HgP_ls_yr.sum('lev'))
    # OLD_Hg2_ls_t = annual_avg(OLD_Hg2_ls_yr.sum('lev'))
    # OLD_HgP_cv_t = annual_avg(OLD_HgP_cv_yr.sum('lev'))
    # OLD_Hg2_cv_t = annual_avg(OLD_Hg2_cv_yr.sum('lev'))

    # NEW_HgP_ls_t = annual_avg(NEW_HgP_ls_yr.sum('lev'))
    # NEW_Hg2_ls_t = annual_avg(NEW_Hg2_ls_yr.sum('lev'))
    # NEW_HgP_cv_t = annual_avg(NEW_HgP_cv_yr.sum('lev'))
    # NEW_Hg2_cv_t = annual_avg(NEW_Hg2_cv_yr.sum('lev'))
    
    # for total wet deposition, sum variables    
    # OLD_Hg_totwdep = OLD_HgP_ls_t + OLD_Hg2_ls_t + OLD_HgP_cv_t + OLD_Hg2_cv_t
    # NEW_Hg_totwdep = NEW_HgP_ls_t + NEW_Hg2_ls_t + NEW_HgP_cv_t + NEW_Hg2_cv_t
    
    return OLD_Hg_totwdep, NEW_Hg_totwdep
            
def MDN_USA(Dataset_OLD, Dataset_NEW, Year = None):
    """Plot the reference and new simulations wet deposition map against observations from the MDN network
    
    Parameters
    ----------
    Dataset_OLD : xarray dataset
        Reference Model dataset (wet deposition)
    Dataset_NEW : xarray dataset
        New Model dataset (wet deposition)
    
    Year : int or list of int, optional
        Optional parameter to only select subset of years    
    
    """

