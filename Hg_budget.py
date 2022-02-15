# -*- coding: utf-8 -*-
import xarray as xr
import numpy as np
import pandas as pd
from helper_functions import ds_sel_yr, annual_avg, round_sig
import matplotlib.pyplot as plt

def budget_calc(Dataset_OLD, Dataset_NEW, Year = None):
    """Main script for calling different routines that calculate budget terms for simulations
    
    Parameters
    ----------
    Dataset_OLD : xarray dataset
        Reference Model dataset (Hg budget)
                
    Dataset_NEW : xarray dataset
        New Model dataset (Hg budget)
            
    Year : int or list of int, optional
        Optional parameter to only select subset of years    
    
    """
    var_names = ["EmisHg0anthro","EmisHg2HgPanthro","EmisHg0geogenic",
                  "EmisHg0soil", "EmisHg0biomass","EmisHg0land","EmisHg0snow",
                  "DryDep_Hg0","DryDep_Hg2","DryDep_HgP",
                  "WetLossTot_Hg2","WetLossTot_HgP",
                  "FluxHg0fromAirToOcean","FluxHg0fromOceanToAir",
                  "LossHg2bySeaSalt","Gross_Hg_Ox","ProdHg2fromHg0"]
    
    # Call script to extract global fluxes based on variable names
    df_budg = glob_vals(Dataset_OLD, Dataset_NEW, var_names, Year)
    
    # Plot fluxes as a pdf
    plot1 = plot_budget(df_budg)
    return plot1

def glob_vals(ds1, ds2, vars_to_use, Year = None):
    """Load globally, annually averaged values of desired variables
    
    Parameters
    ----------
    ds1 : xarray dataset
        Reference Model dataset (Hg budget)
                
    ds2 : xarray dataset
        New Model dataset (Hg budget)
        
    vars_to_use : list of strings
        variables to extract and average from dataset  
            
    Year : int or list of int, optional
        Optional parameter to only select subset of years    
    
    """
    # Create pandas dataframe with columns named after variables
    df = pd.DataFrame(columns=vars_to_use)
    
    # Fill dataframe using list
    for ivar in vars_to_use:
        # Allow subsetting for years, if inputted into the function
        OLD_var_yr = ds_sel_yr(ds1, ivar, Year) 
        NEW_var_yr = ds_sel_yr(ds2, ivar, Year)
        
        # calculate annual average
        OLD_var = annual_avg(OLD_var_yr)
        NEW_var = annual_avg(NEW_var_yr)
        
        # convert units from kg/s to Mg/yr
        s_in_yr = 365.2425 * 24 * 3600 # s in one year
        kg_Mg = 1e-3 # kg in Mg 
    
        unit_conv = s_in_yr * kg_Mg 
        
        OLD_var = OLD_var * unit_conv # Mg/yr
        NEW_var = NEW_var * unit_conv # Mg/yr
        
        # fill dataframe column with data from this variable
        df[ivar] = [OLD_var.values.item(), NEW_var.values.item()]

    # Name rows after simulations
    df.index = ['Ref','New']

    return df

def plot_budget(df):
    """Plot the budget, as a table for now
    
    Parameters
    ----------
    df : pandas DataFrame
        DataFrame with Hg budget terms
                    
    """
    # change table to 5 significant digits, for easier viewing:
    rs = np.vectorize(round_sig) # vectorize function
    sig_data = rs(df.values,5) # data with 5 significant digits
    
    # Plot DataFrame as a table
    # https://stackoverflow.com/questions/32137396/how-do-i-plot-only-a-table-in-matplotlib
    fig, ax =plt.subplots(figsize=(16, 4))
    ax.axis('tight')
    ax.axis('off')
    the_table = ax.table(cellText=sig_data,colLabels=df.columns,
                         rowLabels=df.index, loc='center')
    the_table.set_fontsize(12)
    the_table.scale(1.2, 1.2)  # may help with fontsize

    return fig
