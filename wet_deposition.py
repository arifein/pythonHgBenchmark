import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from load_Hgmodel_data import ds_sel_yr

def wet_dep_plots(Dataset_OLD_LS, Dataset_OLD_CV, Dataset_NEW_LS, Dataset_NEW_CV, Year = None):
    """Main script for calling different routines that produce wet deposition map plots
    
    Parameters
    ----------
    Dataset_OLD_LS : xarray dataset
        Reference Model dataset (wet deposition from large-scale precipitation)
        
    Dataset_OLD_CV : xarray dataset
        Reference Model dataset (wet deposition from convective precipitation)
        
    Dataset_NEW_LS : xarray dataset
        New Model dataset (wet deposition from large-scale precipitation)
        
    Dataset_NEW_CV : xarray dataset
        New Model dataset (wet deposition from convective precipitation)
    
    Year : int or list of int, optional
        Optional parameter to only select subset of years    
    
    """
    
    return
            
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

