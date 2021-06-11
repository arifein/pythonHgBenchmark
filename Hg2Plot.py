import xarray as xr
import numpy as np
from diff_plots_Hg import diff_plots

def SurfaceHg2(Old_Dataset, New_Dataset, Year = None):
    """ Plot the mean surface Hg2 + HgP for the reference and new models.
    
    Parameters
    ----------
    Dataset_OLD : string
        Reference Model xarray dataset
    Dataset_NEW : string
        New Model xarray dataset 
    
    Year : int or list of int, optional
        Optional parameter to only select subset of years
    
    """    
    
    # Make a variable for the unit conversion factor from vmr to  pg/m^3
    # Now more traceable
    R = 8.314462 # m^3 Pa K^-1 mol ^-1
    MW_Hg = 200.59 # g mol^-1
    pg_g = 1e12 # pg/g
    
    stdpressure = 101325 # Pascals
    stdtemp = 273.15 # Kelvins
    
    unit_conv = stdpressure / R / stdtemp * MW_Hg * pg_g # converter from vmr to pg m^-3
    
    # Allow subsetting for years, if inputted into the function
    if Year is not None: # take average over subset of years
        # OLD simulation        
        OLD_HgP_yr = Old_Dataset.SpeciesConc_HgP.sel(time=Old_Dataset.time.dt.year.isin(Year))
        OLD_Hg2_yr = Old_Dataset.SpeciesConc_Hg2.sel(time=Old_Dataset.time.dt.year.isin(Year))
        # NEW simulation        
        NEW_HgP_yr = New_Dataset.SpeciesConc_HgP.sel(time=New_Dataset.time.dt.year.isin(Year))
        NEW_Hg2_yr = New_Dataset.SpeciesConc_Hg2.sel(time=New_Dataset.time.dt.year.isin(Year))
    else: # use all years
        # OLD simulation        
        OLD_HgP_yr = Old_Dataset.SpeciesConc_HgP
        OLD_Hg2_yr = Old_Dataset.SpeciesConc_Hg2
        # NEW simulation                
        NEW_HgP_yr = New_Dataset.SpeciesConc_HgP
        NEW_Hg2_yr = New_Dataset.SpeciesConc_Hg2
   
    # Extract and add together Hg2 and HgP at the surface from both 
    # model simulations, multiplying by the unit conversion factor
    OLD_HgP = OLD_HgP_yr.isel(lev=0).mean('time')              
    OLD_Hg2 = OLD_Hg2_yr.isel(lev=0).mean('time')

    NEW_HgP = NEW_HgP_yr.isel(lev=0).mean('time')                
    NEW_Hg2 = NEW_Hg2_yr.isel(lev=0).mean('time')
                       
    Hg2_tot_Old = (OLD_HgP + OLD_Hg2) * unit_conv # sum of HgP and Hg2
    Hg2_tot_New = (NEW_HgP + NEW_Hg2) * unit_conv # sum of HgP and Hg2
        
    # Create difference plot for Hg2 + HgP at surface 
    diff_MAP = diff_plots(Hg2_tot_Old, Hg2_tot_New, Units="pg/m$^3$",
                          Title="Surface Hg(II)+Hg(P)")
    
    # return plot, for saving in PDF 
    return diff_MAP
