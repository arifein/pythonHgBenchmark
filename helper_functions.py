import xarray as xr
from calendar import monthrange
import numpy as np
from math import log10, floor
import sys

def open_Hg (fn_OLD, fn_NEW):
    """ Open GEOSChem.* netcdf files for Hg species as an xarray dataset.
    Can take either list of monthly-averaged files (including *), or a 
    single time-concatenated file for each simulation.

    Parameters
    ----------
    fn_OLD : string
        Reference model filename(s). If it includes an asterisk, will open as
        multi-file dataset.
    fn_NEW : string
        New model filename(s). If it includes an asterisk, will open as
        multi-file dataset.
    """
    
    # check if dataset contains '*', that means it is multi-file
    if '*' in fn_OLD: # multiple files to open
        ds_OLD = xr.open_mfdataset(fn_OLD)
    else: # single filename
        ds_OLD = xr.open_dataset(fn_OLD)
            
    if '*' in fn_NEW: # multiple files to open
        ds_NEW = xr.open_mfdataset(fn_NEW)
    else: # single filename
        ds_NEW = xr.open_dataset(fn_NEW)
    
    return ds_OLD, ds_NEW

def ds_sel_yr (ds, varname, Year):
    """ If a year is given, then subset load for that year. Otherwise load all data into variable

    Parameters
    ----------
    ds : xarray dataset
        Dataset of simulation to extract data from
    varname : string
        Name of parameter to extract
    Year : int
        Subset of year(s) to analyze model data  
        
    """
    
    if Year is not None: # take average over subset of years
        var_yr = ds[varname].sel(time=ds.time.dt.year.isin(Year))
    else: # use all years
        var_yr = ds[varname]

    return var_yr

def annual_avg (var_to_avg):
    """ Take annual average from monthly data for multi-year data, accounting for the difference in day number
    Parameters
    ----------
    var_to_avg : xarray dataArray
        variable to average
    """
    time_v = var_to_avg.time
    
    #first check that the data is in monthly time resolution
    diff = time_v.dt.month[1]- time_v.dt.month[0]
    
    if diff > 1: # more than monthly time difference
        print("Simulation needs to be in monthly-averaged time resolution. This data is super-monthly resolution")
        sys.exit(1)
    elif diff<1: # less than monthly time difference
        print("Simulation needs to be in monthly-averaged time resolution. This data is sub-monthly resolution")
        sys.exit(1)
    else:
        # calculate number of days per month
        days_in_month = time_v.dt.days_in_month
        
        wgts = days_in_month.groupby("time.year") / days_in_month.groupby("time.year").sum()
        
        # Calculate the numerator for weighted average
        obs_sum = (var_to_avg * wgts).resample(time="AS").sum(dim="time").squeeze()
        
        # Calculate the denominator for weighted average
        ones_out = (wgts).resample(time="AS").sum(dim="time").squeeze()
        
        # Return the weighted average
        return obs_sum / ones_out
    
def round_sig(x, sig=2):
    # Rounding one number to specific number of significant digits
    return round(x, sig-int(floor(log10(abs(x))))-1)