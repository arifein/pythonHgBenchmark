import xarray as xr
def open_Hg (fn_OLD, fn_NEW):
    """ Open GEOSChem.* netcdf files for Hg species as an xarray dataset.
    Can take either list of hourly/monthly/daily files (including *), or a 
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