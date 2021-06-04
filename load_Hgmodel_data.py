import xarray as xr
def open_Hg_spc (fn_OLD, fn_NEW):
    """ Open GEOSChem.SpeciesConc.* netcdf files for Hg species as an xarray dataset.
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

