import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import cartopy.crs as ccrs

# Define a function for comparing new and old maps of a certain variable
def diff_plots (Var_OLD, Var_NEW, Units="ng/m$^3$", Title="Surface TGM"):

    """ Plot the mean values of a certain Hg variable for both the reference and new simulations. 
    Produce the absolute and percent differences between the reference and new simulations.

     
    Args:
    Var_OLD (xarray.DataArray) : Reference Model variable
    Var_NEW (xarray.DataArray) : New Model variable 
    Units (str) : Name of the units the data is in.
    Title (str) : Title of your graph. 
    
    
    """            
    # Find the absolute difference between the reference and new model.
    Abs_diff = Var_NEW - Var_OLD
    # Find the absolute maximum value of the absolute difference. 
    Abs_MaxVal= np.max(np.abs(Abs_diff.values))
    
    # Find the percent difference of the models.  
    Perc_diff = (Abs_diff / Var_OLD)*100
    
    # Find the absolute maximum value of the percent  difference. 
    Perc_diff_non_nan = Perc_diff.values[~np.isnan(Perc_diff.values)] # Only check non_nan values for limit
    Perc_MaxVal= np.max(np.abs(Perc_diff_non_nan)) # for plotting limits
    
    # Set limit to MaxVal as 100%, since can't have negative numbers
    Perc_MaxVal = min(Perc_MaxVal, 100)
    
    #find maxval and minval for plotting (5th and 99th percentile)
    ref_minval = min(np.percentile(Var_OLD,5), np.percentile(Var_NEW,5))
    ref_maxval = max(np.percentile(Var_OLD,99), np.percentile(Var_NEW,99))
    
    # Plot the four graphs as subplots.
    diffGraph,  axes = plt.subplots(2, 2, figsize=[16,12],
                                   subplot_kw=dict(projection=ccrs.PlateCarree()),
                                   gridspec_kw=dict(hspace=0.2, wspace=0.1))
    axes = axes.flatten()
    
    # Plot the reference model and use a geographical map.
    ax = axes[0]
    im = Var_OLD.plot.pcolormesh(x='lon',y='lat',rasterized = True, ax=ax, 
                             vmin=ref_minval, vmax=ref_maxval,
                             transform=ccrs.PlateCarree(), cmap='viridis',                               
                             cbar_kwargs={'orientation':'horizontal',
                                      'label': Units,
                                      'fraction':0.046,
                                      'pad':0.04})
    # Add a title 
    ax.set_title(' Reference Model Version: '+Title)     
    # Add the coastlines.
    ax.coastlines()
    
    
    
    # Plot the new model using a geographical map.       
    ax = axes[1]
    im = Var_NEW.plot.pcolormesh(x='lon',y='lat',rasterized = True, ax=ax,
                              vmin=ref_minval, vmax=ref_maxval,
                              cmap='viridis', transform=ccrs.PlateCarree(),
                              cbar_kwargs={'orientation':'horizontal',
                                      'label': Units,
                                      'fraction':0.046,
                                      'pad':0.04})
    # Add a title.
    ax.set_title('New Model Version: '+ Title)
    # Add the coastlines.
    ax.coastlines()
    
    # Plot the absolute difference using a geograpical map.
    ax = axes[2]
    im = Abs_diff.plot.pcolormesh(x='lon',y='lat', ax=ax,transform=ccrs.PlateCarree(), 
                                  rasterized = True, cmap='RdBu_r',
                             vmin=-Abs_MaxVal, vmax=Abs_MaxVal,
                          cbar_kwargs={'orientation':'horizontal',
                                      'ticklocation':'auto',
                                      'label': Units,
                                      'fraction':0.046,
                                      'pad':0.04})
    # Add a title.
    ax.set_title("Absolute Difference")
    # Add the coastlines.
    ax.coastlines()
     
        
    ax = axes[3]
    im = Perc_diff.plot.pcolormesh(x='lon',y='lat',ax=ax,transform=ccrs.PlateCarree(), 
                                   rasterized = True, cmap='RdBu_r',
                              vmin=(-Perc_MaxVal), vmax=(Perc_MaxVal),
                        cbar_kwargs={'orientation':'horizontal',
                                      'ticklocation':'auto',
                                      'label':"%" ,
                                      'fraction':0.046,
                                      'pad':0.04})
                                        
    # Add a title.
    ax.set_title("Percent Difference (%)")
    # Add the coastlines. 
    ax.coastlines()
     
    # Show the four subplots.
    plt.show()

    # Return the four graphs. 
    return diffGraph

# Define a function for comparing new and old altitude-latitude plots of a certain variable
def diff_profiles (Var_OLD, Var_NEW, Units="ng/m$^3$", Title="Surface TGM"):

    """ Plot the mean zonal profiles of a certain Hg variable for both the reference and new simulations. 
    Produce the absolute and percent differences between the reference and new simulations.

     
    Args:
    Var_OLD (xarray.DataArray) : Reference Model variable
    Var_NEW (xarray.DataArray) : New Model variable 
    Units (str) : Name of the units the data is in.
    Title (str) : Title of your graph. 
    
    
    """            
    # Find the absolute difference between the reference and new model.
    Abs_diff = Var_NEW - Var_OLD
    # Find the absolute maximum value of the absolute difference. 
    Abs_MaxVal= np.max(np.abs(Abs_diff.values))
    
    # Find the percent difference of the models.  
    Perc_diff = (Abs_diff / Var_OLD)*100
    
    # Find the absolute maximum value of the percent  difference. 
    Perc_diff_non_nan = Perc_diff.values[~np.isnan(Perc_diff.values)] # Only check non_nan values for limit
    Perc_MaxVal= np.max(np.abs(Perc_diff_non_nan)) # for plotting limits
    
    # Set limit to MaxVal as 100%, since can't have negative numbers
    Perc_MaxVal = min(Perc_MaxVal, 100)
    
    #find maxval and minval for plotting (5th and 99th percentile)
    ref_minval = min(np.percentile(Var_OLD,5), np.percentile(Var_NEW,5))
    ref_maxval = max(np.percentile(Var_OLD,99), np.percentile(Var_NEW,99))
    
    # For y-axis replace model levels by altitude in a column with surface sea level pressure
    alt_levels = [0.058, 0.189, 0.32, 0.454, 0.589, 0.726, 0.864, 1.004, 1.146,
                  1.29, 1.436, 1.584, 1.759, 1.988, 2.249, 2.517, 2.792, 3.074,
                  3.439, 3.896, 4.375, 4.879, 5.413, 5.98, 6.585, 7.237, 7.943,
                  8.846, 9.936, 11.021, 12.086, 13.134, 14.17, 15.198, 16.222,
                  17.243, 18.727, 20.836, 23.02, 25.307, 28.654, 34.024, 40.166,
                  47.135, 54.834, 63.053, 72.18] 
    
    lat_var = Var_NEW['lat'] # lat, for plotting
    
    # Plot the four graphs as subplots.
    diffGraph_v,  axes = plt.subplots(2, 2, figsize=[16,12],
                                   gridspec_kw=dict(hspace=0.25, wspace=0.25))
    axes = axes.flatten()
    
    # Plot the reference model
    ax = axes[0]
    im = ax.pcolormesh(lat_var,alt_levels, Var_OLD,
                       rasterized = True, 
                       vmin=ref_minval, vmax=ref_maxval, shading='auto')
    # Add a colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(Units)

    # Add a title 
    ax.set_title(' Reference Model Version: '+Title)         
    ax.set_ylim([0, 35]) # restrict y limit
    ax.set_ylabel('Model level ~ altitude (km)')
    ax.set_xlabel('Latitude (°N)')

    # Plot the new model       
    ax = axes[1]
    im = ax.pcolormesh(lat_var, alt_levels, Var_NEW,
                       rasterized = True, 
                       vmin=ref_minval, vmax=ref_maxval, shading='auto')
    # Add a colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(Units)

    # Add a title.
    ax.set_title('New Model Version: '+ Title)
    ax.set_ylim([0, 35]) # restrict y limit
    ax.set_ylabel('Model level ~ altitude (km)')
    ax.set_xlabel('Latitude (°N)')
   
    # Plot the absolute difference using a geograpical map.
    ax = axes[2]
    im = ax.pcolormesh(lat_var, alt_levels, Abs_diff,
                        rasterized = True, cmap='RdBu_r',
                        vmin=(-Abs_MaxVal), vmax=Abs_MaxVal, shading='auto')


    # Add a colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(Units)
    # Add a title.
    ax.set_title("Absolute Difference")
    
    ax.set_ylim([0, 35]) # restrict y limit
    ax.set_ylabel('Model level ~ altitude (km)')
    ax.set_xlabel('Latitude (°N)')

    ax = axes[3]
    im = ax.pcolormesh(lat_var, alt_levels, Perc_diff,
                    rasterized = True, cmap='RdBu_r',
                    vmin=(-Perc_MaxVal), vmax=Perc_MaxVal, shading='auto')
   
    # Add a colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('%')
                                        
    # Add a title.
    ax.set_title("Percent Difference (%)")
    
    ax.set_ylim([0, 35]) # restrict y limit
    ax.set_ylabel('Model level ~ altitude (km)')
    ax.set_xlabel('Latitude (°N)')

    # Show the four subplots.
    plt.show()

    # Return the four graphs. 
    return diffGraph_v

