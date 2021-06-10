import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import cartopy.crs as ccrs
#%matplotlib inline

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
    Perc_MaxVal= np.max(np.abs(Perc_diff.values))
    
    #find maxval and minval for plotting (5th and 99th percentile)
    ref_minval = round_sig(np.percentile(Var_OLD,5))
    ref_maxval = round_sig(np.percentile(Var_OLD,99))
    
    # Plot the four graphs as subplots.
    TGMGraph = plt.figure(figsize=(15,10))
    
    # Plot the reference model and use a geographical map.
    ax = plt.subplot(221, projection=ccrs.PlateCarree())
    im=Var_OLD.plot.contourf(x='lon',y='lat',levels=11, ax=ax, 
                             vmin=ref_minval, vmax=ref_maxval,
                             transform=ccrs.PlateCarree(), cmap='viridis', 
                              
                             cbar_kwargs={'orientation':'horizontal',
                                      'label': Units})  
    # Add a title 
    plt.title(' Reference Model Version: '+Title)     
    # Add the coastlines.
    ax.coastlines()
    
    
    
    # Plot the new model using a geographical map.       
    ax = plt.subplot(222, projection=ccrs.PlateCarree())
    im= Var_NEW.plot.contourf(x='lon',y='lat',levels=11, ax=ax,
                              vmin=ref_minval, vmax=ref_maxval,
                              cmap='viridis', transform=ccrs.PlateCarree(),
                              cbar_kwargs={'orientation':'horizontal',
                                      'label': Units})
    # Add a title.
    plt.title('New Model Version: '+ Title)
    # Add the coastlines.
    ax.coastlines()
    
    # Plot the absolute difference using a geograpical map.
    ax = plt.subplot(223, projection=ccrs.PlateCarree())
    im= Abs_diff.plot.imshow(x='lon',y='lat', ax=ax,transform=ccrs.PlateCarree(),  cmap='RdBu_r',
                             vmin=-Abs_MaxVal, vmax=Abs_MaxVal,
                          cbar_kwargs={'orientation':'horizontal',
                                      'ticklocation':'auto',
                                      'label': Units})
    # Add a title.
    plt.title("Absolute Difference")
    # Add the coastlines.
    ax.coastlines()
     
    # Plot the percent difference using a geographical map.
    ax = plt.subplot(224, projection=ccrs.PlateCarree())
    im= Perc_diff.plot.imshow(x='lon',y='lat',ax=ax,transform=ccrs.PlateCarree(), cmap='RdBu_r',
                              vmin=(-Perc_MaxVal), vmax=(Perc_MaxVal),
                        cbar_kwargs={'orientation':'horizontal',
                                      'ticklocation':'auto',
                                      'label':"%" })
    # Add a title.
    plt.title("Percent Difference (%)")
    # Add the coastlines. 
    ax.coastlines()
 
    plt.tight_layout()
    
    # Show the four subplots.
    plt.show()

    # Return the four graphs. 
    return TGMGraph
    
def round_sig(x, sig=2):  # round to 2 significant figures
   return round(x, sig-int(np.floor(np.log10(abs(x))))-1)

