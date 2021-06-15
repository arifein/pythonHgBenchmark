import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from load_Hgmodel_data import ds_sel_yr, annual_avg
import cartopy.crs as ccrs
import cartopy.feature as cf
from matplotlib import colors
from scipy import stats

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

    # Plot MDN comparison maps 
    plot1, plot2 = MDN_USA(OLD_Hg_totwdep, NEW_Hg_totwdep, Year)
    
    return plot1, plot2
            
def MDN_USA(totwetdep_OLD, totwetdep_NEW, Year):
    """Plot the reference and new simulations wet deposition map against observations from the MDN network
    
    Parameters
    ----------
    totwetdep_OLD : xarray DataArray
        Reference Model dataset (wet deposition)
    totwetdep_NEW : xarray DataArray
        New Model dataset (wet deposition)
        
    """
    # Read observations data file
    MDN_Hg = pd.read_csv('data/MDNannual.dat',skiprows=[0], na_values=(-9999))
    MDN_Hg.columns=['SiteID', 'Lat', 'Lon', 'Hg_Dep_ngm2yr']

    # Set variable names for the longitude and latitude in the dataset.
    Lon_MDN = MDN_Hg['Lon']
    Lat_MDN = MDN_Hg['Lat']
    
    # Set a variable name for the observed values of Hg and find the mean of these values. 
    Value_MDN = MDN_Hg['Hg_Dep_ngm2yr']    
    
    # unit conversion from ng/m2/yr to ug/m^2/yr
    Value_MDN = Value_MDN * 1e-3
    # Convert model data from kg/s to ug/m^2/yr    
    # Load grid cell area for unit conversion of model
    fn_gbox = 'data/GEOSChem_2x25_gboxarea.nc'
    ds_gbox = xr.open_dataset(fn_gbox)
    gbox_GC = ds_gbox.cell_area
    
    s_in_yr = 365.2425 * 24 * 3600 # s in one year
    kg_ug = 1e9 # kg in ug 
    
    unit_conv = s_in_yr * kg_ug / gbox_GC
    
    totwetdep_OLD = totwetdep_OLD * unit_conv # ug/m^2/yr
    totwetdep_NEW = totwetdep_NEW * unit_conv # ug/m^2/yr

    # Calculate statistics to add to plot
    
    # Round the mean of the observations to 1 decimal places. 
    MeObsDP=round(np.mean(Value_MDN),1)
    
    # Find the standard deviation of the observed Hg and round this value to 1 decimal places.
    StdObs= round(np.std(Value_MDN),1)

    # Create an array of numpy zeros to fill for the reference and new models based on the amount of observed values.
    OLDval=np.zeros(len(Lat_MDN))
    NEWval=np.zeros(len(Lat_MDN))
    
    # Create a for loop to extract values from the models based on the latitude and longitude of the observations.
    for i in range(len(Lat_MDN)):
        OLDval[i]= totwetdep_OLD.sel(lat=[Lat_MDN[i]], lon=[Lon_MDN[i]], method='nearest')
        NEWval[i]= totwetdep_NEW.sel(lat=[Lat_MDN[i]], lon=[Lon_MDN[i]], method='nearest')
    
    # Find the standard deviation of the values extracted from the reference and new model and round this value 
    # to 2 decimal places.
    ErrOLD= round(np.std(OLDval),1)
    ErrNEW= round(np.std(NEWval),1)

    # Round the means of the values extracted from the reference and new models to decimal places.
    MeMoOl=round(np.mean(OLDval),1)
    MeMoNe= round(np.mean(NEWval),1)

    # Find the correlation coefficient 
    corrOld, _ = stats.pearsonr(Value_MDN, OLDval)
    corrNew, _ = stats.pearsonr(Value_MDN, NEWval)    
    # Find the coefficient of determination for the reference and new models
    CoeffOld= corrOld ** 2
    CoeffNew= corrNew ** 2


    # Create text strings for relevant information: mean, coefficient of determination (rounding to 3DP),
    textstr1= "Mean Obs. = %s $\pm$ %s \u03BCg m$^{-2}$ yr$^{-1}$ "%(MeObsDP, StdObs)
    textstr2= "Mean Mod. = %s $\pm$ %s \u03BCg m$^{-2}$ yr$^{-1}$ "%(MeMoOl ,ErrOLD)
    textstr3= "Mean Mod. = %s $\pm$ %s \u03BCg m$^{-2}$ yr$^{-1}$ "%(MeMoNe ,ErrNEW)
    textstr4= "$R^2$= %s" %(round(CoeffOld,3))
    textstr5= "$R^2$= %s" %(round(CoeffNew,3))

    
    # Add a figure.
    OLDMAP = plt.figure(figsize=[10,4.5])
    
    # Add a geographical projection on the map.
    ax = OLDMAP.add_subplot(111, projection=ccrs.PlateCarree())
    
    # Set color limits as 0, 20
    norm = colors.Normalize(vmin=0.,vmax=20.0)

    # Plot the reference model on the projection.
    im = totwetdep_OLD.plot.pcolormesh( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                                      norm=norm, cmap='viridis',rasterized = True, zorder=0,
                                 cbar_kwargs={'orientation':'horizontal',
                                           'ticklocation':'auto',
                                      'label':"\u03BCg m$^{-2}$ yr$^{-1}$ ",
                                      'extend': 'max',
                                      'fraction':0.046,
                                      'pad':0.04})
    # Show the coastlines.
    ax.coastlines()
    # Add in state borders
    # Create a feature for States/Admin 1 regions at 1:1km from Natural Earth
    states_provinces = cf.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
    
    ax.add_feature(states_provinces, edgecolor='black', linewidth=0.5,zorder=5)
    ax.add_feature(cf.BORDERS,zorder=10)
     
    # Add the observed values to the plot.
    plt.scatter(Lon_MDN, Lat_MDN,  transform=ccrs.PlateCarree(),marker='o',
                norm=norm,
                linewidths=0.5, edgecolors='black', 
                label=None, c=Value_MDN, cmap='viridis', zorder=15)

    # Add a title.
    plt.title('Hg Wet Deposition, Reference Model (' + str(Year) + '), MDN (2015)' ,fontsize=15)   
  
    # Zoom into contiguous USA area
    ax.set_xlim([-128, -64])
    ax.set_ylim([22, 50])
    
    
    # Add text to the plot.
    plt.text(1.05, 0.1, textstr1, fontsize=13, transform=ax.transAxes)
    plt.text(1.05, 0.25, textstr2, fontsize=13, transform=ax.transAxes)
    plt.text(1.05, 0.4, textstr4, fontsize=13, transform=ax.transAxes)


    # Adjust to give text space
    OLDMAP.subplots_adjust(right=0.6, left=0.05)
    
    # Show the plot.
    OLDMAP.show()

    # Add figure for new simulation
    NEWMAP = plt.figure(figsize=[10,4.5])
    
    # Add a geographical projection on the map.
    ax = NEWMAP.add_subplot(111, projection=ccrs.PlateCarree())
    
    # Set color limits as 0, 20
    norm = colors.Normalize(vmin=0.,vmax=20.0)

    # Plot the reference model on the projection.
    im = totwetdep_NEW.plot.pcolormesh( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                                      norm=norm, cmap='viridis',rasterized = True, zorder=0,
                                 cbar_kwargs={'orientation':'horizontal',
                                           'ticklocation':'auto',
                                      'label':"\u03BCg m$^{-2}$ yr$^{-1}$ ",
                                      'extend': 'max',
                                      'fraction':0.046,
                                      'pad':0.04})
    # Show the coastlines.
    ax.coastlines()
    # Add in state borders
    # Create a feature for States/Admin 1 regions at 1:1km from Natural Earth
    states_provinces = cf.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
    
    ax.add_feature(states_provinces, edgecolor='black', linewidth=0.5,zorder=5)
    ax.add_feature(cf.BORDERS,zorder=10)
     
    # Add the observed values to the plot.
    plt.scatter(Lon_MDN, Lat_MDN,  transform=ccrs.PlateCarree(),marker='o',
                norm=norm,
                linewidths=0.5, edgecolors='black', 
                label=None, c=Value_MDN, cmap='viridis', zorder=15)

    # Add a title.
    plt.title('Hg Wet Deposition, New Model (' + str(Year) + '), MDN (2015)' ,fontsize=15)   
  
    # Zoom into contiguous USA area
    ax.set_xlim([-128, -64])
    ax.set_ylim([22, 50])
    
    
    # Add text to the plot.
    plt.text(1.05, 0.1, textstr1, fontsize=13, transform=ax.transAxes)
    plt.text(1.05, 0.25, textstr3, fontsize=13, transform=ax.transAxes)
    plt.text(1.05, 0.4, textstr5, fontsize=13, transform=ax.transAxes)


    # Adjust to give text space
    NEWMAP.subplots_adjust(right=0.6, left=0.05)
    
    # Show the plot.
    NEWMAP.show()
   
    return OLDMAP, NEWMAP