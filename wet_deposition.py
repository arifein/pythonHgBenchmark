# -*- coding: utf-8 -*-
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from helper_functions import ds_sel_yr, annual_avg
import cartopy.crs as ccrs
import cartopy.feature as cf
from matplotlib import colors
from scipy import stats
from diff_plots_Hg import diff_plots

def wet_dep_plots(Dataset_OLD, Dataset_NEW, Year1 = None, Year2 = None):
    """Main script for calling different routines that produce wet deposition map plots
    
    Parameters
    ----------
    Dataset_OLD : xarray dataset
        Reference Model dataset (total wet deposition)
                
    Dataset_NEW : xarray dataset
        New Model dataset (total wet deposition)
            
    Year1 : int or list of int, optional
        Optional parameter to only select subset of years for old sim
    Year2 : int or list of int, optional
        Optional parameter to only select subset of years for new sim
    
    """
    # Allow subsetting for years, if inputted into the function
    
    OLD_Hg_totwdep_yr = ds_sel_yr(Dataset_OLD, 'WetLossTot_Hg', Year1)

    NEW_Hg_totwdep_yr = ds_sel_yr(Dataset_NEW, 'WetLossTot_Hg', Year2)
    
    
    # calculate annual average
    OLD_Hg_totwdep = annual_avg(OLD_Hg_totwdep_yr)
    NEW_Hg_totwdep = annual_avg(NEW_Hg_totwdep_yr)
    
    # Convert model data from kg/s to ug/m^2/yr for annual average   
    # Load grid cell area for unit conversion of model
    fn_gbox = 'data/GEOSChem_2x25_gboxarea.nc'
    ds_gbox = xr.open_dataset(fn_gbox)
    gbox_GC = ds_gbox.cell_area
    
    s_in_yr = 365.2425 * 24 * 3600 # s in one year
    kg_ug = 1e9 # kg in ug 
    
    unit_conv = s_in_yr * kg_ug / gbox_GC
    
    OLD_Hg_totwdep = OLD_Hg_totwdep * unit_conv # ug/m^2/yr
    NEW_Hg_totwdep = NEW_Hg_totwdep * unit_conv # ug/m^2/yr

    # Plot global maps 
    plot1, plot2 = wdep_global(OLD_Hg_totwdep, NEW_Hg_totwdep, Year1, Year2)
    
    # Plot MDN comparison maps 
    plot3, plot4 = MDN_USA(OLD_Hg_totwdep, NEW_Hg_totwdep, Year1, Year2)
    
    # calculate climatology, if have multi-year model runs
    OLD_Hg_totwdep_clim = OLD_Hg_totwdep_yr.groupby('time.month').mean() # kg/s
    NEW_Hg_totwdep_clim = NEW_Hg_totwdep_yr.groupby('time.month').mean() # kg/s

    # Plot MDN seasonal comparison by latitude 
    plot5 = MDN_USA_Seasonal(OLD_Hg_totwdep_clim, NEW_Hg_totwdep_clim, Year1, Year2)

    # Plot EMEP comparison maps 
    plot6, plot7 = EMEP_map(OLD_Hg_totwdep, NEW_Hg_totwdep, Year1, Year2)
    
    # Plot wet deposition difference plot 
    plot8 = diff_plots(OLD_Hg_totwdep, NEW_Hg_totwdep, 
                       Units="\u03BCg m$^{-2}$ yr$^{-1}$ ",
                       Title="Total Wet Dep")
    
    plotlist = [plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8]

    return plotlist

def wdep_global(totwetdep_OLD, totwetdep_NEW, Year1, Year2):
    """Plot the reference and new simulations wet deposition map against global observations, compiled in Shah et al. (2021)
    
    Parameters
    ----------
    totwetdep_OLD : xarray DataArray
        Reference Model dataset (wet deposition)
    totwetdep_NEW : xarray DataArray
        New Model dataset (wet deposition)
        
    Year1 : int or list of int, optional
        Optional parameter to only select subset of years for old sim
    Year2 : int or list of int, optional
        Optional parameter to only select subset of years for new sim
        
    """
    #---Read observations data file for NADP (== MDN)--- 
    NADP_Hg = pd.read_csv('data/NADP_wetdep.csv', na_values=(-9999))
    NADP_Hg.columns=['SiteID', 'Lat', 'Lon', 'Hg_Dep']

    # Set variable names for the longitude and latitude in the dataset.
    Lon_NADP = NADP_Hg['Lon']
    Lat_NADP = NADP_Hg['Lat']
    
    # Set a variable name for the observed values of Hg 
    Value_NADP = NADP_Hg['Hg_Dep']    
    
    # unit conversion from ng/m2/yr to ug/m^2/yr
    Value_NADP = Value_NADP * 1e-3
    
    #---Read observations data file for EMEP--- 
    EMEP_Hg = pd.read_csv('data/EMEP_wetdep.csv', na_values=(-9999))
    EMEP_Hg.columns=['SiteID', 'Lat', 'Lon', 'Hg_Dep']

    # Set variable names for the longitude and latitude in the dataset.
    Lon_EMEP = EMEP_Hg['Lon']
    Lat_EMEP = EMEP_Hg['Lat']
    
    # Set a variable name for the observed values of Hg 
    Value_EMEP = EMEP_Hg['Hg_Dep']    
    
    # unit conversion from ng/m2/yr to ug/m^2/yr
    Value_EMEP = Value_EMEP * 1e-3
    
    #---Read observations data file for GMOS--- 
    GMOS_Hg = pd.read_csv('data/GMOS_wetdep.csv', na_values=(-9999))
    GMOS_Hg.columns=['SiteID', 'Lat', 'Lon', 'Elevation', 'Hg_Dep']

    # Set variable names for the longitude and latitude in the dataset.
    Lon_GMOS = GMOS_Hg['Lon']
    Lat_GMOS = GMOS_Hg['Lat']
    
    # Set a variable name for the observed values of Hg 
    Value_GMOS = GMOS_Hg['Hg_Dep']   
    
    #---Read observations data file for Asia--- 
    Asia_Hg = pd.read_csv('data/Asia_wetdep.csv', na_values=(-9999))
    Asia_Hg.columns=['SiteID', 'Lat', 'Lon', 'Elevation','Hg_Dep',	'Ref', 'HighElev']

    # Set variable names for the longitude and latitude in the dataset.
    Lon_Asia = Asia_Hg['Lon']
    Lat_Asia = Asia_Hg['Lat']
    
    # Set a variable name for the observed values of Hg 
    Value_Asia = Asia_Hg['Hg_Dep']    
    
    # Concatenate all data files
    Value_all = pd.concat([Value_NADP, Value_EMEP, Value_GMOS, Value_Asia], ignore_index=True)
    Lon_all = pd.concat([Lon_NADP, Lon_EMEP, Lon_GMOS, Lon_Asia], ignore_index=True)
    Lat_all = pd.concat([Lat_NADP, Lat_EMEP, Lat_GMOS, Lat_Asia], ignore_index=True)
    
    # Calculate statistics to add to plot
    
    # Round the mean of the observations to 1 decimal places. 
    MeObsDP=round(np.mean(Value_all),1)
    
    # Find the standard deviation of the observed Hg and round this value to 1 decimal places.
    StdObs= round(np.std(Value_all),1)

    # Create an array of numpy zeros to fill for the reference and new models based on the amount of observed values.
    OLDval=np.zeros(len(Lat_all))
    NEWval=np.zeros(len(Lat_all))
    
    # Create a for loop to extract values from the models based on the latitude and longitude of the observations.
    for i in range(len(Lat_all)):
        OLDval[i]= totwetdep_OLD.sel(lat=[Lat_all[i]], lon=[Lon_all[i]], method='nearest')
        NEWval[i]= totwetdep_NEW.sel(lat=[Lat_all[i]], lon=[Lon_all[i]], method='nearest')
    
    # Find the standard deviation of the values extracted from the reference and new model and round this value 
    # to 2 decimal places.
    ErrOLD= round(np.std(OLDval),1)
    ErrNEW= round(np.std(NEWval),1)

    # Round the means of the values extracted from the reference and new models to decimal places.
    MeMoOl=round(np.mean(OLDval),1)
    MeMoNe= round(np.mean(NEWval),1)
    
    print('OBS model median wdep: ' + str(round(np.median(Value_all),1)))
    print('OLD model median wdep: ' + str(round(np.median(OLDval),1)))
    print('NEW model median wdep: ' + str(round(np.median(NEWval),1)))
    
    # Find the correlation coefficient 
    corrOld, _ = stats.pearsonr(Value_all, OLDval)
    corrNew, _ = stats.pearsonr(Value_all, NEWval)    
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
    totwetdep_OLD.plot.pcolormesh(x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                                  norm=norm, cmap='viridis',rasterized = True, zorder=0,
                                  cbar_kwargs={'orientation':'horizontal',
                                               'ticklocation':'auto',
                                               'label':"\u03BCg m$^{-2}$ yr$^{-1}$ ",
                                               'extend': 'max',
                                               'fraction':0.046,
                                               'pad':0.04})
    # Show the coastlines.
    ax.coastlines()
     
    # Add the observed values to the plot.
    plt.scatter(Lon_all, Lat_all,  transform=ccrs.PlateCarree(),marker='o',
                norm=norm,
                linewidths=0.5, edgecolors='black', 
                label=None, c=Value_all, cmap='viridis', zorder=15)

    # Add a title.
    plt.title('Wet Deposition, Reference Model (' + str(Year1) + '), obs from Shah21' ,fontsize=15)       
    
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
    totwetdep_NEW.plot.pcolormesh( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                                  norm=norm, cmap='viridis',rasterized = True, zorder=0,
                                  cbar_kwargs={'orientation':'horizontal',
                                               'ticklocation':'auto',
                                               'label':"\u03BCg m$^{-2}$ yr$^{-1}$ ",
                                               'extend': 'max',
                                               'fraction':0.046,
                                               'pad':0.04})
    # Show the coastlines.
    ax.coastlines()
     
    # Add the observed values to the plot.
    plt.scatter(Lon_all, Lat_all,  transform=ccrs.PlateCarree(),marker='o',
                norm=norm,
                linewidths=0.5, edgecolors='black', 
                label=None, c=Value_all, cmap='viridis', zorder=15)

    # Add a title.
    plt.title('Wet Deposition, New Model (' + str(Year2) + '),  obs from Shah21' ,fontsize=15)   
      
    
    # Add text to the plot.
    plt.text(1.05, 0.1, textstr1, fontsize=13, transform=ax.transAxes)
    plt.text(1.05, 0.25, textstr3, fontsize=13, transform=ax.transAxes)
    plt.text(1.05, 0.4, textstr5, fontsize=13, transform=ax.transAxes)


    # Adjust to give text space
    NEWMAP.subplots_adjust(right=0.6, left=0.05)
    
    # Show the plot.
    NEWMAP.show()
   
    return OLDMAP, NEWMAP


def MDN_USA(totwetdep_OLD, totwetdep_NEW, Year1, Year2):
    """Plot the reference and new simulations wet deposition map against observations from the MDN network
    
    Parameters
    ----------
    totwetdep_OLD : xarray DataArray
        Reference Model dataset (wet deposition)
    totwetdep_NEW : xarray DataArray
        New Model dataset (wet deposition)

    Year1 : int or list of int, optional
        Optional parameter to only select subset of years for old sim
    Year2 : int or list of int, optional
        Optional parameter to only select subset of years for new sim
        
    """
    # Read observations data file
    MDN_Hg = pd.read_csv('data/MDNannual.dat',skiprows=[0], na_values=(-9999))
    MDN_Hg.columns=['SiteID', 'Lat', 'Lon', 'Hg_Dep_ngm2yr']

    # Set variable names for the longitude and latitude in the dataset.
    Lon_MDN = MDN_Hg['Lon']
    Lat_MDN = MDN_Hg['Lat']
    
    # Set a variable name for the observed values of Hg 
    Value_MDN = MDN_Hg['Hg_Dep_ngm2yr']    
    
    # unit conversion from ng/m2/yr to ug/m^2/yr
    Value_MDN = Value_MDN * 1e-3

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
    totwetdep_OLD.plot.pcolormesh(x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
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
    plt.title('Hg Wet Deposition, Reference Model (' + str(Year1) + '), MDN (2015)' ,fontsize=15)   
  
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
    totwetdep_NEW.plot.pcolormesh( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
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
    plt.title('Hg Wet Deposition, New Model (' + str(Year2) + '), MDN (2015)' ,fontsize=15)   
  
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

def MDN_USA_Seasonal(totwetdep_OLD, totwetdep_NEW, Year1, Year2):
    
    """Plot the reference and new simulations wet deposition map against monthly observations from the MDN network
    
    Parameters
    ----------
    totwetdep_OLD : xarray DataArray
        Reference Model dataset (wet deposition)
    totwetdep_NEW : xarray DataArray
        New Model dataset (wet deposition)
        
    Year1 : int or list of int, optional
        Optional parameter to only select subset of years for old sim
    Year2 : int or list of int, optional
        Optional parameter to only select subset of years for new sim
        
    """
    
    # Convert units from kg/s to ng/m2/d
    # Load grid cell area for unit conversion of model
    fn_gbox = 'data/GEOSChem_2x25_gboxarea.nc'
    ds_gbox = xr.open_dataset(fn_gbox)
    gbox_GC = ds_gbox.cell_area
    
    s_in_d = 24 * 3600 # s in one day
    kg_ng = 1e12 # kg in ng 
    
    unit_conv = s_in_d * kg_ng / gbox_GC
    
    totwetdep_OLD = totwetdep_OLD * unit_conv # ng/m^2/d
    totwetdep_NEW = totwetdep_NEW * unit_conv # ng/m^2/d
    
    # Read observations data file
    MDN_Hg_m = pd.read_csv('data/MDNmonthly.dat',skiprows=[0], na_values=(-9999))
    MDN_Hg_m.columns=['SiteID', 'Lat', 'Lon', 'Month','Year', 'Hg_Dep_ngm2d']
    
    # Only look at Eastern USA
    min_lon = -92.5
    max_lon = -72.5
    
    # midpoint latitudes to average
    mid_lat = [30, 34, 38, 42, 46]
    
    # initialize mean and standard deviation arrays
    MDN_Hg_clim = np.zeros((len(mid_lat), 12)) # mean
    MDN_Hg_clim_sd = np.zeros((len(mid_lat), 12)) # standard deviation
    n_sites = np.zeros(len(mid_lat), int) # number of sites
    
    # take average at same model sites - changed from IDL benchmark
    OLD_Hg_clim = np.zeros((len(mid_lat), 12)) # mean
    NEW_Hg_clim = np.zeros((len(mid_lat), 12)) # mean
    
    # loop over latitudes, average sites with these criteria
    for ii, ilat in enumerate(mid_lat):
        # Check whether within lat/lon criteria
        lat_bool = abs(MDN_Hg_m['Lat'] - ilat) <=2 # within two lat units
        lon_bool = (MDN_Hg_m['Lon'] >= min_lon) & (MDN_Hg_m['Lon'] <= max_lon) #within Eastern USA
        
        # Filter dataFrame for these criteria
        MDN_Hg_m_i = MDN_Hg_m.loc[lat_bool & lon_bool, :]
        
        # Find lat and lon of each site
        lon_MDN_i = MDN_Hg_m_i['Lon'][::12].to_numpy()
        lat_MDN_i = MDN_Hg_m_i['Lat'][::12].to_numpy()
        
        # Calculate number of sites
        n_sites[ii] = int(len(MDN_Hg_m_i)/12)
        
        # Calculate monthly mean over lat range
        MDN_Hg_clim[ii,:] = MDN_Hg_m_i.groupby(['Month'])['Hg_Dep_ngm2d'].mean()
        MDN_Hg_clim_sd[ii,:] = MDN_Hg_m_i.groupby(['Month'])['Hg_Dep_ngm2d'].std()    
        
        # Calculate model values at sites
        
        OLD_sites_Hg = np.zeros((n_sites[ii], 12)) # initialize before loop
        NEW_sites_Hg = np.zeros((n_sites[ii], 12)) # initialize before loop

        # Fill in monthly values at grid cells corresponding to MDN sites
        for jj in range(n_sites[ii]):            
            OLD_sites_Hg[jj,:] = totwetdep_OLD.sel(lat=lat_MDN_i[jj], lon=lon_MDN_i[jj], method='nearest')
            NEW_sites_Hg[jj,:] = totwetdep_NEW.sel(lat=lat_MDN_i[jj], lon=lon_MDN_i[jj], method='nearest')
        
        # take modelled average over sites close to the midpoint latitude
        OLD_Hg_clim[ii, :] = np.mean(OLD_sites_Hg, axis=0)
        NEW_Hg_clim[ii, :] = np.mean(NEW_sites_Hg, axis=0)
        
    
    # Plot the zonal climatologies from MDN and simulations (at MDN sites)
    plot = plt.figure(figsize=(8, 9))
    # First add map figure in the background
    ax = plot.add_axes([0.05, 0, 1, .95], projection=ccrs.PlateCarree())
    ax.axis('off') # remove axis border
    
    # Add in state borders
    # Create a feature for States/Admin 1 regions at 1:1km from Natural Earth
    states_provinces = cf.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='110m',
        facecolor='none')
    
    # Show the coastlines, borders, state borders, and lakes
 #   ax.coastlines(edgecolor='gray', linewidth=0.1) 
    ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.2, zorder=5)
    ax.add_feature(cf.BORDERS,zorder=10, edgecolor='gray', linewidth=0.2)
    ax.add_feature(cf.LAKES,zorder=11, edgecolor='gray', linewidth=0.2,facecolor='none')
    ax.add_feature(cf.COASTLINE, edgecolor='gray', linewidth=0.2, zorder=6)
    
    # Zoom into Eastern USA area
    ax.set_xlim([-92.5, -72.5])
    ax.set_ylim([28, 48])
    
    plt.xticks([])
    plt.yticks([])
    
    # Plot inset axes as a loop
    y_pos = np.linspace(0.05, 0.79, num=len(mid_lat))
    
    for ii_r in reversed(range(len(mid_lat))): # reverse index so go from north to south
        iax = plt.axes([0.08, y_pos[ii_r], 0.9, .13], facecolor='none')
    
        # Add the data from the observations, the reference model and the new model
        iax.errorbar(range(12), MDN_Hg_clim[ii_r,:], MDN_Hg_clim_sd[ii_r,:],
                     color= "k", label='MDN 2015')
        iax.plot(range(12),OLD_Hg_clim[ii_r, :],color='blue', label='Reference Model ' + str(Year1))
        iax.plot(range(12),NEW_Hg_clim[ii_r, :],color='red', label='New Model' + str(Year2))
    
    
        # Label the axes, add a legend and add a title
        iax.set_ylabel('ng m$^{-2}$ d$^{-1}$')
        iax.set_title('{0} °N ({1} sites)'.format(mid_lat[ii_r], n_sites[ii_r]), fontsize=12)
    
    
        # Set ticks to every month 
        iax.set_xticks(range(12))
       
        # Set tick labels to month names
        mn = ['J','F','M','A','M','J','J','A','S','O','N','D']            
        iax.set_xticklabels(mn)
        # Add legend to first plot
        if ii_r==len(mid_lat)-1:
            iax.legend(loc='upper left')
            
    # Add overall title to plot
    plot.suptitle('Wet deposition fluxes, Eastern USA',fontweight='bold')
    
    return plot

    
def EMEP_map(totwetdep_OLD, totwetdep_NEW, Year1, Year2):
    """Plot the reference and new simulations wet deposition map against observations from the MDN network
    
    Parameters
    ----------
    totwetdep_OLD : xarray DataArray
        Reference Model dataset (wet deposition)
    totwetdep_NEW : xarray DataArray
        New Model dataset (wet deposition)
        
    Year1 : int or list of int, optional
        Optional parameter to only select subset of years for old sim
    Year2 : int or list of int, optional
        Optional parameter to only select subset of years for new sim
        
    """
    # Read observations data file
    EMEP_Hg = pd.read_csv('data/emep_ann_hg_dep_2013-15.dat',skiprows=[0])
    EMEP_Hg.columns=['SiteID', 'Lat', 'Lon', 'Alt', 'Hg_Dep_ngm2yr']

    # Set variable names for the longitude and latitude in the dataset.
    Lon_EMEP = EMEP_Hg['Lon']
    Lat_EMEP = EMEP_Hg['Lat']
    
    # Set a variable name for the observed values of Hg 
    Value_EMEP = EMEP_Hg['Hg_Dep_ngm2yr']    
    
    # unit conversion from ng/m2/yr to ug/m^2/yr
    Value_EMEP = Value_EMEP * 1e-3
    
    # Calculate statistics to add to plot
    
    # Round the mean of the observations to 1 decimal places. 
    MeObsDP=round(np.mean(Value_EMEP),1)
    
    # Find the standard deviation of the observed Hg and round this value to 1 decimal places.
    StdObs= round(np.std(Value_EMEP),1)

    # Create an array of numpy zeros to fill for the reference and new models based on the amount of observed values.
    OLDval=np.zeros(len(Lat_EMEP))
    NEWval=np.zeros(len(Lat_EMEP))
    
    # Create a for loop to extract values from the models based on the latitude and longitude of the observations.
    for i in range(len(Lat_EMEP)):
        OLDval[i]= totwetdep_OLD.sel(lat=[Lat_EMEP[i]], lon=[Lon_EMEP[i]], method='nearest')
        NEWval[i]= totwetdep_NEW.sel(lat=[Lat_EMEP[i]], lon=[Lon_EMEP[i]], method='nearest')
    
    # Find the standard deviation of the values extracted from the reference and new model and round this value 
    # to 2 decimal places.
    ErrOLD= round(np.std(OLDval),1)
    ErrNEW= round(np.std(NEWval),1)

    # Round the means of the values extracted from the reference and new models to decimal places.
    MeMoOl=round(np.mean(OLDval),1)
    MeMoNe= round(np.mean(NEWval),1)

    # Find the correlation coefficient 
    corrOld, _ = stats.pearsonr(Value_EMEP, OLDval)
    corrNew, _ = stats.pearsonr(Value_EMEP, NEWval)    
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
    norm = colors.Normalize(vmin=0.,vmax=15.0)

    # Plot the reference model on the projection.
    im = totwetdep_OLD.plot.pcolormesh(x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                                       norm=norm, cmap='viridis',rasterized = True, zorder=0,
                                       cbar_kwargs={'orientation':'horizontal',
                                                    'ticklocation':'auto',
                                                    'label':"\u03BCg m$^{-2}$ yr$^{-1}$ ",
                                                    'extend': 'max',
                                                    'fraction':0.046,
                                                    'pad':0.04})
    # Show the coastlines.
    ax.coastlines()
    ax.add_feature(cf.BORDERS,zorder=10)
     
    # Add the observed values to the plot.
    plt.scatter(Lon_EMEP, Lat_EMEP,  transform=ccrs.PlateCarree(),marker='o',
                norm=norm,
                linewidths=0.5, edgecolors='black', 
                label=None, c=Value_EMEP, cmap='viridis', zorder=15)

    # Add a title.
    plt.title('Hg Wet Deposition, Reference Model (' + str(Year1) + '), EMEP (2013-15)' ,fontsize=15)   
  
    # Zoom into Europe area
    ax.set_xlim([-25, 36])
    ax.set_ylim([31, 71])
    
    
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
    norm = colors.Normalize(vmin=0.,vmax=15.0)

    # Plot the reference model on the projection.
    totwetdep_NEW.plot.pcolormesh( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                                  norm=norm, cmap='viridis',rasterized = True, zorder=0,
                                  cbar_kwargs={'orientation':'horizontal',
                                               'ticklocation':'auto',
                                               'label':"\u03BCg m$^{-2}$ yr$^{-1}$ ",
                                               'extend': 'max',
                                               'fraction':0.046,
                                               'pad':0.04})
    # Show the coastlines.
    ax.coastlines()
    ax.add_feature(cf.BORDERS,zorder=10)
     
    # Add the observed values to the plot.
    plt.scatter(Lon_EMEP, Lat_EMEP,  transform=ccrs.PlateCarree(),marker='o',
                norm=norm,
                linewidths=0.5, edgecolors='black', 
                label=None, c=Value_EMEP, cmap='viridis', zorder=15)

    # Add a title.
    plt.title('Hg Wet Deposition, New Model (' + str(Year2) + '), EMEP (2013-15)' ,fontsize=15)   
  
    # Zoom into Europe area
    ax.set_xlim([-25, 36])
    ax.set_ylim([31, 71])
    
    
    # Add text to the plot.
    plt.text(1.05, 0.1, textstr1, fontsize=13, transform=ax.transAxes)
    plt.text(1.05, 0.25, textstr3, fontsize=13, transform=ax.transAxes)
    plt.text(1.05, 0.4, textstr5, fontsize=13, transform=ax.transAxes)


    # Adjust to give text space
    NEWMAP.subplots_adjust(right=0.6, left=0.05)
    
    # Show the plot.
    NEWMAP.show()
   
    return OLDMAP, NEWMAP

