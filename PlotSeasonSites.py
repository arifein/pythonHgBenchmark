# Create a for loop in order to graph values for each unique site
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SiteLevels import levels   
from helper_functions import ds_sel_yr

def PlotSeasonSites(Dataset_OLD, Dataset_NEW, Year1 = None, Year2 = None):
    """ Plot the seasonal cycle of reference and new models against the TGM observations made at each site
    
    Parameters
    ----------
    Dataset_OLD : xarray dataset
        Reference Model dataset
    Dataset_NEW : xarray dataset
        New Model dataset 
    
    Year1 : int or list of int, optional
        Optional parameter to only select subset of years for old sim
    Year2 : int or list of int, optional
        Optional parameter to only select subset of years for new sim
        
    """
    # Import the observed data from the sites     
    Hgobs = pd.read_csv('data/TGMSiteMonthly.csv',  skiprows=[0], na_values=(-9999))
    Hgobs.columns=['SiteID', 'Lat', 'Lon','Month', 'Year', 'Concentration', 'Standard deviation']
    Site= Hgobs.SiteID
    
    # Arrange the data by order of latitude and ensure when graphs are plotted the data is sequential
    Graph_order=Hgobs.sort_values(by=['Lat'], ascending=0)
    HgobsOrder=Graph_order.sort_values(by=['Month'])

    # Calculate constant for the unit conversion factor from vmr to  ng/m^3
    R = 8.314462 # m^3 Pa K^-1 mol ^-1
    MW_Hg = 200.59 # g mol^-1
    ng_g = 1e9 # ng/g
    
    stdpressure = 101325 # Pascals
    stdtemp = 273.15 # Kelvins
    
    unit_conv = stdpressure / R / stdtemp * MW_Hg * ng_g # converter from vmr to ng m^-3
    
    # Create plot handles (42 sites for now...)
    SeasonGraph, axes = plt.subplots(6, 7, figsize=[16,12], 
                                     gridspec_kw=dict(hspace=0.6, wspace=0.7))
    axes = axes.flatten()
      
    # Create a loop that specifies unique site IDs and makes subplots
    for ii, SiteID in enumerate(Graph_order[ 'SiteID'].unique()):
        
        # Select all values with the same Site ID
        Dataset = HgobsOrder[HgobsOrder.SiteID == SiteID].reset_index()
        
    
        # Specify the first latitude and longitude of the Site ID 
        Lat=Dataset.Lat[0]
        Lon=Dataset.Lon[0]
        
        # Load TGM fields from model

        # Allow subsetting for years of the simulation, if inputted into the function
        OLD_Hg0_yr = ds_sel_yr(Dataset_OLD, 'vmrhg0', Year1)
        NEW_Hg0_yr = ds_sel_yr(Dataset_NEW, 'vmrhg0', Year2)
                  
        # Select level, lat, and longitude where data should be extracted from for both model runs
        OLD_Hg0_site = OLD_Hg0_yr.isel(lev=levels(SiteID)).\
            sel(lat=Lat, lon=Lon, method='nearest').squeeze()
        NEW_Hg0_site = NEW_Hg0_yr.isel(lev=levels(SiteID)).\
            sel(lat=Lat, lon=Lon, method='nearest').squeeze()
        
        # Consider TGM as only GEM and do unit conversion 
        OLD_TGM_site = (OLD_Hg0_site) * unit_conv
        NEW_TGM_site = (NEW_Hg0_site) * unit_conv
                    
        # Make sure that it is a climatology, if have multi-year model runs
        OLD_mod = OLD_TGM_site.groupby('time.month').mean()
        NEW_mod = NEW_TGM_site.groupby('time.month').mean()
    
        # Convert the time data from a float to a string, specifying months for graph labels
        Dataset.index=pd.to_datetime(Dataset.Month, format='%m')

        # Add a subplot 
        iax = axes[ii]
        
        # Add the data from the observations, the reference model and the new model
        iax.errorbar(Dataset.Month, Dataset.Concentration, Dataset['Standard deviation'],
                     color= "k")
        iax.plot(Dataset.Month,OLD_mod.data,color='blue')
        iax.plot(Dataset.Month,NEW_mod.data,color='red')
    
    
        # Label the axes, add a legend and add a title
        iax.set_ylabel('TGM (ng/m$^3$)')
        iax.set_title('{0} ({1}, {2})'.format(SiteID, Lat, Lon), fontsize=12)
    

        # Set ticks to every month 
        iax.set_xticks(Dataset.Month)
       
        # Set tick labels to month names
        mn = ['J','F','M','A','M','J','J','A','S','O','N','D']            
        iax.set_xticklabels(mn, fontsize=8)
        
        if ii>34:
            iax.set_xlabel('Month')
    
    SeasonGraph.legend(['Reference Model','New Model', 'Observations'],
                       loc = 'upper left') # add one legend to figure
    return SeasonGraph
