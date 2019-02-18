# Create a for loop in order to graph values for each unique site
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import pandas as pd
import xbpch
import cartopy.crs as ccrs
from matplotlib import colorbar, colors
import statistics
from sklearn.metrics import r2_score
from SiteLevels import levels   
def PlotSeasonSites(Dataset_OLD, Dataset_NEW):
    """ Plot the reference and new models against the observations made at each site for a year. 
    
    Args:
    Dataset_OLD (str) : Reference Model bpch file
    Dataset_NEW (str) : New Model bpch file     
    
    """
    # Import the observed data from the sites     
    Hgobs = pd.read_csv('data/TGMSiteMonthly.csv',  skiprows=[0], na_values=(-9999))
    Hgobs.columns=['SiteID', 'Lat', 'Lon','Month', 'Year', 'Concentration', 'Standard deviation']
    Site= Hgobs.SiteID
    
    # Arrange the data by order of latitude and ensure when graphs are plotted the data is sequential
    Graph_order=Hgobs.sort_values(by=['Lat'], ascending=0)
    HgobsOrder=Graph_order.sort_values(by=['Month'])

    # Make a variable for the unit conversion factor to obtain ng/m^3
    Unit_Conversion= 8.93
    
    # Create a loop that specifies unique site IDs
    for SiteID in Graph_order[ 'SiteID'].unique():
        
        # Select all values with the same Site ID
        Dataset = HgobsOrder[HgobsOrder.SiteID == SiteID].reset_index()
        
        # Choose the first Site ID that is unique
        Site_ID= Dataset.SiteID[0]
        
        
        # Extract and add together Hg0 and Hg2 at the surface from the reference model multiplying by the unit converion factor 
        # to obtain values for Total Gaseous Mercury.
        OLD_Hg0 =((Dataset_OLD['IJ_AVG_S_Hg0'].isel(lev=levels(SiteID))) * Unit_Conversion)                              
        OLD_Hg2 =((Dataset_OLD['IJ_AVG_S_Hg2'].isel(lev=levels(SiteID))) * Unit_Conversion)                
        TGM_Old = (OLD_Hg0 + OLD_Hg2)




   # Extract and add together Hg0 and Hg2 at the surface from the new model multiplying by the unit converion factor 
        # to obtain values for Total Gaseous Mercury.
        NEW_Hg0 =((Dataset_NEW['IJ_AVG_S_Hg0'].isel(lev=levels(SiteID)) * Unit_Conversion))                         
        NEW_Hg2 =((Dataset_NEW['IJ_AVG_S_Hg2'].isel(lev=levels(SiteID)) * Unit_Conversion))
        TGM_New =( NEW_Hg0 + NEW_Hg2)
    
    
        # Specify the first latitude and longitude of the Site ID 
        Lat=Dataset.Lat[0]
        Lon=Dataset.Lon[0]
        

        # Specify the latitude and longitude where data should be extracted from for both reference and new models.
        OLD_mod= (TGM_Old.sel(lat=Lat, lon=Lon,  method='nearest'))
        NEW_mod= (TGM_New.sel(lat=Lat, lon=Lon,  method='nearest'))   
    
        # Convert the time data from a float to a string, specifying months for graph labels
        Dataset.index=pd.to_datetime(Dataset.Month, format='%m')

        
        # Add a plot 
        SeasonGraph= plt.figure()
        
    
        # Add the data from the observations, the reference model and the new model
        ax=Dataset.plot(x='Month', y='Concentration',yerr='Standard deviation' ,color= "k")
        ax.plot(Dataset.Month,OLD_mod.data,color='blue')
        ax.plot(Dataset.Month,NEW_mod.data,color='red')
    
    
        # Label the axes, add a legend and add a title
        plt.xlabel('Month')
        plt.ylabel('TGM (ng/m3)')
        plt.legend([ 'Reference Model','New Model', 'Observations'])
        plt.title('{0} ({1}, {2})'.format(SiteID, Lat, Lon), fontsize=15)
    
    
        # Set ticks to every month 
        ax.set_xticks(Dataset.Month)
    
    
        # Set tick labels to month names
        ax.set_xticklabels(Dataset.index.strftime('%b'))
    
    
        # Show the plot
        plt.show()
    return SeasonGraph
