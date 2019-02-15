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
def SurfaceObsTGM(Old_Dataset, New_Dataset):
    
    # Read in the data for the observed sites
    AnHgObs= pd.read_csv('data/TGMSiteAnnual.csv',skiprows=[0], na_values=(-9999))
    AnHgObs.columns=['SiteID', 'Lat', 'Lon','Alt', 'TGM', 'Hg0']
    # Set levels for the colorbar in order to have a nonlinear scale.
    Levels= (0.75, 0.95, 1.15, 1.35, 1.55, 1.75, 2.30, 2.90, 3.50)
    SiteID=AnHgObs.SiteID
    
    # Make a variable for the unit conversion factor to obtain ng/m^3
    Unit_Conversion= 8.93
    SiteID=AnHgObs.SiteID
        # Extract and add together Hg0 and Hg2 at the surface from the reference model multiplying by the unit converion factor 
        # to obtain values for Total Gaseous Mercury.
    for i in range (len(SiteID)):
        OLD_Hg0 =((Old_Dataset['IJ_AVG_S_Hg0'].isel(lev=levels(SiteID[i])).mean('time')) * Unit_Conversion)                              
        OLD_Hg2 =((Old_Dataset['IJ_AVG_S_Hg2'].isel(lev=levels(SiteID[i])).mean('time')) * Unit_Conversion)                
        TGM_Old = (OLD_Hg0 + OLD_Hg2)
    
    
        # Extract and add together Hg0 and Hg2 at the surface from the new model multiplying by the unit converion factor 
        # to obtain values for Total Gaseous Mercury.
        NEW_Hg0 =((New_Dataset['IJ_AVG_S_Hg0'].isel(lev=levels(SiteID[i])).mean('time') * Unit_Conversion))                         
        NEW_Hg2 =((New_Dataset['IJ_AVG_S_Hg2'].isel(lev=levels(SiteID[i])).mean('time') * Unit_Conversion))
        TGM_New = NEW_Hg0 + NEW_Hg2
    
        
    # Find the absolute difference between the reference and new model.
    Abs_diff = TGM_New - TGM_Old
    # Find the absolute maximum value of the absolute difference. 
    Abs_MaxVal= np.max(np.abs(Abs_diff))
        
    
    # Find the percent difference of the models.  
    Perc_diff = (Abs_diff / TGM_Old)*100
    # Find the absolute maximum value of the percent  difference. 
    Perc_MaxVal= np.max(np.abs(Perc_diff))
    
    
    
    # Set variable names for the longitude and latitude in the dataset.
    Long=(AnHgObs['Lon'])
    Lati=(AnHgObs['Lat'])
    
    # Set a variable name for the observed values of Hg and find the mean of these values. 
    Value= AnHgObs['Hg0']
    Meanobvs=np.mean(Value)
    
    # Round the mean of the observations to 2 decimal places. 
    MeObsDP=round(Meanobvs,2)
    
    # Find the standard deviation of the observed Hg and round this value to 2 decimal places.
    ErrObs= round(statistics.stdev(Value),2)

    # Create an array of numpy zeros to fill for the reference and new models based on the amount of observed values.
    OLDval=np.zeros(len(Lati))
    NEWval=np.zeros(len(Lati))
    
    # Create a for loop to extract values from the models based on the latitude and longitude of the observations.
    for i in range(len(Lati)):
        OLDval[i]= (TGM_Old.sel(lat=[Lati[i]], lon=[Long[i]], method='nearest'))
        NEWval[i]= (TGM_New.sel(lat=[Lati[i]], lon=[Long[i]], method='nearest'))


    # Take the mean of the values extracted from the reference and new models. 
    MeanModOld=(np.mean(OLDval))
    MeanModNew=(np.mean(NEWval))
    
    # Find the standard deviation of the values extracted from the reference and new model and round this value 
    # to 2 decimal places.
    ErrOLD= round(statistics.stdev(OLDval),2)
    ErrNEW= round(statistics.stdev(NEWval),2)

    # Round the means of the values extracted from the reference and new models to decimal places.
    MeMoOl=round(MeanModOld,2)
    MeMoNe= round(MeanModNew,2)

    # Find the coefficient of determination for the reference and new models/
    CoeffOld= r2_score(Value, OLDval)
    CoeffNew= r2_score(Value, NEWval)


    # Create text strings for relevant information: mean, coefficient of determination (rounding to 3DP),
    textstr1= "Mean Obs. = %s +- %s $ng/m^3$ "%(MeObsDP, ErrObs)
    textstr2= "Mean Mod. = %s +- %s $ng/m^3$ "%(MeMoOl ,ErrOLD)
    textstr3= "Mean Mod. = %s +- %s $ng/m^3$ "%(MeMoNe ,ErrNEW)
    textstr4= "Terrestrial $R^2$= %s" %(round(CoeffOld,3))
    textstr5= "Terrestrial $R^2$= %s" %(round(CoeffNew,3))



    # Add a figure.
    OLDMAP = plt.figure()
    
    # Add a geographical projection on the map.
    ax = OLDMAP.add_subplot(111, projection=ccrs.PlateCarree())
    
    # Plot the reference model on the projection.
    im=TGM_Old.plot.contourf( x='lon',y='lat', ax=ax, levels= Levels, transform=ccrs.PlateCarree(), cmap='viridis', 
                                 cbar_kwargs={'orientation':'horizontal',
                                           'ticklocation':'auto',
                                      'label':"Not Linear $ng/m^3$ "})
    
    # Add text to the plot.
    plt.text(200,-50,textstr1, fontsize=14)
    plt.text(200,-75,textstr2, fontsize=14)
    plt.text(200,-100,textstr4, fontsize=14)

    # Add the observed values to the plot.
    plt.scatter(Long, Lati,  transform=ccrs.PlateCarree(),marker='D',
                norm=colors.BoundaryNorm(boundaries=Levels, ncolors=256), 
                linewidths=0.5, edgecolors='black',
                label=None, c=Value, cmap='viridis')

    # Add a title.
    plt.title(' Reference Model Version: Surface TGM',fontsize=15)   
    
    # Show the coastlines.
    ax.coastlines()
    
    # Show the plot.
    plt.show()
    

    
    # Add a figure.
    NEWMAP = plt.figure()
    
    # Add a geographical projection on the map.
    ax = NEWMAP.add_subplot(111, projection=ccrs.PlateCarree())
    
    # Plot the new model on the projection.
    im=TGM_New.plot.contourf(x='lon',y='lat',levels=Levels, ax=ax,transform=ccrs.PlateCarree(), cmap='viridis', 
                         cbar_kwargs={'orientation':'horizontal',
                                      'ticklocation':'auto',
                                      'label':"Not Linear $ng/m^3$ "})
    
    # Add text to the plot.
    plt.text(200,-50,textstr1, fontsize=14)
    plt.text(200,-75,textstr3, fontsize=14)
    plt.text(200,-100,textstr5, fontsize=14)
    
    # Add the observed values to the plot.
    plt.scatter(Long, Lati,  transform=ccrs.PlateCarree(),marker='D',
                norm=colors.BoundaryNorm(boundaries=Levels, ncolors=256), 
                linewidths=0.75, edgecolors='black',
                label=None, c=Value, cmap='viridis')

    # Add a title.
    plt.title(' New Model Version: Surface TGM', fontsize=15)       
    
    # Show the coastlines.
    ax.coastlines()
    
    # Show the plot.
    plt.show()  
    
    return (OLDMAP, NEWMAP)

    
    # Add a figure.
    NEWMAP = plt.figure()
    
    # Add a geographical projection on the map.
    ax = NEWMAP.add_subplot(111, projection=ccrs.PlateCarree())
    
    # Plot the new model on the projection.
    im=TGM_New.plot.contourf(x='lon',y='lat',levels=Levels, ax=ax,transform=ccrs.PlateCarree(), cmap='viridis', 
                         cbar_kwargs={'orientation':'horizontal',
                                      'ticklocation':'auto',
                                      'label':"Not Linear $ng/m^3$ "})
    
    # Add text to the plot.
    plt.text(200,-50,textstr1, fontsize=14)
    plt.text(200,-75,textstr3, fontsize=14)
    plt.text(200,-100,textstr5, fontsize=14)
    
    # Add the observed values to the plot.
    plt.scatter(Long, Lati,  transform=ccrs.PlateCarree(),marker='D',
            norm=colors.BoundaryNorm(boundaries=Levels, ncolors=256), 
            linewidths=0.75, edgecolors='black',
            label=None, c=Value, cmap='viridis')

    # Add a title.
    plt.title(' New Model Version: Surface TGM', fontsize=15)       
    
    # Show the coastlines.
    ax.coastlines()
    
    # Show the plot.
    plt.show()  
    
    return (OLDMAP, NEWMAP)
