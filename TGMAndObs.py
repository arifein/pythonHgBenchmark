import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import pandas as pd
from scipy import stats
from SiteLevels import levels
from diff_plots_Hg import diff_plots
from matplotlib import colors

def SurfaceObsTGM(Old_Dataset, New_Dataset, Year = None):
    """ Plot the mean surface TGM for mercury against different sites for the reference and new models. Also calculate
    the mean for both models, the mean of the observations and the coefficient of determination. 
    
    Parameters
    ----------
    Dataset_OLD : string
        Reference Model xarray dataset
    Dataset_NEW : string
        New Model xarray dataset 
    
    Year : int or list of int, optional
        Optional parameter to only select subset of years
    
    """    
    # Read in the data for the observed sites
    AnHgObs= pd.read_csv('data/TGMSiteAnnual.csv',skiprows=[0], na_values=(-9999))
    AnHgObs.columns=['SiteID', 'Lat', 'Lon','Alt', 'TGM', 'Hg0']
    # Set levels for the colorbar in order to have a nonlinear scale.
    Levels= (0.75, 0.95, 1.15, 1.35, 1.55, 1.75, 2.30, 2.90, 3.50)
    SiteID=AnHgObs.SiteID
    
    # Make a variable for the unit conversion factor from vmr to  ng/m^3
    # Now more traceable
    R = 8.314462 # m^3 Pa K^-1 mol ^-1
    MW_Hg = 200.59 # g mol^-1
    ng_g = 1e9 # ng/g
    
    stdpressure = 101325 # Pascals
    stdtemp = 273.15 # Kelvins
    
    unit_conv = stdpressure / R / stdtemp * MW_Hg * ng_g # converter from vmr to ng m^-3
    
    SiteID=AnHgObs.SiteID
    
    # Allow subsetting for years, if inputted into the function
    if Year is not None: # take average over subset of years
        # OLD simulation        
        OLD_Hg0_yr = Old_Dataset.SpeciesConc_Hg0.sel(time=Old_Dataset.time.dt.year.isin(Year))
        OLD_Hg2_yr = Old_Dataset.SpeciesConc_Hg2.sel(time=Old_Dataset.time.dt.year.isin(Year))
        # NEW simulation        
        NEW_Hg0_yr = New_Dataset.SpeciesConc_Hg0.sel(time=New_Dataset.time.dt.year.isin(Year))
        NEW_Hg2_yr = New_Dataset.SpeciesConc_Hg2.sel(time=New_Dataset.time.dt.year.isin(Year))
    else: # use all years
        # OLD simulation        
        OLD_Hg0_yr = Old_Dataset.SpeciesConc_Hg0
        OLD_Hg2_yr = Old_Dataset.SpeciesConc_Hg2
        # NEW simulation                
        NEW_Hg0_yr = New_Dataset.SpeciesConc_Hg0
        NEW_Hg2_yr = New_Dataset.SpeciesConc_Hg2
   
    # Extract and add together Hg0 and Hg2 at the surface from the reference model multiplying by the unit converion factor 
    # to obtain values for Total Gaseous Mercury.
    OLD_Hg0 = OLD_Hg0_yr.isel(lev=0).mean('time') * unit_conv                 
    OLD_Hg2 = OLD_Hg2_yr.isel(lev=0).mean('time') * unit_conv

    NEW_Hg0 = NEW_Hg0_yr.isel(lev=0).mean('time') * unit_conv                 
    NEW_Hg2 = NEW_Hg2_yr.isel(lev=0).mean('time') * unit_conv
                       
    TGM_Old = (OLD_Hg0 + OLD_Hg2) # TGM is sum of Hg0 and Hg2
    TGM_New = (NEW_Hg0 + NEW_Hg2) # TGM is sum of Hg0 and Hg2
        
        
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
    StdObs= round(np.std(Value),2)

    # Create an array of numpy zeros to fill for the reference and new models based on the amount of observed values.
    OLDval=np.zeros(len(Lati))
    NEWval=np.zeros(len(Lati))
    
    # Create a for loop to extract values from the models based on the latitude and longitude of the observations.
    for i in range(len(Lati)):
        lev_site = levels(SiteID[i]) # get level of site
        if lev_site == 0: # surface site
            OLDval[i]= TGM_Old.sel(lat=[Lati[i]], lon=[Long[i]], method='nearest')
            NEWval[i]= TGM_New.sel(lat=[Lati[i]], lon=[Long[i]], method='nearest')
        else: # other level, need to reload TGM
            OLD_Hg0_lv = OLD_Hg0_yr.isel(lev=lev_site).\
                sel(lat=[Lati[i]], lon=[Long[i]], method='nearest').mean('time')\
                * unit_conv
            OLD_Hg2_lv = OLD_Hg2_yr.isel(lev=lev_site).\
                sel(lat=[Lati[i]], lon=[Long[i]], method='nearest').mean('time')\
                * unit_conv
            NEW_Hg0_lv = NEW_Hg0_yr.isel(lev=lev_site).\
                sel(lat=[Lati[i]], lon=[Long[i]], method='nearest').mean('time')\
                * unit_conv
            NEW_Hg2_lv = NEW_Hg2_yr.isel(lev=lev_site).\
                sel(lat=[Lati[i]], lon=[Long[i]], method='nearest').mean('time')\
                * unit_conv  
            # Calculate TGM values as sum of Hg0 and Hg2     
            OLDval[i]= OLD_Hg0_lv + OLD_Hg2_lv
            NEWval[i]= NEW_Hg0_lv + NEW_Hg2_lv

    # Take the mean of the values extracted from the reference and new models. 
    MeanModOld=(np.mean(OLDval))
    MeanModNew=(np.mean(NEWval))
    
    # Find the standard deviation of the values extracted from the reference and new model and round this value 
    # to 2 decimal places.
    ErrOLD= round(np.std(OLDval),2)
    ErrNEW= round(np.std(NEWval),2)

    # Round the means of the values extracted from the reference and new models to decimal places.
    MeMoOl=round(MeanModOld,2)
    MeMoNe= round(MeanModNew,2)

    # Find the correlation coefficient 
    corrOld, _ = stats.pearsonr(Value, OLDval)
    corrNew, _ = stats.pearsonr(Value, NEWval)    
    # Find the coefficient of determination for the reference and new models
    CoeffOld= corrOld ** 2
    CoeffNew= corrNew ** 2


    # Create text strings for relevant information: mean, coefficient of determination (rounding to 3DP),
    textstr1= "Mean Obs. = %s $\pm$ %s ng m$^{-3}$ "%(MeObsDP, StdObs)
    textstr2= "Mean Mod. = %s $\pm$ %s ng m$^{-3}$ "%(MeMoOl ,ErrOLD)
    textstr3= "Mean Mod. = %s $\pm$ %s ng m$^{-3}$ "%(MeMoNe ,ErrNEW)
    textstr4= "Terrestrial $R^2$= %s" %(round(CoeffOld,3))
    textstr5= "Terrestrial $R^2$= %s" %(round(CoeffNew,3))

    # Add a figure.
    OLDMAP = plt.figure(figsize=[10,4.5])
    
    # Add a geographical projection on the map.
    ax = OLDMAP.add_subplot(111, projection=ccrs.PlateCarree())
    
    # Plot the reference model on the projection.
    im = TGM_Old.plot.contourf( x='lon',y='lat', ax=ax, levels= Levels, transform=ccrs.PlateCarree(), cmap='viridis', 
                                 cbar_kwargs={'orientation':'horizontal',
                                           'ticklocation':'auto',
                                      'label':"Not Linear ng m$^{-3}$ ",
                                      'fraction':0.046,
                                      'pad':0.04})
     
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

    # Adjust to give text space
    OLDMAP.subplots_adjust(right=0.6, left=0.05)
    plt.tight_layout()
    
    # Show the plot.
    OLDMAP.show()
    

    
    # Add a figure.
    NEWMAP = plt.figure(figsize=[10,4.5])
    # Add a geographical projection on the map.
    ax = NEWMAP.add_subplot(111, projection=ccrs.PlateCarree())
    
    # Plot the new model on the projection.
    im=TGM_New.plot.contourf(x='lon',y='lat',levels=Levels, ax=ax,transform=ccrs.PlateCarree(), cmap='viridis', 
                         cbar_kwargs={'orientation':'horizontal',
                                      'ticklocation':'auto',
                                      'label':"Not Linear ng m$^{-3}$ ",
                                      'fraction':0.046,
                                      'pad':0.04})
    
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
    
    # Adjust to give text space
    NEWMAP.subplots_adjust(right=0.6, left=0.05)
    plt.tight_layout()
    # Show the plot.
    NEWMAP.show()  
    
    # Create difference plot for TGM at surface 
    diff_MAP = diff_plots(TGM_Old, TGM_New)
    
    # return all plots, for saving in PDF 
    return OLDMAP, NEWMAP, diff_MAP
