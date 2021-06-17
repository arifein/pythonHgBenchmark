import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SiteLevels import levels
from scipy.io import readsav
import sys
from helper_functions import ds_sel_yr, annual_avg

def Seasonal_Lat_Regions(Dataset_OLD, Dataset_NEW, Year = None):
    """Plot observational seasonal cycle against the model for different 
    latitudes (Southern Mid Latitiude, North Mid Latitude,
    Arctic, Antarctic).
    
    Parameters
    ----------
    Dataset_OLD : xarray dataset
        Reference Model dataset
    Dataset_NEW : xarray dataset
        New Model dataset 
    
    Year : int or list of int, optional
        Optional parameter to only select subset of years    
    
    """

    # Import the observed data from the sites     
    Hgobs = pd.read_csv('data/TGMSiteMonthly.csv',  skiprows=[0], na_values=(-9999))
    Hgobs.columns=['SiteID', 'Lat', 'Lon','Month', 'Year', 'Concentration', 'Standard deviation']
 
    # Make arrays of SiteIDs for the Arctic, Antarctic and Northern and Souther Mid Latitudes, extracting the sites
    # in the data sets. 
    Arctic = ['ALT', 'VRS', 'ZEP', 'AND', 'PAL','AMD']
    Arctic = [e for e in Arctic if e in list(Hgobs.SiteID)]
    
    SouthMidLat=['CPT', 'AMS', 'BAR']
    SouthMidLat = [e for e in SouthMidLat if e in list(Hgobs.SiteID)]
    
    Antarctic= ['TRO', 'DDU', 'DMC']
    Antarctic = [e for e in Antarctic if e in list(Hgobs.SiteID)]
    
    NorthMidLat= ['MHD', 'UDH', 'KEJ',  'HTW', 'PNY', 'ATN', 'YKV', 'GRB','BIR', 'WAL', 'BRA', 'SAT', 'THOMPFARM', 'SCO', 'STIWELL', 'EBG'] 
    NorthMidLat = [e for e in NorthMidLat if e in list(Hgobs.SiteID)]
    
    # Calculate mean, std for model and obs in each region
    
    # Arctic
    Arc_df = filter_sites_region(Arctic, Hgobs, Dataset_OLD, Dataset_NEW, Year)
    # Antarctic                  
    Ant_df = filter_sites_region(Antarctic, Hgobs, Dataset_OLD, Dataset_NEW, Year)
    # Northern Mid Latitudes                  
    NML_df = filter_sites_region(NorthMidLat, Hgobs, Dataset_OLD, Dataset_NEW, Year)
    # Southern Mid Latitudes
    SML_df = filter_sites_region(SouthMidLat, Hgobs, Dataset_OLD, Dataset_NEW, Year)

    # Create a list of all regions for looped plots
    Regions_df_all = [Arc_df, Ant_df, NML_df, SML_df] 
    Region_names = ['Arctic', 'Antarctic','Northern Mid Latitudes', 'Southern Mid Latitudes']
     
    # Plot the four graphs as subplots.
    RegPlot,  axes = plt.subplots(2, 2, figsize=[16,12],
                                   gridspec_kw=dict(hspace=0.3, wspace=0.2))
    axes = axes.flatten()
    
    # Loop over regions and plot
    for ii, iax in enumerate(axes):
        
        Reg_df = Regions_df_all[ii] # dataframe with data
        # Plot the observations and their error.
        iax.errorbar(Reg_df.index, Reg_df['Obs_mean'], yerr=Reg_df['Obs_std'], 
                     color='k', capsize=4)
        # Plot the reference and new models on the same graph with their errors.
        iax.errorbar(Reg_df.index, Reg_df['OLD_mean'], yerr=Reg_df['OLD_std'], 
                     color='Blue', capsize=4)
        iax.errorbar(Reg_df.index, Reg_df['NEW_mean'], yerr=Reg_df['NEW_std'], 
                     color='Red', capsize=4)
        # Label the x and y axis. 
        iax.set_xlabel('Month')
        iax.set_ylabel('TGM (ng/m$^3$)')
        if ii==0: # add legend only for first plot
            iax.legend([ 'Observations','Reference Model','New Model' ])
        # Add a title.
        iax.set_title(Region_names[ii])
        # Set ticks to every month 
        iax.set_xticks(Reg_df.index)
        # Set tick labels to month names
        mn = ['J','F','M','A','M','J','J','A','S','O','N','D']
        iax.set_xticklabels(mn)
    
    return RegPlot

def filter_sites_region(Region, Hgobs, Dataset_OLD, Dataset_NEW, Year = None):
    """ Calculate the regional seasonal cycle in observations and the two model simulations.
    
    Parameters
    ----------
    Region : list
        List of strings giving the site names to average over    
    Hgobs : DataFrame
        Observational dataset for seasonal cycle
    Dataset_OLD : xarray dataset
        Reference Model dataset
    Dataset_NEW : xarray dataset
        New Model dataset 
    """

    # Calculate constant for the unit conversion factor from vmr to  ng/m^3
    R = 8.314462 # m^3 Pa K^-1 mol ^-1
    MW_Hg = 200.59 # g mol^-1
    ng_g = 1e9 # ng/g
    
    stdpressure = 101325 # Pascals
    stdtemp = 273.15 # Kelvins
    
    unit_conv = stdpressure / R / stdtemp * MW_Hg * ng_g # converter from vmr to ng m^-3

        
    # Extract the data from each site in the region, creating a new DataFrame
    for i, isite in enumerate(Region):
      site_df = Hgobs[Hgobs.SiteID==isite].reset_index()
      if i==0:
          All_region_df = site_df
      else:
          All_region_df = pd.concat([All_region_df, site_df])
          
    # Calculate the mean and stanadard deviation of observations for each month.
    obs_mean = np.asarray(All_region_df.groupby('Month').mean().Concentration)
    obs_std = np.asarray(All_region_df.groupby('Month').std().Concentration)
    
    # Select all unique latitudes and longitudes from the dataset.
    obs_lat = All_region_df.Lat.unique()
    obs_lon = All_region_df.Lon.unique()
    
    # Allow subsetting for years of the simulation, if inputted into the function
    OLD_Hg0_yr = ds_sel_yr(Dataset_OLD, 'SpeciesConc_Hg0', Year)
    OLD_Hg2_yr = ds_sel_yr(Dataset_OLD, 'SpeciesConc_Hg2', Year)
    NEW_Hg0_yr = ds_sel_yr(Dataset_NEW, 'SpeciesConc_Hg0', Year)
    NEW_Hg2_yr = ds_sel_yr(Dataset_NEW, 'SpeciesConc_Hg2', Year)
    
    # Create datasets for seasonal TGM at each site for the ref and new models     
    for i in range (len(Region)): 
      OLD_Hg0_site = OLD_Hg0_yr.isel(lev=levels(Region[i])).\
          sel(lat=[obs_lat[i]], lon=[obs_lon[i]], method='nearest').squeeze()
      OLD_Hg2_site = OLD_Hg2_yr.isel(lev=levels(Region[i])).\
          sel(lat=[obs_lat[i]], lon=[obs_lon[i]], method='nearest').squeeze()
      NEW_Hg0_site = NEW_Hg0_yr.isel(lev=levels(Region[i])).\
          sel(lat=[obs_lat[i]], lon=[obs_lon[i]], method='nearest').squeeze()
      NEW_Hg2_site = NEW_Hg2_yr.isel(lev=levels(Region[i])).\
          sel(lat=[obs_lat[i]], lon=[obs_lon[i]], method='nearest').squeeze()  
          
      # Calculate TGM values as sum of Hg0 and Hg2     
      Reg_OLD_mod = (OLD_Hg0_site + OLD_Hg2_site) * unit_conv
      Reg_NEW_mod = (NEW_Hg0_site + NEW_Hg2_site) * unit_conv

      # calculate climatology (needed if more than one year are averaged)
      Reg_OLD_clim = Reg_OLD_mod.groupby('time.month').mean() 
      Reg_NEW_clim = Reg_NEW_mod.groupby('time.month').mean() 
   
      if i==0:
          Reg_DS_OLD = Reg_OLD_clim
          Reg_DS_NEW = Reg_NEW_clim
      else: # concatenate site values together
          Reg_DS_OLD= xr.concat([Reg_DS_OLD,Reg_OLD_clim], dim='concat_dims')
          Reg_DS_NEW= xr.concat([Reg_DS_NEW,Reg_NEW_clim], dim='concat_dims')
          
    # Calculate the mean and standard deviations for the reference and new models.
    OLD_mean = np.asarray(Reg_DS_OLD.mean('concat_dims'))
    OLD_std = np.asarray(Reg_DS_OLD.std('concat_dims'))
      
    NEW_mean = np.asarray(Reg_DS_NEW.mean('concat_dims'))
    NEW_std = np.asarray(Reg_DS_NEW.std('concat_dims'))
    
    # Save results in a Pandas DataFrame
    data_dic = {'Obs_mean': obs_mean,'Obs_std': obs_std, 
                'NEW_mean': NEW_mean,'NEW_std': NEW_std,
                'OLD_mean': OLD_mean,'OLD_std': OLD_std}
    Out_df = pd.DataFrame(data_dic)
    
    return Out_df

def plot_gradient_TGM(Dataset_OLD, Dataset_NEW, Year = None):
    """Plot meridional gradient of TGM at surface from netcdf files of model simulation
    
    Parameters
    ----------
    Dataset_OLD : xarray dataset
        Reference Model dataset
    Dataset_NEW : xarray dataset
        New Model dataset 
    
    Year : int or list of int, optional
        Optional parameter to only select subset of years        
    """
    # First load surface TGM from model simulations
    # Make a variable for the unit conversion factor from vmr to  ng/m^3
    # Now more traceable
    R = 8.314462 # m^3 Pa K^-1 mol ^-1
    MW_Hg = 200.59 # g mol^-1
    ng_g = 1e9 # ng/g
    
    stdpressure = 101325 # Pascals
    stdtemp = 273.15 # Kelvins
    
    unit_conv = stdpressure / R / stdtemp * MW_Hg * ng_g # converter from vmr to ng m^-3
        
    # Allow subsetting for years, if inputted into the function
    OLD_Hg0_yr = ds_sel_yr(Dataset_OLD, 'SpeciesConc_Hg0', Year)
    OLD_Hg2_yr = ds_sel_yr(Dataset_OLD, 'SpeciesConc_Hg2', Year)
    NEW_Hg0_yr = ds_sel_yr(Dataset_NEW, 'SpeciesConc_Hg0', Year)
    NEW_Hg2_yr = ds_sel_yr(Dataset_NEW, 'SpeciesConc_Hg2', Year)
   
    # Extract and add together Hg0 and Hg2 at the surface from both 
    # model simulations, multiplying by the unit conversion factor 
    # to obtain values for Total Gaseous Mercury.
    OLD_Hg0 = annual_avg(OLD_Hg0_yr.isel(lev=0))               
    OLD_Hg2 = annual_avg(OLD_Hg2_yr.isel(lev=0))

    NEW_Hg0 = annual_avg(NEW_Hg0_yr.isel(lev=0))                 
    NEW_Hg2 = annual_avg(NEW_Hg2_yr.isel(lev=0))
                       
    TGM_Old = (OLD_Hg0 + OLD_Hg2) * unit_conv # TGM is sum of Hg0 and Hg2
    TGM_New = (NEW_Hg0 + NEW_Hg2) * unit_conv # TGM is sum of Hg0 and Hg2
    
    # Calulate zonal mean of the model results
    TGM_Old_z = TGM_Old.mean('lon')
    TGM_New_z = TGM_New.mean('lon')

    # Load observational datasets
    
    # Read in the data for the Land-Based TGM
    AnHgObs= pd.read_csv('data/TGMSiteAnnual.csv',skiprows=[0], na_values=(-9999))
    AnHgObs.columns=['SiteID', 'Lat', 'Lon','Alt', 'TGM', 'Hg0']
    
    # Read in the cruise data
    fn_cruise = 'data/tgm_gradient.sav'
    temme_lat, temme_tgm = load_sav_data(fn_cruise, 'temme')
    lamborg_lat, lamborg_tgm = load_sav_data(fn_cruise, 'lamborg')
    laurier_lat, laurier_tgm = load_sav_data(fn_cruise, 'laurier2007')
    soerenson_lat, soerenson_tgm = load_sav_data(fn_cruise, 'galathea') # Soerenson et al. 2010, ES&T
    
    # Plot zonal mean of model vs. observations
    Grad_fig = plt.figure(figsize=[8,9])
    
    # Model lines
    plt.plot(Dataset_OLD.lat, TGM_Old_z, color='blue')
    plt.plot(Dataset_NEW.lat, TGM_New_z, color='red')
    
    # Observational sites
    plt.plot(AnHgObs['Lat'], AnHgObs['Hg0'], 'gs')
    plt.plot(temme_lat, temme_tgm, 'k+')
    plt.plot(laurier_lat, laurier_tgm, 'kd', mfc="None")
    plt.plot(lamborg_lat, lamborg_tgm, 'ks', mfc="None")
    plt.plot(soerenson_lat, soerenson_tgm, 'ko', mfc="None")
    
    # Add a title and axes labels
    plt.title('Surface TGM', fontsize=15)       
    plt.xlabel('Latitude', fontsize=13)       
    plt.ylabel('TGM (ng/m$^{3}$)', fontsize=13)       

    # Add a legend
    plt.legend(['Reference Model Zonal Mean','New Model Zonal Mean', 'Land-based stations', 
                'Temme et al. 2003', 'Laurier et al. 2007', 'Lamborg et al. 1999',
                'Soerenson et al. 2010'],
               fontsize=13)
    
    return Grad_fig

def load_sav_data(fn, ds_name):
    """Read .sav files and return the latitudinal gradient of TGM for a specific dataset
    
    Parameters
    ----------
    fn : str
        Filename of .sav file with TGM data
    ds_name : str
        Dataset to return (Options: "temme", "lamborg", "laurier2007", "galathea") 
    
    """
    possible_ds = ["temme", "lamborg", "laurier2007", "galathea", 
                   "TEMME", "LAMBORG", "LAURIER2007", "GALATHEA"]
    if ds_name not in possible_ds:
        print("Dataset name is not a possible option in .sav file")
        sys.exit(1)
    sav_data = readsav(fn,  python_dict=True) # read .sav file that has cruise data
    ds_data = sav_data['gradient'][ds_name][0]
    
    if ds_name.lower() == 'galathea': # different name for TGM
        ds_data_tgm = ds_data['hg0'][0]
        ds_data_lat = ds_data['lat'][0]
        # Lat only has up to size 44, filter ds_data_tgm array up to 44
        n_lat = len(ds_data_lat)
        ds_data_tgm = ds_data_tgm[:n_lat]
    else:
        ds_data_tgm = ds_data['tgm'][0]
        ds_data_lat = ds_data['lat'][0]
        # check where have invalid results for tgm and remove from dataset
        bool_z = ds_data_tgm!=0 
        ds_data_tgm = ds_data_tgm[bool_z]
        ds_data_lat = ds_data_lat[bool_z]
    
    # Calculate constant for the unit conversion factor from pptv to  ng/m^3
    R = 8.314462 # m^3 Pa K^-1 mol ^-1
    MW_Hg = 200.59 # g mol^-1
    ng_g = 1e9 # ng/g
    ppt_vmr = 1e-12 # ppt to vmr

    stdpressure = 101325 # Pascals
    stdtemp = 273.15 # Kelvins

    unit_conv = stdpressure / R / stdtemp * MW_Hg * ng_g * ppt_vmr # converter from ppt to vmr

    ds_data_tgm = ds_data_tgm * unit_conv # convert to ng/m^3

    return ds_data_lat, ds_data_tgm



