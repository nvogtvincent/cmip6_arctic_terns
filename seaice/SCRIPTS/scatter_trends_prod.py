# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 16:39:28 2022

@author: IsoldeGlissenaar

Creates Supplementary Figure 8
"""

import pandas as pd
import xarray as xr
import numpy as np
from scipy import stats, spatial
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m')
import matplotlib.gridspec as gridspec

# Set some plotting defaults
plt.rcParams['figure.figsize'] = (6, 4)
plt.rcParams['figure.dpi'] = 200

#%%
''' **Read data** '''
# **Read observed (OISST) sea ice concentration data**
data_directory = '../DATA/'
dataset = xr.open_dataset(data_directory+'icec.mnmean.nc')
x,y = np.meshgrid(dataset.lon.values, dataset.lat.values)

# Productivity
prod = xr.open_dataset(data_directory+'global-reanalysis-bio-intnpp_011993-112019.nc')
prod_lat, prod_lon = np.meshgrid(prod.LATITUDE.values, prod.LONGITUDE.values)
coord_prod = np.transpose(np.array([prod_lat.flatten(), prod_lon.flatten()]))
intnpp = prod.INTNPP.values.reshape(323,289440)
time_prod = pd.DatetimeIndex(prod.TIME.values)
    
#%%
#Plot decadal sea ice edge for given month 
months = [10,11,12,1,2,3]
prod_sie_save = np.zeros((len(months),27))

for k in range(len(months)):

    month = months[k]
    if month==12:
        years = np.arange(1993,2019,1)
    else:
        years = np.arange(1993,2020,1)
    prod_sie = np.zeros(len(years))
    
    for yr in range(len(years)):
        #Get contour sea ice edge
        fig1=plt.figure(2, dpi=200)
        ax = plt.axes(projection=ccrs.Orthographic(central_longitude=0, central_latitude=-90, globe=None))
        contour = ax.contour(x, y, dataset.icec.sel(time=(str(years[yr])+'-'+str(month)+'-01')), levels=[15], linewidths=0.5, transform=ccrs.PlateCarree())
        
        # Get ice edge 
        p1 = contour.collections[0].get_paths()  # grab the 1st path
        coor_p1 = np.zeros((0,2))
        for i in range(len(p1)):
            coor_p1 = np.append(p1[i].vertices, coor_p1, axis=0)
        
        coor_p1 = coor_p1[coor_p1[:,1]<-57.5]
        
        edge = coor_p1[(coor_p1[:,0]<121)|(coor_p1[:,0]>300)]
        edge[(edge[:,0]>180),0] = edge[(edge[:,0]>180),0]-360
        edge = edge[edge[:,0].argsort()] 
        
        #Get productivity at ice edge
        idx_prod_yr = np.where((time_prod.year==years[yr])&(time_prod.month==month))
        
        coord = np.transpose(np.array([edge[:,1], edge[:,0]]))
        MdlKDT = spatial.KDTree(coord_prod)
        dist, closest = MdlKDT.query(coord, k=1, distance_upper_bound=0.5)
        closest[closest==len(coord_prod)] = 0
        
        prod_sie[yr] = np.nanmean(intnpp[idx_prod_yr[0][0],closest])
        
    if month==12:
        prod_sie_save[k,:-1] = prod_sie
        prod_sie_save[k,-1] = np.nan
    else:
        prod_sie_save[k,:] = prod_sie
        
#%% 
#Make figure (Supplementary figure 8)
letter = ['a','b','c','d','e','f']
fig = plt.figure(figsize=(16,10))
gs1 = gridspec.GridSpec(2, 3)
gs1.update(wspace=0.05, hspace=0.05)

for k in range(len(months)):
    #Get trend in productivity over time
    if k==2:
        slope, intercept, r_value, p_value, std_err = stats.linregress(years[:-1], prod_sie_save[k,:-1])
    else:
        slope, intercept, r_value, p_value, std_err = stats.linregress(years, prod_sie_save[k,:])
    
    ax1 = plt.subplot(gs1[k])
    ax1.scatter(years, prod_sie_save[k,:])
    ax1.plot(years, slope*years+intercept, color='black')
    ax1.set_ylim([200,410])
    if (k==0)|(k==3):
        ax1.set_ylabel('mean productivity at ice edge [g C m$^{-2}$ day$^{-1}$]')
    else:
        ax1.axes.yaxis.set_ticklabels([])
    if k<3:
        ax1.axes.xaxis.set_ticklabels([])
    ax1.grid(linestyle=':')
    
    ax1.text(1992.5,210,letter[k],fontsize=18, weight='bold')

plt.show()

