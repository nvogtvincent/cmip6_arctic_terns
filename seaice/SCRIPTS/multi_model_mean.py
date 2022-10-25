# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 14:59:43 2022

@author: zq19140
"""

#Plots the change in decadal mean sea ice extent (1981-2020) from observations and CMIP6.

import glob
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m')

month=1

#%%

''' **Observed sea ice extent** '''
### **Read ice concentration data**
#NOAA OISST ice concentration 1981-2021, retrieved from:
#https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html

data_directory = '../DATA/'

dataset = xr.open_dataset(data_directory+'icec.mnmean.nc').sel(time=slice('1981','2014'))
x,y = np.meshgrid(dataset.lon.values, dataset.lat.values)

#%%
''' **Historical Modelled sea ice edge** '''
#CMIP6 data too large to add to GitHub

direc = '../DATA/SImon_historical/'
filenames = glob.glob(direc+'*.nc')

siconc_hist = np.full((29,34,180,360),np.nan)
i=0
for f in filenames:
    cmip = xr.open_dataset(f).sel(time=slice('1981','2014'))
    cmip = cmip.sel(time=cmip.time.dt.month.isin([month]))
    try:
        siconc_hist[i,:,:,:] = cmip.siconc
    except ValueError:
        siconc_hist[i,:-1,:,:] = cmip.siconc
    i = i+1
    
model_mean_hist = np.nanmean(siconc_hist,axis=0)
model_mean_hist[model_mean_hist==0] = np.nan

year_hist = np.zeros(len(cmip.time.values))
for i in range(len(cmip.time.values)):
    year_hist[i] = cmip.time.values[i].year
    
x_hist,y_hist = np.meshgrid(cmip.ETOPO60X, cmip.ETOPO60Y)  


#%%

''' Plot ice extent change '''

fig=plt.figure(figsize=(5,5), dpi=600)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=0, central_latitude=-90, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-180,180,-50,-90],crs=ccrs.PlateCarree()) 
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')

data_month = dataset.groupby('time.month').mean('time')
obs = plt.contour(x, y, data_month.icec.sel(month=month), levels=[15], colors='C0', linewidths=0.8, label='observed 1981-2014', transform=ccrs.PlateCarree())
for i in range(siconc_hist.shape[0]):
    plt.contour(x_hist, y_hist, np.nanmean(siconc_hist[i,:,:,:],axis=0), levels=[15], colors='grey', linewidths=0.2, transform=ccrs.PlateCarree(),zorder=0)
histo = plt.contour(x_hist, y_hist, np.nanmean(model_mean_hist[:,:,:],axis=0), levels=[15], colors='black', linewidths=0.8, label='CMIP historical 1981-2014', transform=ccrs.PlateCarree(),zorder=0)

ax.scatter(0,-90,s=5000,c='white', transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')

ax.text(-135, -32,'c',weight='bold',fontsize=12,transform=ccrs.PlateCarree())

legends = []
legends.extend(obs.collections)
legends.extend(histo.collections)
labels = ['observed 1981-2014', 'historical 1981-2014']#, 'SSP5 8.5 2081-2100']
ax.legend(legends, labels, bbox_to_anchor=(0.5, -0.02), loc='upper center',fontsize=10,ncol=3,frameon=False)


