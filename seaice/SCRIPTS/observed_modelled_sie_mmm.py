# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 14:07:33 2021

@author: zq19140
"""

''' **Observed and modelled sea ice extent** '''

#Plots the change in decadal mean sea ice extent (1981-2020) from observations and CMIP6.
#CMIP6 data too large to add to GitHub

import glob
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m')
import matplotlib
import matplotlib.gridspec as gridspec
import cmocean

# Set some plotting defaults
plt.rcParams['figure.dpi'] = 200

#%%

''' **Observed sea ice extent** '''
### **Read ice concentration data**
#NOAA OISST ice concentration 1981-2021, retrieved from:
#https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html

data_directory = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/CMIP6_Hackathon/'

dataset = xr.open_dataset(data_directory+'icec.mnmean.nc')
x,y = np.meshgrid(dataset.lon.values, dataset.lat.values)

#%%

''' Plot ice extent change '''

month = 1

fig = plt.figure(figsize=(10,5))
gs1 = gridspec.GridSpec(1, 2)
gs1.update(wspace=0.04, hspace=0.05)

cmap = matplotlib.cm.get_cmap(cmocean.cm.phase)
lines = []
labels = []
i = 0

ax = plt.subplot(gs1[0], projection=ccrs.Orthographic(central_longitude=0, central_latitude=-90, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-180,180,-50,-90],crs=ccrs.PlateCarree()) 
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')

for i in range(4):
    yr = str(1981+i*10)
    yr2 = str(1990+i*10)
    data = dataset.sel(time=slice(yr+"-"+str("{:02d}".format(month)),yr2+"-"+str("{:02d}".format(month))))
    data_month = data.groupby('time.month').mean('time')
    rgba = cmap(i/4)
    contour = plt.contour(x, y, data_month.icec.sel(month=1), levels=[15], colors=[rgba], linewidths=0.8, transform=ccrs.PlateCarree())
    
    lines.extend(contour.collections)
    if yr==1981:
        labels.extend(['1982-1990'])
    else:
        labels.extend([yr+'-'+yr2])
    i += len(contour.collections)
ax.add_feature(land_50m, facecolor='#eeeeee')
ax.text(-135, -32,'a',weight='bold',fontsize=12,transform=ccrs.PlateCarree())

#%%

''' **Modelled sea ice edge** '''
### **Read modelled sea ice concentration**
# **Multi-model SSP5 8.5 (2015-2020)**
direc = '../DATA/SImon_ssp585/'
f = 'ETOPO_siconc_SImon_ACCESS-CM2_ssp585_201501-210012.nc'
acces = xr.open_dataset(direc+f).sel(time=slice('2015-01','2100-01'))
acces = acces.sel(time=acces.time.dt.month.isin([month]))
filenames = glob.glob(direc+'*.nc')

siconc = np.full((29,86,180,360),np.nan)
i=0
for f in filenames:
    cmip = xr.open_dataset(f).sel(time=slice('2015','2100'))
    cmip = cmip.sel(time=cmip.time.dt.month.isin([month]))
    try:
        siconc[i,:,:,:] = cmip.siconc
    except ValueError:
        siconc[i,:-1,:,:] = cmip.siconc
    i = i+1
    
model_mean = np.nanmean(siconc,axis=0)
model_mean[model_mean==0] = np.nan

year = np.zeros(len(cmip.time.values))
for i in range(len(cmip.time.values)):
    year[i] = cmip.time.values[i].year
    
x,y = np.meshgrid(cmip.ETOPO60X, cmip.ETOPO60Y)

#%%
# **Multi-model Historical**
direc = '../DATA/SImon_historical/'
f = 'ETOPO_siconc_SImon_ACCESS-CM2_historical_185001-201412.nc'
acces = xr.open_dataset(direc+f).sel(time=slice('1900','2014'))
acces = acces.sel(time=acces.time.dt.month.isin([month]))
filenames = glob.glob(direc+'*.nc')

siconc_hist = np.full((29,115,180,360),np.nan)
i=0
for f in filenames:
    cmip = xr.open_dataset(f).sel(time=slice('1900','2014'))
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
    
x,y = np.meshgrid(cmip.ETOPO60X, cmip.ETOPO60Y)  

#%%
# Plot CMIP6 multi-model mean sea ice extent

data = model_mean_hist[-4:,:,:]
data_proj = model_mean[:6,:,:]
data_comb = np.concatenate([data, data_proj])

cmap = matplotlib.cm.get_cmap(cmocean.cm.phase)
lines = []
labels = []

ax = plt.subplot(gs1[1], projection=ccrs.Orthographic(central_longitude=0, central_latitude=-90, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-180,180,-50,-90],crs=ccrs.PlateCarree()) 
ax.add_feature(land_50m, facecolor='#eeeeee')
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')

#Plot historical ice edge (1982-2010)
for i in range(3):
    yr = str(1981+i*10)
    yr2 = str(1990+i*10)
    print(yr+'-'+yr2)
    idx = np.where((year_hist>=int(yr))&(year_hist<=int(yr2)))
    rgba = cmap((i)/4)
    contour = plt.contour(x, y, np.nanmean(model_mean_hist[idx[0],:,:],axis=0), levels=[15], colors=[rgba], linewidths=0.8, transform=ccrs.PlateCarree(),zorder=0)
    lines.extend(contour.collections)
    labels.extend([yr+'-'+yr2])

#Plot 2011-2020
yr = str(2011)
yr2 = str(2020)
rgba = cmap((3)/4)
contour = plt.contour(x, y, np.nanmean(data_comb,axis=0), levels=[15], colors=[rgba], linewidths=0.8, transform=ccrs.PlateCarree(),zorder=0)
lines.extend(contour.collections)
labels.extend([yr+'-'+yr2])

ax.scatter(0,-90,s=7000,c='white', transform=ccrs.PlateCarree())
ax.text(-135, -32,'b',weight='bold',fontsize=12,transform=ccrs.PlateCarree())
plt.legend(lines, labels, bbox_to_anchor=(-0.02, -0.02), loc='upper center', fontsize=10, ncol=4,frameon=False)








