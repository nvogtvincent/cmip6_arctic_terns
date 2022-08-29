# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 14:07:33 2021

@author: Isolde Glissenaar
"""

''' **Observed and modelled sea ice extent** '''

#Plots the change in decadal mean sea ice extent (1981-2020) from observations and UKESM.
#Observations (icec.mnmean.nc) can be retrieved from https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html (file to big to post on GitHub).

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import calendar
import cartopy.crs as ccrs
import cartopy.feature as cfeature
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m')
from matplotlib import cm
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

data_directory = '../DATA/'
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
# **UKESM Historical**
hist = xr.load_dataset(data_directory+'ETOPO_siconc_SImon_UKESM1-0-LL_historical_r1i1p1f2_198201-201411_SH.nc').sel(time=slice('1982','2014-11'))
x,y = np.meshgrid(hist.ETOPO60X, hist.ETOPO60Y)

project = xr.load_dataset(data_directory+'ETOPO_siconc_SImon_UKESM1-0-LL_ssp585_r1i1p1f2_201501-202101_SH.nc').sel(time=slice('2015','2021-01'))
x,y = np.meshgrid(project.ETOPO60X, project.ETOPO60Y)

#%%
data = hist.sel(time=slice("2011-01","2014-01"))
data_proj = project.sel(time=slice("2015-01","2020-01"))
data_comb = xr.concat([data, data_proj], "time")

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
    data = hist.sel(time=slice(yr+"-01",yr2+"-01"))
    data_month = data.groupby('time.month').mean('time')
    rgba = cmap((i)/4)
    contour = plt.contour(x, y, data_month.siconc.sel(month=1), levels=[0.15], colors=[rgba], linewidths=0.8, transform=ccrs.PlateCarree())
    lines.extend(contour.collections)
    labels.extend([yr+'-'+yr2])

#Plot 2011-2020
yr = str(2011)
yr2 = str(2020)
data_comb["siconc"].groupby('time.month').mean('time')
rgba = cmap((3)/4)
contour = plt.contour(x, y, data_comb.siconc[0,:,:], levels=[0.15], colors=[rgba], linewidths=0.8, transform=ccrs.PlateCarree())
lines.extend(contour.collections)
labels.extend([yr+'-'+yr2])

ax.text(-135, -32,'b',weight='bold',fontsize=12,transform=ccrs.PlateCarree())

plt.legend(lines, labels, bbox_to_anchor=(-0.02, -0.02), loc='upper center', fontsize=10, ncol=4,frameon=False)
plt.show()

