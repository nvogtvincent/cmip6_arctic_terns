# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:55:42 2021

@author: zq19140
"""

''' **Observed sea ice extend, decadal trend**'''
#Plots the observed sea ice extend change (1981-2020) for Nov-Mar. 

import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m')
import matplotlib.path as mpath
import matplotlib.gridspec as gridspec
import cmocean

# Set some plotting defaults
plt.rcParams['figure.figsize'] = (6, 4)
plt.rcParams['figure.dpi'] = 200

#%%
''' **Read ice concentration data** '''

#NOAA OISST ice concentration 1981-2021, retrieved from:
#https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html

data_directory = '../DATA/'
dataset = xr.open_dataset(data_directory+'icec.mnmean.nc')
x,y = np.meshgrid(dataset.lon.values, dataset.lat.values)

#%%

# Plot sea ice edge Nov-Mar

cmap = matplotlib.cm.get_cmap(cmocean.cm.phase) 

#Make circle for plot
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)


months = [10,11,12,1,2,3]
letter = ['a','b','c','d','e','f']

fig = plt.figure(figsize=(15,10), dpi=600)
projection=ccrs.Orthographic(central_longitude=0, central_latitude=-90, globe=None)

gs1 = gridspec.GridSpec(2, 3)
gs1.update(wspace=0, hspace=0.05)

for j in range(len(months)):
    lines = []
    labels = []
    i = 0
    
    ax = plt.subplot(gs1[j], projection=ccrs.Orthographic(central_longitude=0, central_latitude=-90, globe=None))
    ax.coastlines(resolution='50m',linewidth=0.5)
    ax.set_extent([-180,180,-50,-90],crs=ccrs.PlateCarree()) 
    ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')

    for i in range(4):
        yr = str(1981+i*10)
        yr2 = str(1990+i*10)
        data = dataset.sel(time=slice(yr+"-"+str("{:02d}".format(months[j])),yr2+"-"+str("{:02d}".format(months[j]))))
        data_month = data.groupby('time.month').mean('time')
        rgba = cmap(i/4)
        contour = plt.contour(x, y, data_month.icec.sel(month=months[j]), levels=[15], colors=[rgba], linewidths=0.8, transform=ccrs.PlateCarree())
    
        lines.extend(contour.collections)
        if yr==1981:
            labels.extend(['1982-1990'])
        else:
            labels.extend([yr+'-'+yr2])
        i += len(contour.collections)
    
    ax.add_feature(land_50m, facecolor='#eeeeee')
    ax.text(-135, -32, letter[j], size=16, weight='bold', transform=ccrs.PlateCarree())
    
    ax.set_boundary(circle, transform=ax.transAxes)
    
plt.legend(lines, labels, bbox_to_anchor=(-0.55, -0.1), loc='upper center', fontsize=14, ncol=4, frameon=False)



