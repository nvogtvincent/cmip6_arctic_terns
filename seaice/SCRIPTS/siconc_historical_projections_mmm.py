# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 11:03:10 2022

@author: zq19140
"""

import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m')

# Set some plotting defaults
plt.rcParams['figure.dpi'] = 600

months = [11,1,3]
scenarios = ['ssp245', 'ssp585']
letter = ['a','b','c','d','e','f']

fig = plt.figure(figsize=(15,10))
gs1 = gridspec.GridSpec(2, 3)
gs1.update(wspace=0.04, hspace=0.05)

p=0

for scen in scenarios:
    for month in months:
        
        #%%
        direc = '../DATA/SImon_'+scen+'/'
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
        
        direc = '../DATA/SImon_historical/'
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
        
        cmap = matplotlib.cm.get_cmap('viridis')
        lines = []
        labels = []

        ax = plt.subplot(gs1[p], projection=ccrs.Orthographic(central_longitude=0, central_latitude=-90, globe=None))
        ax.coastlines(resolution='50m',linewidth=0.5)
        ax.set_extent([-180,180,-50,-90],crs=ccrs.PlateCarree()) 
        ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
        
        # #Plot historical ice edge (1901-2010)
        for i in range(11):
            yr = str(1901+i*10)
            yr2 = str(1900+10+i*10)
            print(yr+'-'+yr2)
            idx = np.where((year_hist>=int(yr))&(year_hist<=int(yr2)))
            rgba = cmap((i)/20)
            contour = plt.contour(x, y, np.nanmean(model_mean_hist[idx[0],:,:],axis=0), levels=[15], colors=[rgba], linestyles='--', linewidths=0.5, transform=ccrs.PlateCarree(),zorder=0)
            
            lines.extend(contour.collections)
            labels.extend([yr+'-'+yr2])
            
        
        #Plot 2011-2020
        yr = str(2011)
        yr2 = str(2020)
        idx = np.where((year>=int(yr))&(year<=int(yr2)))
        idx2 = np.where((year_hist>=int(yr))&(year_hist<=int(yr2)))
        data = np.nanmean(np.concatenate([model_mean[idx[0],:,:],model_mean_hist[idx2[0],:,:]]),axis=0)
        rgba = cmap((11)/20)
        contour = plt.contour(x, y, data, levels=[15], colors=[rgba], linestyles='-.', linewidths=0.5, transform=ccrs.PlateCarree(),zorder=0)
            
        lines.extend(contour.collections)
        labels.extend([yr+'-'+yr2])
        
        # Plot projections (2021-2100)    
        for i in range(8):
            yr = str(2021+i*10)
            yr2 = str(2020+10+i*10)
            print(yr+'-'+yr2)
            idx = np.where((year>=int(yr))&(year<=int(yr2)))
            rgba = cmap((i+12)/20)
            contour = plt.contour(x, y, np.nanmean(model_mean[idx[0],:,:],axis=0), levels=[15], colors=[rgba], linewidths=0.5, transform=ccrs.PlateCarree(),zorder=0)
            
            lines.extend(contour.collections)
            labels.extend([yr+'-'+yr2])
            
        ax.scatter(0,-90,s=7000,c='white', transform=ccrs.PlateCarree())
        ax.text(-135, -32,letter[p],weight='bold',fontsize=12,transform=ccrs.PlateCarree())
        ax.add_feature(land_50m, facecolor='#eeeeee')
    
        if p==2:
            plt.legend(lines, labels, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, frameon=False)
    
        p=p+1
    