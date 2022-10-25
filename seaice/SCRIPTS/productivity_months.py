# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 13:13:39 2022

@author: zq19140
"""

import pandas as pd
import xarray as xr
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m')
ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m')
import matplotlib
import matplotlib.path as mpath
import matplotlib.gridspec as gridspec

# Set some plotting defaults
plt.rcParams['figure.figsize'] = (6, 4)
plt.rcParams['figure.dpi'] = 100

#%%
''' **Read data** '''
# **Read observed (OISST) sea ice concentration data**
data_directory = '../DATA/'
dataset = xr.open_dataset(data_directory+'icec.mnmean.nc')
x,y = np.meshgrid(dataset.lon.values, dataset.lat.values)


# Productivity 
#To large a file to upload to GitHub
prod = xr.open_dataset('../DATA/global-reanalysis-bio-intnpp_011993-112019.nc')
prod_lon, prod_lat = np.meshgrid(prod.LONGITUDE.values, prod.LATITUDE.values)
intnpp = prod.INTNPP.values
time_prod = pd.DatetimeIndex(prod.TIME.values)

#Make circle for plot
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)


#%%

months = [10,11,12,1,2,3]

fig = plt.figure(figsize=(16,10),dpi=600)
projection=ccrs.Orthographic(central_longitude=0, central_latitude=-90, globe=None)

gs1 = gridspec.GridSpec(2, 3)
gs1.update(wspace=0.04, hspace=0.05)

for j in range(len(months)):
    month = months[j]
    if month==12:
        yrs = np.arange(1993,2019,1)
    else:
        yrs = np.arange(1993,2020,1)
    
    idx_prod_yr = np.where((time_prod.month==month))[0]
    intnpp_m = intnpp[idx_prod_yr,:,:]
    
    trend_prod = np.zeros((prod_lat.shape)); trend_prod[:,:] = np.nan
    
    for i in range(prod_lat.shape[0]):
        for k in range(prod_lat.shape[1]):
            trend_prod[i,k], intercept, r_value, p_value, std_err = stats.linregress(yrs, intnpp_m[:,i,k])
    
    #%%
    ax = plt.subplot(gs1[j], projection=ccrs.Orthographic(central_longitude=0, central_latitude=-90, globe=None))
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.coastlines(resolution='50m',linewidth=0.5, zorder=6)
    ax.set_extent([-180,180,-50,-90],crs=ccrs.PlateCarree()) 
    ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':', zorder=15)
    pc=ax.pcolormesh(prod.LONGITUDE, prod.LATITUDE[prod.LATITUDE<-50], trend_prod[prod.LATITUDE<-50,:], 
                       cmap='RdBu_r', vmin=-20, vmax=20, clip_path=(circle, ax.transAxes), transform=ccrs.PlateCarree(),zorder=2)

    ax.add_feature(land_50m, facecolor='#eeeeee', zorder=7) 
    ax.add_feature(ocean_50m, facecolor='lightgrey', zorder=1)    
    
    #Plot bird locations
    data_directory = '../DATA/'
    migration = pd.read_csv(data_directory+'Collated migration data.csv')
    migr_temp = migration[pd.to_datetime(migration['Date'], dayfirst=True).dt.month==month]
    migr_temp = migr_temp.reset_index(drop=True)
    migr_temp = migr_temp.drop(np.where(migr_temp['Lat']>-50)[0]).reset_index(drop=True)
    scat = ax.scatter(migr_temp['Long'], migr_temp['Lat'], c='black', s=5, marker='.', 
                            edgecolors='none', transform=ccrs.PlateCarree(), zorder=15)
    
    
    #Plot observed ice edge
    dataset = xr.open_dataset(data_directory+'icec.mnmean.nc')
    x,y = np.meshgrid(dataset.lon.values, dataset.lat.values)
    cmap = matplotlib.cm.get_cmap('viridis')
    ls = [':', '-.', '--', '-']
    lines = []
    labels = []
    for i in range(0,4):
        yr = str(1981+i*10)
        yr2 = str(1990+i*10)
        data = dataset.sel(time=slice(yr+"-"+str("{:02d}".format(month)),yr2+"-"+str("{:02d}".format(month))))
        data_month = data.groupby('time.month').mean('time')
        rgba = cmap(i/3)
        contour = plt.contour(x, y, data_month.icec.sel(month=month), levels=[15], linestyles=ls[i], colors='#14165c', linewidths=0.5, transform=ccrs.PlateCarree())
        if yr=='1981':
            labels.extend(['1982-1990'])
            lines.extend(contour.collections)
        else:
            labels.extend([yr+'-'+yr2])
            lines.extend(contour.collections)
            
    ax.set_boundary(circle, transform=ax.transAxes)
    
    

plt.legend(lines, labels, bbox_to_anchor=(-0.5, -0.15), loc='center',frameon=False, ncol=4, fontsize=18)   

fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8,
                    wspace=0.02, hspace=0.02)
cb_ax = fig.add_axes([0.83, 0.1, 0.015, 0.8])
cbar = fig.colorbar(pc, cax=cb_ax, extend='both', extendfrac=0.02)
cbar.ax.tick_params(labelsize=15)
cbar.set_label(label='NPP [g C m$^{-2}$ day$^{-1}$] / yr', fontsize=15)

fig.text(0.1,0.53,'a', fontsize=20, weight='bold')
fig.text(0.34,0.53,'b', fontsize=20, weight='bold')
fig.text(0.58,0.53,'c', fontsize=20, weight='bold')
fig.text(0.1,0.12,'d', fontsize=20, weight='bold')
fig.text(0.34,0.12,'e', fontsize=20, weight='bold')
fig.text(0.58,0.12,'f', fontsize=20, weight='bold')

fig.tight_layout()







        