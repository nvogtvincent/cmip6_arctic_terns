# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 14:30:25 2021

Purpose
-------
    Plot the tern area-use density 
    


@author: pearseb
"""

#%% imports
 
import os
import numpy as np
import netCDF4 as nc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import matplotlib.cm as cm
import seaborn as sb
sb.set(style='ticks')


#%% data

#os.chdir('C://Users/pearseb/Dropbox/PostDoc/CMIP6 hackathon/data')
os.chdir('/Users/pbuchanan/Dropbox/PostDoc/collaborations/CMIP6 hackathon 2021/data')

data = nc.Dataset('bird_heatmap.nc','r')
density = data.variables['density'][...]
lon = data.variables['longitude'][...]
lat = data.variables['latitude'][...]

lon = np.ma.concatenate((lon[:,:], lon[:,0, np.newaxis]+360.0), axis=1)
lat = np.ma.concatenate((lat[:,:], lat[:,0, np.newaxis]), axis=1)
density = np.ma.concatenate((density[:,:], density[:,0, np.newaxis]), axis=1)


tracks = pd.read_csv('Collated migration data.csv')
bird_ids = tracks['Bird ID'].unique()


#%% plot

proj = ccrs.Robinson(central_longitude=20)

levs1 = np.arange(0,21,1)*0.1
colmap1 = cm.viridis

fstic = 13
fslab = 15
lw = 1
ls = '-'
alf = 0.75
col = 'firebrick'


fig = plt.figure(figsize=(7,9))
gs = GridSpec(2,1)

ax1 = plt.subplot(gs[0], projection=proj)
gl1 = ax1.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=0.75, zorder=3, draw_labels=True, xlocs=np.arange(-180,181,60), ylocs=np.arange(-60,61,20))
gl1.rotate_labels = False
for iid in bird_ids[0:17]:
    if iid != bird_ids[17]:
        plt.plot(tracks[tracks['Bird ID'] == iid]['Long'], tracks[tracks['Bird ID'] == iid]['Lat'], transform=ccrs.PlateCarree(), linestyle=ls, linewidth=lw, alpha=alf, color=col, zorder=3)
ax1.add_feature(cfeature.LAND, color='w', zorder=2)
ax1.coastlines(zorder=2)
ax1.set_global()

ax2 = plt.subplot(gs[1], projection=proj)
gl2 = ax2.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=0.75, zorder=3, draw_labels=True, xlocs=np.arange(-180,181,60), ylocs=np.arange(-60,61,20))
gl2.rotate_labels = False
#gl2.left_labels = False
p2 = plt.contourf(lon, lat, density, transform=ccrs.PlateCarree(), cmap=colmap1, levels=levs1, vmin=np.min(levs1), vmax=np.max(levs1), extend='max')
ax2.add_feature(cfeature.LAND, color='w', zorder=2)
ax2.coastlines(zorder=2)


plt.subplots_adjust(right=0.95, top=0.95, bottom=0.175, left=0.05)


cbax1 = fig.add_axes([0.2, 0.1, 0.6, 0.03])
cbar1 = plt.colorbar(p2, cax=cbax1, orientation='horizontal', ticks=levs1[::2])
cbar1.ax.set_xlabel('Total time spent at location by terns (days)', fontsize=fslab)
cbar1.ax.tick_params(labelsize=fstic)

xx = 0.05; yy = 0.95
plt.text(xx,yy, 'a', transform=ax1.transAxes, va='center', ha='center', fontsize=fslab+2, fontweight='bold')
plt.text(xx,yy, 'b', transform=ax2.transAxes, va='center', ha='center', fontsize=fslab+2, fontweight='bold')


#%% save

#os.chdir('C://Users/pearseb/Dropbox/PostDoc/CMIP6 hackathon/figures')
os.chdir('/Users/pbuchanan/Dropbox/PostDoc/collaborations/CMIP6 hackathon 2021/figures')
fig.savefig('fig1_ternDensity.png', dpi=300, bbox_inches='tight')