# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 10:51:44 2021

@author: zq19140
"""

'''Sea ice edge and bird location'''

# Plots the sea ice edge observed OISST and locations of tracked Arctic Terns

import pandas as pd
import xarray as xr
import numpy as np
import calendar
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m')
ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m')
import matplotlib
import cmocean

# Set some plotting defaults
plt.rcParams['figure.figsize'] = (6, 4)
plt.rcParams['figure.dpi'] = 200

#%%
''' **Read data** '''
# **Read observed (OISST) sea ice concentration data**

data_directory = '../DATA/'
dataset = xr.open_dataset(data_directory+'icec.mnmean.nc')
x,y = np.meshgrid(dataset.lon.values, dataset.lat.values)

# Read bird tracking data
data_directory = '../DATA/'
migration = pd.read_csv(data_directory+'Collated migration data.csv')


#%%

# Plot sea ice edge and bird locations 

# Bird IDs with data all boreal summer
# 'A (Y1)', 'A (Y2)', 'ARTE_370', 'ARTE_371'
# 'ARTE_373', 'ARTE_376', 'ARTE_390', 'ARTE_395', 'ARTE_406'
# 'ARTE_408', 'ARTE_410', 'B', 'C', 'D', 'E (Y1)', 'E (Y2)', 'F'
# 'G (Y1)'

# 2007: ARTE_370, ARTE_371, ARTE_373, ARTE_376, 
#       ARTE_390, ARTE_395, ARTE_406, ARTE_408, 
#       ARTE_410
# 2008: A (Y1), B
# 2009: A (Y2), C, D
# 2011: E (Y1), F
# 2012: E (Y2)
# 2013: G (Y1)

bird = 'ARTE_370'
birds = ['ARTE_370', 'ARTE_371', 'ARTE_373', 'ARTE_376', 'ARTE_390', 
         'ARTE_395', 'ARTE_406', 'ARTE_408', 'ARTE_410']
markers = ['o', 'v', '>', '<', '^', 'd', 's', 'h', 'p']

migration_1 = migration[migration['Bird ID']==bird]
year = pd.to_datetime(migration_1['Date']).dt.year.iloc[0]
cmap = matplotlib.cm.get_cmap(cmocean.cm.phase)

fig=plt.figure(dpi=600)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-10, central_latitude=-60, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-50,40,-50,-80],crs=ccrs.PlateCarree()) 
ax.add_feature(land_50m, facecolor='#eeeeee') 
# ax.add_feature(ocean_50m, facecolor='lightgrey')    
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')

months = [11,12,1,2,3]
years = [year,year,year+1,year+1,year+1]
lines = []
labels = []

markers_leg = []
j=0
for i in range(5):
    for j in range(len(birds)):
        migration_1 = migration[migration['Bird ID']==birds[j]]
        migr_temp = migration_1[pd.to_datetime(migration_1['Date'], dayfirst=True).dt.month==months[i]]
        scat = plt.scatter(migr_temp['Long'], migr_temp['Lat'], c=[cmap(i/5)], s=5, marker=markers[j], edgecolors='none', transform=ccrs.PlateCarree())
        scat2 = plt.scatter(np.nan, np.nan, c='black', s=5, marker=markers[j], edgecolors='none', transform=ccrs.PlateCarree())
        if i==0:
            markers_leg.append(scat2)
            
    contour = plt.contour(x, y, dataset.icec.sel(time=str(years[i])+'-'+str("{:02d}".format(months[i])))[0,:,:], colors=[cmap(i/5)], levels=[15], linewidths=1, transform=ccrs.PlateCarree())
    lines.extend(contour.collections)
    labels.extend([calendar.month_name[months[i]]])
    j += len(contour.collections)
  
ax.text(0.02,0.94,'a', fontsize=11, weight='bold',transform=ax.transAxes)

legend1 = plt.legend(lines, labels, bbox_to_anchor=(1.02, 1), loc='upper left', title='Sea ice edge', title_fontsize=6, fontsize=6, frameon=False)
legend1._legend_box.align = "left"
legend2 = plt.legend(markers_leg, birds, bbox_to_anchor=(1.02,0.7), loc='upper left', title='Bird ID', title_fontsize=6, fontsize=6, frameon=False)
legend2._legend_box.align = "left"
plt.gca().add_artist(legend1)

