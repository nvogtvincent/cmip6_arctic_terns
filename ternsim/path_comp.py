#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compares flight routes across CMIP6 ensemble to observations
@author: Noam Vogt-Vincent
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import os
from netCDF4 import Dataset
from glob import glob
import cmasher as cmr
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

##############################################################################
# DIRECTORIES & PARAMETERS                                                   #
##############################################################################

# PARAMETERS
param = {'terns_per_release': 2000,
         'fig_fn': 'tern_routes_comparison.png',

         # Plotting variables
         'res'           : 1, # (degrees)
         'bounds'        : {'lon_min' : -90,
                            'lon_max' : +50,
                            'lat_min' : -75,
                            'lat_max' : +70},
         'fig_fn_1'        : 'histogram.png',
         'fig_fn_2'        : 'limits.png',

         'source_bnds'   : {'lon_min' : -60,
                            'lon_max' : -40,
                            'lat_min' : -70,
                            'lat_max' : -70},
         'sink_lat'      : 60.,
         'target_coords' : {'lon'     : -20,
                            'lat'     : 65},
         'hours_dt'      : 2,
         'cmap1'         : cmr.flamingo_r,
         'cmap2'         : cmr.fusion_r
         }


# DIRECTORIES
dirs = {'script': os.path.dirname(os.path.realpath(__file__)),
        'traj': os.path.dirname(os.path.realpath(__file__)) + '/TRAJ_DATA/TRAJ/',
        'obs': os.path.dirname(os.path.realpath(__file__)) + '/TRAJ_DATA/OBS/',
        'figs': os.path.dirname(os.path.realpath(__file__)) + '/FIGURES/',
        'grid': os.path.dirname(os.path.realpath(__file__)) + '/GRID_DATA/'}

# FILE HANDLES
fh = {'traj': sorted(glob(dirs['traj'] + '*.nc')),
      'obs': dirs['obs'] + 'collated_migration_data.xlsx',
      'fig': dirs['figs'] + param['fig_fn'],
      'lsm': dirs['grid'] + 'GSHHS_L/GSHHS_l_L1.shp'}

# Make a list of model names
model_names = []

get_model_name = lambda file_name : file_name.split('/')[-1].split('_')[0]
get_scen_name = lambda file_name : file_name.split('/')[-1].split('_')[1]

for i in range(len(fh['traj'])):
    model_name = get_model_name(fh['traj'][i])
    if model_name not in model_names:
        model_names.append(model_name)

# Check all files are present
for i in range(len(model_names)):
    for scenario in ['HISTORICAL', 'SSP245', 'SSP585']:
        check_name = dirs['traj'] + model_names[i] + '_' + scenario + '_TRAJ.nc'
        if check_name not in fh['traj']:
            raise FileNotFoundError(check_name + ' does not exist!')

print('Model names:')
print(model_names)

##############################################################################
# PROCESS MODEL DATA                                                           #
##############################################################################

# Define years for averages
hist_years = np.arange(1950, 2000)
scen_years = np.arange(2050, 2100)

# Create grids for binning
x_bnd = np.arange(param['bounds']['lon_min'],
                  param['bounds']['lon_max']+param['res'],
                  param['res'])

y_bnd = np.arange(param['bounds']['lat_min'],
                  param['bounds']['lat_max']+param['res'],
                  param['res'])

x_bnd_mp = 0.5*(x_bnd[1:] + x_bnd[:-1])
y_bnd_mp = 0.5*(y_bnd[1:] + y_bnd[:-1])

grid = np.zeros((len(x_bnd_mp), len(y_bnd_mp)), dtype=np.float64) # HIST/245/585
n_traj_tot = np.array([0, 0, 0])

gshhs = None
lsm = None

# Grid the model output
for model_i, model in enumerate(model_names):
    scenario_i = 0
    scenario = 'HISTORICAL'

    output_name = dirs['traj'] + model + '_' + scenario + '_TRAJ.nc'
    print('Gridding ' + scenario + ' ' + model + '...')

    with Dataset(output_name, mode='r') as nc:
        t0 = nc.variables['time'][:, 0]
        nyr = int(len(t0)/(param['terns_per_release']*3))

        if scenario_i == 0:
            release_year = 2015-nyr

        else:
            if nyr != 86:
                raise ValueError('Unexpected number of years in scenario data!')
            else:
                release_year = 2101-nyr

        year_list = np.arange(release_year, release_year+nyr)
        year_list = np.repeat(year_list, param['terns_per_release']*3)

        assert len(year_list) == len(t0)

        if scenario_i == 0:
            traj_idx = np.where(np.isin(year_list, hist_years))[0]
        else:
            traj_idx = np.where(np.isin(year_list, scen_years))[0]

        # Load and grid data
        lon = nc.variables['lon'][traj_idx, :]
        lat = nc.variables['lat'][traj_idx, :]
        n_traj = np.shape(lon)[0]+1 # Number of unique terns gathered
        n_traj_tot[scenario_i] += n_traj

        tid = np.zeros_like(lon) # Trajectory identifier
        tid[:] = np.arange(np.shape(lon)[0]).reshape(-1, 1)+1
        tid = np.ma.masked_array(tid, mask=lon.mask).compressed()
        lon = lon.compressed()
        lat = lat.compressed()

        # Split arrays if necessary
        n_splits = n_traj*len(x_bnd)*len(y_bnd)/5e8
        n_splits = 1 if n_splits < 1 else int(np.ceil(n_splits))

        lon = np.array_split(lon, n_splits)
        lat = np.array_split(lat, n_splits)
        tid = np.array_split(tid, n_splits)

        # Now calculate the number of unique trajectories passing through each cell, and normalise
        for i in range(len(lon)):
            grid_ = np.histogramdd(np.array([lon[i], lat[i], tid[i]]).T, bins=[x_bnd, y_bnd, np.arange(tid[i].min(), tid[i].max()+2, 1)-0.5])[0]

            grid_[grid_ > 1] = 1
            grid[:, :] += np.sum(grid_, axis=2)

grid /= n_traj_tot[0]
# Now calculate differences
# grid[0] = HIST, grid[1] = SSP245, grid[2] = SSP585, grid[3] = SSP245-HIST, grid[4] = SSP585-HIST

##############################################################################
# PROCESS OBS DATA                                                           #
##############################################################################

obs_data = pd.read_excel(fh['obs'])
tern_list = pd.unique(obs_data['Bird ID'])

# Filter by date
obs_data['month'] = pd.to_datetime(obs_data['Date']).dt.month
obs_data = obs_data.loc[(obs_data['month'] >= 3)*(obs_data['month'] <= 6)]

# Replace bird IDs
for i in range(len(tern_list)):
    obs_data = obs_data.replace(tern_list[i], i)

tern_list = np.unique(obs_data['Bird ID'])
bid_bnd = np.append(tern_list, tern_list[-1]+1) - 0.5

# Create grids for binning
x_bnd_b = np.arange(param['bounds']['lon_min'],
                    param['bounds']['lon_max']+param['res'],
                    5)

y_bnd_b = np.arange(param['bounds']['lat_min'],
                    param['bounds']['lat_max']+param['res'],
                    5)

# Grid
b_grid = np.histogramdd(np.array([obs_data['Long'], obs_data['Lat'], obs_data['Bird ID']]).T,
                        bins=[x_bnd_b, y_bnd_b, bid_bnd])[0]
b_grid[b_grid > 1] = 1
b_grid = np.mean(b_grid, axis=2)
b_grid = np.ma.masked_where(b_grid == 0, b_grid)


##############################################################################
# PLOT DATA                                                                  #
##############################################################################

##############################################################################
# HISTORICAL TERN LOCATIONS ##################################################
##############################################################################

f = plt.figure(constrained_layout=True, figsize=(25, 20))
gs = GridSpec(1, 3, figure=f, width_ratios=[1, 1, 0.07])
ax = []
ax.append(f.add_subplot(gs[0, 0], projection = ccrs.Robinson(central_longitude=-30)))
ax.append(f.add_subplot(gs[0, 1], projection = ccrs.Robinson(central_longitude=-30)))
ax.append(f.add_subplot(gs[0, 2]))
label_list = ['a', 'b', 'c', 'd', 'e']

data_crs = ccrs.PlateCarree()

f.subplots_adjust(hspace=0.08, wspace=0.04)

gl = []
hist = []

land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor='gray',
                                        zorder=1)

for i in range(2):
    ax[i].set_aspect(1)
    ax[i].set_extent([-85, 25, -90, 90], crs=data_crs)

    # Add start/finish line
    if i == 0:
        start_line = ax[i].plot(np.array([param['source_bnds']['lon_min'],
                                          param['source_bnds']['lon_max']]),
                                np.array([param['source_bnds']['lat_min'],
                                          param['source_bnds']['lat_max']]),
                                'k-', linewidth=2, transform=data_crs,
                                zorder=100)

        start_bnds = ax[i].scatter(np.array([param['source_bnds']['lon_min'],
                                              param['source_bnds']['lon_max']]),
                                    np.array([param['source_bnds']['lat_min'],
                                              param['source_bnds']['lat_max']]),
                                    c='k', marker='.', s = 100,
                                    zorder=100,
                                    transform=data_crs)

        target_pnt = ax[i].scatter(param['target_coords']['lon'],
                                    param['target_coords']['lat'],
                                    c='k', marker='+', s=100,
                                    transform=data_crs,
                                    zorder=100)

    # Add cartographic features
    gl.append(ax[i].gridlines(crs=data_crs, draw_labels=True,
                              linewidth=0.5, color='black', linestyle='--', zorder=11))
    gl[i].xlocator = mticker.FixedLocator(np.arange(-210, 210, 60))
    gl[i].ylocator = mticker.FixedLocator(np.arange(-90, 120, 30))
    gl[i].xlabels_top = False
    gl[i].ylabels_right = False
    gl[i].ylabels_left = False if i > 0 else True
    gl[i].ylabel_style = {'size': 24}
    gl[i].xlabel_style = {'size': 24}

    ax[i].add_feature(land_10m)
    ax[i].coastlines(resolution='10m', zorder=12)

    # Plot the colormesh
    cmap = param['cmap1']
    if i == 0:
        hist.append(ax[i].contourf(x_bnd_mp, y_bnd_mp, 100*grid.T,
                                    cmap=cmap, levels=np.linspace(1, 7, num=13), transform=data_crs,
                                    zorder=2, extend='max'))
        ax[i].title.set_text('vTerns')
        ax[i].title.set_size(28)
    else:
        hist.append(ax[i].pcolormesh(x_bnd_b, y_bnd_b, 100*b_grid.T, cmap=cmap, transform=data_crs,
                                     zorder=2))
        ax[i].title.set_text('Real terns')
        ax[i].title.set_size(28)

    xpos = -122 if i == 0 else -120

    ax[i].text(xpos, -85, label_list[i], transform=data_crs, fontsize=30, va='bottom', ha='left', fontweight='bold')

cb = plt.colorbar(hist[0], cax=ax[2])
cb.set_label('Percentage of (v)Terns passing through cell on northbound migration', size=26)
ax[2].tick_params(axis='y', labelsize=24)

plt.savefig(fh['fig'], dpi=300, bbox_inches='tight')


