#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compares flight time across CMIP6 ensemble
@author: Noam Vogt-Vincent
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import os
import pandas as pd
from netCDF4 import Dataset
from glob import glob
from scipy.ndimage import gaussian_filter1d
from scipy.stats import linregress
import cmasher as cmr

##############################################################################
# DIRECTORIES & PARAMETERS                                                   #
##############################################################################

# PARAMETERS
param = {'tern_per_yr': 6000,
         'fig_fn': 'tern_times.png',
         }

# DIRECTORIES
dirs = {'script': os.path.dirname(os.path.realpath(__file__)),
        'traj': os.path.dirname(os.path.realpath(__file__)) + '/TRAJ_DATA/TIME/',
        'figs': os.path.dirname(os.path.realpath(__file__)) + '/FIGURES/'}

# FILE HANDLES
fh = {'traj': sorted(glob(dirs['traj'] + '*.nc')),
      'fig': dirs['figs'] + param['fig_fn']}

##############################################################################
# SORT TRAJECTORY FILES                                                      #
##############################################################################

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

hist_years = np.arange(1850, 2015)
scen_years = np.arange(2015, 2101)

# Load data for all files
ensemble = {'HISTORICAL': {k: [] for k in hist_years},
            'SSP245': {k: [] for k in scen_years},
            'SSP585': {k: [] for k in scen_years},}

mean = {}
stdv = {}
time = {}
# time = {'ensemble': {'HISTORICAL': np.arange(1850, 2015),
#                      'SSP245': np.arange(2015, 2101),
#                      'SSP585': np.arange(2015, 2101)}}

global_max = 0
global_min = 1e9

for traj_file in fh['traj']:
    model_name = get_model_name(traj_file)
    scen_name = get_scen_name(traj_file)

    if model_name not in mean.keys():
        mean[model_name] = {}
        stdv[model_name] = {}
        time[model_name] = {}

    with Dataset(traj_file, mode='r') as nc:
        data_in = np.array(nc.variables['flight_time'][:])/(3600*24)
        data_in = np.ma.masked_where(data_in < 10, data_in)

        new_global_min = np.min([np.min(data_in), global_min])
        new_global_max = np.max([np.max(data_in), global_min])
        if not np.isnan(new_global_min):
            global_min = new_global_min

        if not np.isnan(new_global_max):
            global_max = new_global_max

        n_yr = len(data_in)/param['tern_per_yr']

        if scen_name == 'HISTORICAL':
            if n_yr % 100 != 65:
                raise ValueError('Unexcepted number of trajectories for ' + traj_file)
            else:
                yr0 = 2015-n_yr
                time[model_name][scen_name] = np.arange(yr0, yr0+n_yr, dtype=int)
        else:
            if n_yr != 86:
                raise ValueError('Unexcepted number of trajectories for ' + traj_file)
            else:
                time[model_name][scen_name] = np.arange(2015, 2015+n_yr, dtype=int)

        # Calculate the annual means
        data_in = data_in.reshape((-1, param['tern_per_yr']))
        mean[model_name][scen_name] = np.mean(data_in, axis=1)
        stdv[model_name][scen_name] = np.std(data_in, axis=1)

        for i, year in enumerate(time[model_name][scen_name]):
            ensemble[scen_name][year].append(data_in[i, :].data)

# Create histograms
time_axis = np.arange(15, 51)
hist = np.zeros((len(time_axis)-1, 165))
ssp245 = np.zeros((len(time_axis)-1, 86))
ssp585 = np.zeros((len(time_axis)-1, 86))

for i, year in enumerate(hist_years):
    data = np.array(ensemble['HISTORICAL'][year]).flatten()
    data = data[~np.isnan(data)]
    data = data[data > 10]
    ensemble['HISTORICAL'][year] = data
    hist[:, i] = np.histogram(data, bins=time_axis, density=True)[0]

for i, year in enumerate(scen_years):
    data = np.array(ensemble['SSP245'][year]).flatten()
    data = data[~np.isnan(data)]
    data = data[data > 10]
    ensemble['SSP245'][year] = data
    ssp245[:, i] = np.histogram(data, bins=time_axis, density=True)[0]

    data = np.array(ensemble['SSP585'][year]).flatten()
    data = data[~np.isnan(data)]
    data = data[data > 10]
    ensemble['SSP585'][year] = data
    ssp585[:, i] = np.histogram(data, bins=time_axis, density=True)[0]

# Calculate a decade-ly mean
decades = np.arange(1850, 2100, 10)
all_years = np.arange(1850, 2101, 1)
decade_mp = decades + 5
all_years_mp = all_years + 0.5

ssp245_decmean = np.zeros_like(decades, dtype=float)
ssp585_decmean = np.zeros_like(decades, dtype=float)

ssp245_annmean = np.zeros_like(all_years, dtype=float)
ssp585_annmean = np.zeros_like(all_years, dtype=float)

for i, decade in enumerate(decades):
    if decade < 2010:
        data = []
        for j in range(10):
            data.append(ensemble['HISTORICAL'][decade+j])
            ssp245_annmean[decade+j-1850] = np.mean(ensemble['HISTORICAL'][decade+j])
            ssp585_annmean[decade+j-1850] = np.mean(ensemble['HISTORICAL'][decade+j])

        data = np.array([item for sublist in data for item in sublist])

        ssp245_decmean[i] = np.mean(data)
        ssp585_decmean[i] = np.mean(data)
    elif decade > 2010:
        data1 = []
        data2 = []

        ypd = 11 if decade == 2090 else 10
        for j in range(ypd):
            data1.append(ensemble['SSP245'][decade+j])
            data2.append(ensemble['SSP585'][decade+j])
            ssp245_annmean[decade+j-1850] = np.mean(ensemble['SSP245'][decade+j])
            ssp585_annmean[decade+j-1850] = np.mean(ensemble['SSP585'][decade+j])

        data1 = np.array([item for sublist in data1 for item in sublist])
        data2 = np.array([item for sublist in data2 for item in sublist])

        ssp245_decmean[i] = np.mean(data1)
        ssp585_decmean[i] = np.mean(data2)
    elif decade == 2010:
        data1 = []
        data2 = []
        for j in np.arange(2010, 2015):
            data1.append(ensemble['HISTORICAL'][j])
            data2.append(ensemble['HISTORICAL'][j])
            ssp245_annmean[j-1850] = np.mean(ensemble['HISTORICAL'][j])
            ssp585_annmean[j-1850] = np.mean(ensemble['HISTORICAL'][j])
        for j in np.arange(2015, 2020):
            data1.append(ensemble['SSP245'][j])
            data2.append(ensemble['SSP585'][j])
            ssp245_annmean[j-1850] = np.mean(ensemble['SSP245'][j])
            ssp585_annmean[j-1850] = np.mean(ensemble['SSP585'][j])

        data1 = np.array([item for sublist in data1 for item in sublist])
        data2 = np.array([item for sublist in data2 for item in sublist])

        ssp245_decmean[i] = np.mean(data1)
        ssp585_decmean[i] = np.mean(data2)


f, ax = plt.subplots(nrows=2, ncols=1, figsize=(27, 15))
f.subplots_adjust(hspace=0.12)

hist1a = ax[0].pcolormesh(np.arange(1850, 1951)-0.5, time_axis, hist[:, :100],
                          cmap=cmr.neutral_r, vmin=0, vmax=np.max(hist), alpha=0.5)
hist1b = ax[1].pcolormesh(np.arange(1850, 1951)-0.5, time_axis, hist[:, :100],
                          cmap=cmr.neutral_r, vmin=0, vmax=np.max(hist), alpha=0.5)

hist2a = ax[0].pcolormesh(np.arange(1950, 2016)-0.5, time_axis, hist[:, 100:],
                          cmap=cmr.neutral_r, vmin=0, vmax=np.max(hist))
hist2b = ax[1].pcolormesh(np.arange(1950, 2016)-0.5, time_axis, hist[:, 100:],
                          cmap=cmr.neutral_r, vmin=0, vmax=np.max(hist))

scen_a = ax[0].pcolormesh(np.arange(2015, 2102)-0.5, time_axis, ssp245,
                          cmap=cmr.arctic_r, vmin=0, vmax=np.max(hist))
scen_b = ax[1].pcolormesh(np.arange(2015, 2102)-0.5, time_axis, ssp585,
                          cmap=cmr.flamingo_r, vmin=0, vmax=np.max(hist))

for i in range(2):
    ax[i].plot([1950, 1950], [15, 45], 'k--', linewidth=1)
    ax[i].plot([2015, 2015], [15, 45], 'k--', linewidth=1)
    ax[i].set_ylim([15, 45])
    ax[i].set_xlim([1850, 2100])
    # ax[i].set_aspect('equal')
    ax[i].set_xticks([1850, 1900, 1950, 2000, 2015, 2050, 2100])
    ax[i].tick_params(axis='x', labelsize=22)
    ax[i].tick_params(axis='y', labelsize=22)
    ax[i].text(1852, 41, ['a', 'b'][i], fontsize=24, fontweight='bold')
    ax[i].text(1852, 16, 'Historical (5 models)', fontsize=24)
    ax[i].text(1952, 16, 'Historical (10 models)', fontsize=24)
    ax[i].text(2017, 16, ['SSP245 (10 models)' , 'SSP585 (10 models)'][i], fontsize=24)
    ax[i].set_ylabel('Migration time (days)', fontsize=24)

ax[0].plot(decade_mp, ssp245_decmean, 'k--', linewidth=1, marker='s', label='Decadal mean')
ax[1].plot(decade_mp, ssp585_decmean, 'k--', linewidth=1, marker='s')
ax[0].legend(frameon=False, fontsize=24)

ax[1].set_xlabel('Year', fontsize=24)

plt.savefig(fh['fig'], dpi=300)

# # Do linear regression
ssp245a_linreg = linregress(all_years_mp[165:], ssp245_annmean[165:])
ssp245d_linreg = linregress(decade_mp[16:], ssp245_decmean[16:])

ssp585a_linreg = linregress(all_years_mp[165:], ssp585_annmean[165:])
ssp585d_linreg = linregress(decade_mp[16:], ssp585_decmean[16:])

hista_linreg = linregress(all_years_mp[:160], ssp585_annmean[:160])
histb_linreg = linregress(decade_mp[:15], ssp585_decmean[:15])

hista_short_linreg = linregress(all_years_mp[100:160], ssp585_annmean[100:160])
histb_short_linreg = linregress(decade_mp[10:15], ssp585_decmean[10:15])

print('SSP245')
print('Annual mean trend:' + str(ssp245a_linreg[0]) + ' d yr-1 (p=' + str(ssp245a_linreg[3]) + ')')
print('Decadal mean trend:' + str(ssp245d_linreg[0]) + ' d yr-1 (p=' + str(ssp245d_linreg[3]) + ')')
print('')
print('SSP585')
print('Annual mean trend:' + str(ssp585a_linreg[0]) + ' d yr-1 (p=' + str(ssp585a_linreg[3]) + ')')
print('Decadal mean trend:' + str(ssp585d_linreg[0]) + ' d yr-1 (p=' + str(ssp585d_linreg[3]) + ')')
print('')
print('Historical 1850-2010')
print('Annual mean trend:' + str(hista_linreg[0]) + ' d yr-1 (p=' + str(hista_linreg[3]) + ')')
print('Decadal mean trend:' + str(histb_linreg[0]) + ' d yr-1 (p=' + str(histb_linreg[3]) + ')')
print('')
print('Historical 1950-2010')
print('Annual mean trend:' + str(hista_short_linreg[0]) + ' d yr-1 (p=' + str(hista_short_linreg[3]) + ')')
print('Decadal mean trend:' + str(histb_short_linreg[0]) + ' d yr-1 (p=' + str(histb_short_linreg[3]) + ')')

# Get bulk statistics
bulk_historical = np.concatenate([ensemble['HISTORICAL'][x] for x in np.arange(1950,2000)])
bulk_SSP245 = np.concatenate([ensemble['SSP245'][x] for x in np.arange(2050,2100)])
bulk_SSP585 = np.concatenate([ensemble['SSP585'][x] for x in np.arange(2050,2100)])

print('')
print('Mean migration times:')
print('Historical: ' + str(np.mean(bulk_historical)) +  '±' + str(np.std(bulk_historical)) + ' days')
print('SSP245: ' + str(np.mean(bulk_SSP245)) +  '±' + str(np.std(bulk_SSP245)) + ' days')
print('SSP585: ' + str(np.mean(bulk_SSP585)) +  '±' + str(np.std(bulk_SSP585)) + ' days')

