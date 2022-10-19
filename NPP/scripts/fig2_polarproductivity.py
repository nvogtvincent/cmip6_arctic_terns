# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 21:05:15 2021

Plot productivity projections for polar regions



@author: pearseb
"""

#%% imports

import sys
print("python version =",sys.version[:5])

import os
import numpy as np
import netCDF4 as nc
import xarray as xr
from xarrayutils.utils import linear_trend
import dask as dsk
from dask.diagnostics import ProgressBar
import numba as nb
from numba import jit

# plotting packages
import seaborn as sb
sb.set(style='ticks')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
import cmocean.cm as cmo
from cmocean.tools import lighten


print("numpy version =", np.__version__)
print("netCDF4 version =", nc.__version__)
print("xarray version =", xr.__version__)
print("dask version =", dsk.__version__)
print("numba version =", nb.__version__)
print("seaborn version =", sb.__version__)
print("matplotlib version =", sys.modules[plt.__package__].__version__)
print("cartopy version =", sys.modules[ccrs.__package__].__version__)
print("cmocean version =", sys.modules[cmo.__package__].__version__)


#%% tern area use data

#os.chdir('C://Users/pearseb/Dropbox/PostDoc/collaborations/CMIP6 hackathon 2021/data')
os.chdir('/Users/pbuchanan/Dropbox/PostDoc/collaborations/CMIP6 hackathon 2021/data')

data = nc.Dataset('bird_heatmap.nc','r')
density = data.variables['density'][...]
lon = data.variables['longitude'][...]
lat = data.variables['latitude'][...]

lon = np.ma.concatenate((lon[:,200::], lon[:,0:200]+360.0), axis=1)
density = np.ma.concatenate((density[:,200::], density[:,0:200]), axis=1)


#%% Earth System Model NPP projections 

#os.chdir('C://Users/pearseb/Dropbox/PostDoc/collaborations/CMIP6 hackathon 2021/data')
os.chdir('/Users/pbuchanan/Dropbox/PostDoc/collaborations/CMIP6 hackathon 2021/data')

access = xr.open_dataset('ETOPO_intpp_Omon_ACCESS-ESM1-5_r1i1p1f1_185001-210012_yearmonths.nc').intpp
canesm = xr.open_dataset('ETOPO_intpp_Omon_CanESM5_r1i1p2f1_185001-210012_yearmonths.nc').intpp
cesm2 = xr.open_dataset('ETOPO_intpp_Omon_CESM2_r4i1p1f1_185001-210012_yearmonths.nc').intpp
cnrm = xr.open_dataset('ETOPO_intpp_Omon_CNRM-ESM2-1_r1i1p1f2_185001-210012_yearmonths.nc').intpp
gfdlcm4 = xr.open_dataset('ETOPO_intpp_Omon_GFDL-CM4_r1i1p1f1_185001-210012_yearmonths.nc').intpp
gfdlesm = xr.open_dataset('ETOPO_intpp_Omon_GFDL-ESM4_r1i1p1f1_185001-210012_yearmonths.nc').intpp
ipsl = xr.open_dataset('ETOPO_intpp_Omon_IPSL-CM6A-LR_r1i1p1f1_185001-210012_yearmonths.nc').intpp
miroc = xr.open_dataset('ETOPO_intpp_Omon_MIROC-ES2L_r1i1p1f2_185001-210012_yearmonths.nc').intpp
mpiesm = xr.open_dataset('ETOPO_intpp_Omon_MPI-ESM1-2-HR_r1i1p1f1_185001-210012_yearmonths.nc').intpp
mriesm = xr.open_dataset('ETOPO_intpp_Omon_MRI-ESM2-0_r1i2p1f1_185001-210012_yearmonths.nc').intpp
noresm = xr.open_dataset('ETOPO_intpp_Omon_NorESM2-LM_r1i1p1f1_185001-210012_yearmonths.nc').intpp
ukesm = xr.open_dataset('ETOPO_intpp_Omon_UKESM1-0-LL_r1i1p1f2_185001-210012_yearmonths.nc').intpp

access_ssp2 = xr.open_dataset('ETOPO_intpp_Omon_ACCESS-ESM1-5_ssp245_r1i1p1f1_185001-210012_yearmonths.nc').intpp
canesm_ssp2 = xr.open_dataset('ETOPO_intpp_Omon_CanESM5_ssp245_r1i1p2f1_185001-210012_yearmonths.nc').intpp
cnrm_ssp2 = xr.open_dataset('ETOPO_intpp_Omon_CNRM-ESM2-1_ssp245_r1i1p1f2_185001-210012_yearmonths.nc').intpp
gfdlcm4_ssp2 = xr.open_dataset('ETOPO_intpp_Omon_GFDL-CM4_ssp245_r1i1p1f1_185001-210012_yearmonths.nc').intpp
gfdlesm_ssp2 = xr.open_dataset('ETOPO_intpp_Omon_GFDL-ESM4_ssp245_r1i1p1f1_185001-210012_yearmonths.nc').intpp
ipsl_ssp2 = xr.open_dataset('ETOPO_intpp_Omon_IPSL-CM6A-LR_ssp245_r1i1p1f1_185001-210012_yearmonths.nc').intpp
miroc_ssp2 = xr.open_dataset('ETOPO_intpp_Omon_MIROC-ES2L_ssp245_r1i1p1f2_185001-210012_yearmonths.nc').intpp
mpiesm_ssp2 = xr.open_dataset('ETOPO_intpp_Omon_MPI-ESM1-2-HR_ssp245_r1i1p1f1_185001-210012_yearmonths.nc').intpp
noresm_ssp2 = xr.open_dataset('ETOPO_intpp_Omon_NorESM2-LM_ssp245_r1i1p1f1_185001-210012_yearmonths.nc').intpp
ukesm_ssp2 = xr.open_dataset('ETOPO_intpp_Omon_UKESM1-0-LL_ssp245_r1i1p1f2_185001-210012_yearmonths.nc').intpp

years2100 = np.arange(1850.5,2100.6,1)
months = np.arange(1,13,1)

### assign coordinate to record dimension, redefine the months dimension, and chunk the dataset along months
chunky = {'record':251, 'time':1, 'ETOPO60Y':180, 'ETOPO60X':360}

access = access.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
canesm = canesm.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
cesm2 = cesm2.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
cnrm = cnrm.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
gfdlcm4 = gfdlcm4.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
gfdlesm = gfdlesm.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
ipsl = ipsl.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
miroc = miroc.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
mpiesm = mpiesm.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
mriesm = mriesm.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
noresm = noresm.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
ukesm = ukesm.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)

access_ssp2 = access_ssp2.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
canesm_ssp2 = canesm_ssp2.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
cnrm_ssp2 = cnrm_ssp2.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
gfdlcm4_ssp2 = gfdlcm4_ssp2.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
gfdlesm_ssp2 = gfdlesm_ssp2.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
ipsl_ssp2 = ipsl_ssp2.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
miroc_ssp2 = miroc_ssp2.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
mpiesm_ssp2 = mpiesm_ssp2.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
noresm_ssp2 = noresm_ssp2.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
ukesm_ssp2 = ukesm_ssp2.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)


# check size of the arrays in Gb and their chunks in Mb
print("Total size of array in Gb", ukesm.nbytes * 1e-9)
print("Total size of chunk in Mb", ukesm.sel(time=1).nbytes * 1e-6)


#%% calculate the multi-model mean and standard deviation

ProgressBar().register()

# check models first
    # concatenate them all into one xarray by adding model as a new dimension and chunck along that dimension
npp_multi = xr.concat((access,cesm2,canesm,cnrm,gfdlcm4,gfdlesm,ipsl,miroc,mpiesm,mriesm,noresm,ukesm), dim='model')
ssp2_multi = xr.concat((access_ssp2,canesm_ssp2,cnrm_ssp2,gfdlcm4_ssp2,gfdlesm_ssp2,ipsl_ssp2,miroc_ssp2,mpiesm_ssp2,noresm,ukesm_ssp2), dim='model')

models = ['access', 'canesm', 'cesm2', 'cnrm', 'gfdl-cm4', 'gfdl-esm', 'ipsl', 'miroc', 'mpi-esm', 'mri-esm', 'noresm', 'ukesm']
models_ssp2 = ['access', 'canesm', 'cnrm', 'gfdl-cm4', 'gfdl-esm', 'ipsl', 'miroc', 'mpi-esm', 'noresm', 'ukesm']
chunky = {'model':1, 'record':251, 'time':1, 'ETOPO60Y':180, 'ETOPO60X':360}

npp_multi = npp_multi.assign_coords(model=('model', models)).chunk(chunks=chunky)
ssp2_multi = ssp2_multi.assign_coords(model=('model', models_ssp2)).chunk(chunks=chunky)

# take average of the models and the inter-model standard deviation and load into memory
#npp_mean = npp_multi.mean(dim='model').load().values
#npp_std = npp_multi.std(dim='model').load().values


#%% calculate stuff means, changes and how many models agree on the direction of change (NH and SH)

# find pre industrial mean
npp_past_sh_mean = npp_multi.sel(record = slice(1850.5,1950.5)).mean(dim='record').sel(time = slice(1,3)).mean(dim='time').mean(dim='model')
npp_past_nh_mean = npp_multi.sel(record = slice(1850.5,1950.5)).mean(dim='record').sel(time = slice(7,9)).mean(dim='time').mean(dim='model')

# find the mean change at now (1995-2014) and in future (2081-2100) compared with preindustrial period from Jan-Mar
npp_hist_change = npp_multi.sel(record = slice(1995.5,2014.5)).mean(dim='record') - npp_multi.sel(record = slice(1850.5,1950.5)).mean(dim='record')
npp_futu_change = npp_multi.sel(record = slice(2081.5,2100.5)).mean(dim='record') - npp_multi.sel(record = slice(1850.5,1950.5)).mean(dim='record')
npp_hist_sh_mean = npp_hist_change.mean(dim='model').sel(time = slice(1,3)).mean(dim='time')
npp_futu_sh_mean = npp_futu_change.mean(dim='model').sel(time = slice(1,3)).mean(dim='time')
npp_hist_nh_mean = npp_hist_change.mean(dim='model').sel(time = slice(7,9)).mean(dim='time')
npp_futu_nh_mean = npp_futu_change.mean(dim='model').sel(time = slice(7,9)).mean(dim='time')

# determine how many models agree in the direction of change
npp_hist_sh_pos = npp_hist_change.sel(time = slice(1,3)).mean(dim='time') > 0.0
npp_hist_sh_neg = npp_hist_change.sel(time = slice(1,3)).mean(dim='time') < 0.0
npp_hist_sh_possum = npp_hist_sh_pos.sum(dim='model')
npp_hist_sh_negsum = npp_hist_sh_neg.sum(dim='model')

npp_futu_sh_pos = npp_futu_change.sel(time = slice(1,3)).mean(dim='time') > 0.0
npp_futu_sh_neg = npp_futu_change.sel(time = slice(1,3)).mean(dim='time') < 0.0
npp_futu_sh_possum = npp_futu_sh_pos.sum(dim='model')
npp_futu_sh_negsum = npp_futu_sh_neg.sum(dim='model')

npp_hist_nh_pos = npp_hist_change.sel(time = slice(7,9)).mean(dim='time') > 0.0
npp_hist_nh_neg = npp_hist_change.sel(time = slice(7,9)).mean(dim='time') < 0.0
npp_hist_nh_possum = npp_hist_nh_pos.sum(dim='model')
npp_hist_nh_negsum = npp_hist_nh_neg.sum(dim='model')

npp_futu_nh_pos = npp_futu_change.sel(time = slice(7,9)).mean(dim='time') > 0.0
npp_futu_nh_neg = npp_futu_change.sel(time = slice(7,9)).mean(dim='time') < 0.0
npp_futu_nh_possum = npp_futu_nh_pos.sum(dim='model')
npp_futu_nh_negsum = npp_futu_nh_neg.sum(dim='model')


#%% do for SSP2-4.5

ssp2_past_sh_mean = ssp2_multi.sel(record = slice(1850.5,1950.5)).mean(dim='record').sel(time = slice(1,3)).mean(dim='time').mean(dim='model')
ssp2_past_nh_mean = ssp2_multi.sel(record = slice(1850.5,1950.5)).mean(dim='record').sel(time = slice(7,9)).mean(dim='time').mean(dim='model')

# find the mean change at now (1995-2014) and in future (2081-2100) compared with preindustrial period from Jan-Mar
ssp2_hist_change = ssp2_multi.sel(record = slice(1995.5,2014.5)).mean(dim='record') - ssp2_multi.sel(record = slice(1850.5,1950.5)).mean(dim='record')
ssp2_futu_change = ssp2_multi.sel(record = slice(2081.5,2100.5)).mean(dim='record') - ssp2_multi.sel(record = slice(1850.5,1950.5)).mean(dim='record')
ssp2_hist_sh_mean = ssp2_hist_change.mean(dim='model').sel(time = slice(1,3)).mean(dim='time')
ssp2_futu_sh_mean = ssp2_futu_change.mean(dim='model').sel(time = slice(1,3)).mean(dim='time')
ssp2_hist_nh_mean = ssp2_hist_change.mean(dim='model').sel(time = slice(7,9)).mean(dim='time')
ssp2_futu_nh_mean = ssp2_futu_change.mean(dim='model').sel(time = slice(7,9)).mean(dim='time')

# determine how many models agree in the direction of change
ssp2_hist_sh_pos = ssp2_hist_change.sel(time = slice(1,3)).mean(dim='time') > 0.0
ssp2_hist_sh_neg = ssp2_hist_change.sel(time = slice(1,3)).mean(dim='time') < 0.0
ssp2_hist_sh_possum = ssp2_hist_sh_pos.sum(dim='model')
ssp2_hist_sh_negsum = ssp2_hist_sh_neg.sum(dim='model')

ssp2_futu_sh_pos = ssp2_futu_change.sel(time = slice(1,3)).mean(dim='time') > 0.0
ssp2_futu_sh_neg = ssp2_futu_change.sel(time = slice(1,3)).mean(dim='time') < 0.0
ssp2_futu_sh_possum = ssp2_futu_sh_pos.sum(dim='model')
ssp2_futu_sh_negsum = ssp2_futu_sh_neg.sum(dim='model')

ssp2_hist_nh_pos = ssp2_hist_change.sel(time = slice(7,9)).mean(dim='time') > 0.0
ssp2_hist_nh_neg = ssp2_hist_change.sel(time = slice(7,9)).mean(dim='time') < 0.0
ssp2_hist_nh_possum = ssp2_hist_nh_pos.sum(dim='model')
ssp2_hist_nh_negsum = ssp2_hist_nh_neg.sum(dim='model')

ssp2_futu_nh_pos = ssp2_futu_change.sel(time = slice(7,9)).mean(dim='time') > 0.0
ssp2_futu_nh_neg = ssp2_futu_change.sel(time = slice(7,9)).mean(dim='time') < 0.0
ssp2_futu_nh_possum = ssp2_futu_nh_pos.sum(dim='model')
ssp2_futu_nh_negsum = ssp2_futu_nh_neg.sum(dim='model')



#%% prepare things for plotting

proj_nh = ccrs.Orthographic(central_longitude=-20.0, central_latitude=90.0, globe=None)
proj_sh = ccrs.Orthographic(central_longitude=-20.0, central_latitude=-90.0, globe=None)
proj = ccrs.Robinson(central_longitude=-20.0)

levs1 = np.arange(0,101,5)*0.01
levs2 = np.arange(-20,21,2)*0.01

colmap1 = lighten(cmo.haline, 0.8)
colmap2 = lighten(cmo.delta, 0.8)
colmap3 = lighten(cmo.amp, 0.8)

fstic = 13
fslab = 15

hatching = ['....']

scaler = 86400*12   # convert mol C m-2 s-1 to g C m-2 day-1

### wrap longtiudes
density, lons = add_cyclic_point(density, coord=lon[0,:])
npp_past_nh_mean_vals, lons = add_cyclic_point(npp_past_nh_mean, coord=lon[0,:])
npp_past_sh_mean_vals, lons = add_cyclic_point(npp_past_sh_mean, coord=lon[0,:])
npp_hist_nh_mean_vals, lons = add_cyclic_point(npp_hist_nh_mean, coord=lon[0,:])
npp_hist_nh_possum_vals, lons = add_cyclic_point(npp_hist_nh_possum, coord=lon[0,:])
npp_hist_nh_negsum_vals, lons = add_cyclic_point(npp_hist_nh_negsum, coord=lon[0,:])
npp_hist_sh_mean_vals, lons = add_cyclic_point(npp_hist_sh_mean, coord=lon[0,:])
npp_hist_sh_possum_vals, lons = add_cyclic_point(npp_hist_sh_possum, coord=lon[0,:])
npp_hist_sh_negsum_vals, lons = add_cyclic_point(npp_hist_sh_negsum, coord=lon[0,:])
npp_futu_nh_mean_vals, lons = add_cyclic_point(npp_futu_nh_mean, coord=lon[0,:])
npp_futu_nh_possum_vals, lons = add_cyclic_point(npp_futu_nh_possum, coord=lon[0,:])
npp_futu_nh_negsum_vals, lons = add_cyclic_point(npp_futu_nh_negsum, coord=lon[0,:])
npp_futu_sh_mean_vals, lons = add_cyclic_point(npp_futu_sh_mean, coord=lon[0,:])
npp_futu_sh_possum_vals, lons = add_cyclic_point(npp_futu_sh_possum, coord=lon[0,:])
npp_futu_sh_negsum_vals, lons = add_cyclic_point(npp_futu_sh_negsum, coord=lon[0,:])

ssp2_past_nh_mean_vals, lons = add_cyclic_point(ssp2_past_nh_mean, coord=lon[0,:])
ssp2_past_sh_mean_vals, lons = add_cyclic_point(ssp2_past_sh_mean, coord=lon[0,:])
ssp2_hist_nh_mean_vals, lons = add_cyclic_point(ssp2_hist_nh_mean, coord=lon[0,:])
ssp2_hist_nh_possum_vals, lons = add_cyclic_point(ssp2_hist_nh_possum, coord=lon[0,:])
ssp2_hist_nh_negsum_vals, lons = add_cyclic_point(ssp2_hist_nh_negsum, coord=lon[0,:])
ssp2_hist_sh_mean_vals, lons = add_cyclic_point(ssp2_hist_sh_mean, coord=lon[0,:])
ssp2_hist_sh_possum_vals, lons = add_cyclic_point(ssp2_hist_sh_possum, coord=lon[0,:])
ssp2_hist_sh_negsum_vals, lons = add_cyclic_point(ssp2_hist_sh_negsum, coord=lon[0,:])
ssp2_futu_nh_mean_vals, lons = add_cyclic_point(ssp2_futu_nh_mean, coord=lon[0,:])
ssp2_futu_nh_possum_vals, lons = add_cyclic_point(ssp2_futu_nh_possum, coord=lon[0,:])
ssp2_futu_nh_negsum_vals, lons = add_cyclic_point(ssp2_futu_nh_negsum, coord=lon[0,:])
ssp2_futu_sh_mean_vals, lons = add_cyclic_point(ssp2_futu_sh_mean, coord=lon[0,:])
ssp2_futu_sh_possum_vals, lons = add_cyclic_point(ssp2_futu_sh_possum, coord=lon[0,:])
ssp2_futu_sh_negsum_vals, lons = add_cyclic_point(ssp2_futu_sh_negsum, coord=lon[0,:])

lon, lat = np.meshgrid(lons,lat[:,0])


#%% plot the multi-model mean in productivity changes in the polar regions, hashing where < 75% of models agree (9/12)

fig = plt.figure(figsize=(6,10))
gs = GridSpec(4,2)

ax1 = plt.subplot(gs[0,0], projection=proj_nh)
p1 = plt.contourf(lon, lat, npp_past_nh_mean_vals*scaler, transform=ccrs.PlateCarree(), cmap=colmap1, levels=levs1, vmin=np.min(levs1), vmax=np.max(levs1), extend='max')
c1 = plt.contour(lon, lat, density, transform=ccrs.PlateCarree(), levels=np.arange(0.5,1.51,0.5), linewidths=0.75)
ax1.add_feature(cfeature.LAND, color='silver', zorder=2)
ax1.coastlines(zorder=2)

ax2 = plt.subplot(gs[0,1], projection=proj_sh)
p2 = plt.contourf(lon, lat, npp_past_sh_mean_vals*scaler, transform=ccrs.PlateCarree(), cmap=colmap1, levels=levs1, vmin=np.min(levs1), vmax=np.max(levs1), extend='max')
c2 = plt.contour(lon, lat, density, transform=ccrs.PlateCarree(), levels=np.arange(0.5,1.51,0.5), linewidths=0.75)
ax2.add_feature(cfeature.LAND, color='silver', zorder=2)
ax2.coastlines(zorder=2)

ax3 = plt.subplot(gs[1,0], projection=proj_nh)
p3 = plt.contourf(lon, lat, ssp2_futu_nh_mean_vals*scaler, transform=ccrs.PlateCarree(), cmap=colmap2, levels=levs2, vmin=np.min(levs2), vmax=np.max(levs2), extend='both')
c3 = plt.contourf(lon, lat, ssp2_futu_nh_possum_vals, transform=ccrs.PlateCarree(), colors='none', levels=[9,12], hatches=hatching)
c3 = plt.contourf(lon, lat, ssp2_futu_nh_negsum_vals, transform=ccrs.PlateCarree(), colors='none', levels=[9,12], hatches=hatching)
ax3.add_feature(cfeature.LAND, color='silver', zorder=2)
ax3.coastlines(zorder=2)

ax4 = plt.subplot(gs[1,1], projection=proj_sh)
p4 = plt.contourf(lon, lat, ssp2_futu_sh_mean_vals*scaler, transform=ccrs.PlateCarree(), cmap=colmap2, levels=levs2, vmin=np.min(levs2), vmax=np.max(levs2), extend='both')
c4 = plt.contourf(lon, lat, ssp2_futu_sh_possum_vals, transform=ccrs.PlateCarree(), colors='none', levels=[9,12], hatches=hatching)
c4 = plt.contourf(lon, lat, ssp2_futu_sh_negsum_vals, transform=ccrs.PlateCarree(), colors='none', levels=[9,12], hatches=hatching)
ax4.add_feature(cfeature.LAND, color='silver', zorder=2)
ax4.coastlines(zorder=2)

ax5 = plt.subplot(gs[2,0], projection=proj_nh)
p5 = plt.contourf(lon, lat, npp_futu_nh_mean_vals*scaler, transform=ccrs.PlateCarree(), cmap=colmap2, levels=levs2, vmin=np.min(levs2), vmax=np.max(levs2), extend='both')
c5 = plt.contourf(lon, lat, npp_futu_nh_possum_vals, transform=ccrs.PlateCarree(), colors='none', levels=[9,12], hatches=hatching)
c5 = plt.contourf(lon, lat, npp_futu_nh_negsum_vals, transform=ccrs.PlateCarree(), colors='none', levels=[9,12], hatches=hatching)
ax5.add_feature(cfeature.LAND, color='silver', zorder=2)
ax5.coastlines(zorder=2)

ax6 = plt.subplot(gs[2,1], projection=proj_sh)
p6 = plt.contourf(lon, lat, npp_futu_sh_mean_vals*scaler, transform=ccrs.PlateCarree(), cmap=colmap2, levels=levs2, vmin=np.min(levs2), vmax=np.max(levs2), extend='both')
c6 = plt.contourf(lon, lat, npp_futu_sh_possum_vals, transform=ccrs.PlateCarree(), colors='none', levels=[9,12], hatches=hatching)
c6 = plt.contourf(lon, lat, npp_futu_sh_negsum_vals, transform=ccrs.PlateCarree(), colors='none', levels=[9,12], hatches=hatching)
ax6.add_feature(cfeature.LAND, color='silver', zorder=2)
ax6.coastlines(zorder=2)

ax7 = plt.subplot(gs[3,0], projection=proj_nh)
p7 = plt.contourf(lon, lat, np.abs(npp_futu_nh_mean_vals) / npp_past_nh_mean_vals*100, transform=ccrs.PlateCarree(), cmap=colmap3, levels=levs1*100, vmin=np.min(levs1)*100, vmax=np.max(levs1)*100, extend='max')
c7 = plt.contour(lon, lat, np.abs(npp_futu_nh_mean_vals) / npp_past_nh_mean_vals*100, transform=ccrs.PlateCarree(), levels=np.arange(0,101,20), colors='k', linewidths=0.75)
ax7.add_feature(cfeature.LAND, color='silver', zorder=2)
ax7.coastlines(zorder=2)

ax8 = plt.subplot(gs[3,1], projection=proj_sh)
p8 = plt.contourf(lon, lat, np.abs(npp_futu_nh_mean_vals) / npp_past_nh_mean_vals*100, transform=ccrs.PlateCarree(), cmap=colmap3, levels=levs1*100, vmin=np.min(levs1)*100, vmax=np.max(levs1)*100, extend='max')
c8 = plt.contour(lon, lat, np.abs(npp_futu_nh_mean_vals) / npp_past_nh_mean_vals*100, transform=ccrs.PlateCarree(), levels=np.arange(0,101,20), colors='k', linewidths=0.75)
ax8.add_feature(cfeature.LAND, color='silver', zorder=2)
ax8.coastlines(zorder=2)

plt.subplots_adjust(right=0.835, top=0.95, bottom=0.05, left=0.01, hspace=0.02, wspace=-0.3)



cbax1 = fig.add_axes([0.78, 0.75, 0.035, 0.175])
cbar1 = plt.colorbar(p1, cax=cbax1, orientation='vertical', ticks=levs1[::4])
cbar1.ax.set_ylabel('Historical NPP\n(g C m$^{-2}$ day$^{-1}$)', fontsize=fslab)
cbar1.ax.tick_params(labelsize=fstic)

cbax2 = fig.add_axes([0.78, 0.525, 0.035, 0.175])
cbar2 = plt.colorbar(p3, cax=cbax2, orientation='vertical', ticks=levs2[::4])
cbar2.ax.set_ylabel('$\Delta$NPP\n(SSP2-4.5)', fontsize=fslab)
cbar2.ax.tick_params(labelsize=fstic)

cbax3 = fig.add_axes([0.78, 0.3, 0.035, 0.175])
cbar3 = plt.colorbar(p5, cax=cbax3, orientation='vertical', ticks=levs2[::4])
cbar3.ax.set_ylabel('$\Delta$NPP\n(SSP5-8.5)', fontsize=fslab)
cbar3.ax.tick_params(labelsize=fstic)

cbax4 = fig.add_axes([0.78, 0.075, 0.035, 0.175])
cbar4 = plt.colorbar(p7, cax=cbax4, orientation='vertical', ticks=levs1[::4]*100)
cbar4.ax.set_ylabel('$\Delta$NPP (%)\n(SSP5-8.5)', fontsize=fslab)
cbar4.ax.tick_params(labelsize=fstic)


xx = 0.05; yy = 0.95
ax1.text(xx,yy,'a', transform=ax1.transAxes, fontsize=fslab+2, fontweight='bold', ha='center', va='center')
ax2.text(xx,yy,'b', transform=ax2.transAxes, fontsize=fslab+2, fontweight='bold', ha='center', va='center')
ax3.text(xx,yy,'c', transform=ax3.transAxes, fontsize=fslab+2, fontweight='bold', ha='center', va='center')
ax4.text(xx,yy,'d', transform=ax4.transAxes, fontsize=fslab+2, fontweight='bold', ha='center', va='center')
ax5.text(xx,yy,'e', transform=ax5.transAxes, fontsize=fslab+2, fontweight='bold', ha='center', va='center')
ax6.text(xx,yy,'f', transform=ax6.transAxes, fontsize=fslab+2, fontweight='bold', ha='center', va='center')
ax7.text(xx,yy,'g', transform=ax7.transAxes, fontsize=fslab+2, fontweight='bold', ha='center', va='center')
ax8.text(xx,yy,'h', transform=ax8.transAxes, fontsize=fslab+2, fontweight='bold', ha='center', va='center')

#xx = 1.05; yy = 1.075
#ax1.text(xx,yy,'Historical', transform=ax1.transAxes, fontsize=fstic, fontweight='bold', ha='center', va='center')
#ax3.text(xx,yy,'SSP2-4.5', transform=ax3.transAxes, fontsize=fstic, fontweight='bold', ha='center', va='center')
#ax5.text(xx,yy,'SSP5-8.5', transform=ax5.transAxes, fontsize=fstic, fontweight='bold', ha='center', va='center')
#ax7.text(xx,yy,'SSP5-8.5', transform=ax7.transAxes, fontsize=fstic, fontweight='bold', ha='center', va='center')


#%% save figure

#os.chdir('C://Users/pearseb/Dropbox/PostDoc/CMIP6 hackathon/figures')
os.chdir('/Users/pbuchanan/Dropbox/PostDoc/collaborations/CMIP6 hackathon 2021/figures')
fig.savefig('polar_productivity.png', dpi=300, bbox_inches='tight')

