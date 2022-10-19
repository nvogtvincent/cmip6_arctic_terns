# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 14:30:25 2021

Purpose
-------
    Calculate the signal:noise ratio of 12 state-of-the-art Earth System Models
    weighted by the area use of Arctic Terns in four major foraging areas


Steps 
-----
    1. Load tern area-use data and ESM productivity projections
    2. Calculate the multi-model average and standard deviation for each month
    3. Define the signal:noise ratio calculation
    4. Calculate the signal:noise ratio for all months on the multi-model mean and 1 standard deviation around this mean


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

access_85 = xr.open_dataset('ETOPO_intpp_Omon_ACCESS-ESM1-5_r1i1p1f1_185001-210012_yearmonths.nc').intpp
canesm_85 = xr.open_dataset('ETOPO_intpp_Omon_CanESM5_r1i1p2f1_185001-210012_yearmonths.nc').intpp
cesm2_85 = xr.open_dataset('ETOPO_intpp_Omon_CESM2_r4i1p1f1_185001-210012_yearmonths.nc').intpp
cnrm_85 = xr.open_dataset('ETOPO_intpp_Omon_CNRM-ESM2-1_r1i1p1f2_185001-210012_yearmonths.nc').intpp
gfdlcm4_85 = xr.open_dataset('ETOPO_intpp_Omon_GFDL-CM4_r1i1p1f1_185001-210012_yearmonths.nc').intpp
gfdlesm_85 = xr.open_dataset('ETOPO_intpp_Omon_GFDL-ESM4_r1i1p1f1_185001-210012_yearmonths.nc').intpp
ipsl_85 = xr.open_dataset('ETOPO_intpp_Omon_IPSL-CM6A-LR_r1i1p1f1_185001-210012_yearmonths.nc').intpp
miroc_85 = xr.open_dataset('ETOPO_intpp_Omon_MIROC-ES2L_r1i1p1f2_185001-210012_yearmonths.nc').intpp
mpiesm_85 = xr.open_dataset('ETOPO_intpp_Omon_MPI-ESM1-2-HR_r1i1p1f1_185001-210012_yearmonths.nc').intpp
mriesm_85 = xr.open_dataset('ETOPO_intpp_Omon_MRI-ESM2-0_r1i2p1f1_185001-210012_yearmonths.nc').intpp
noresm_85 = xr.open_dataset('ETOPO_intpp_Omon_NorESM2-LM_r1i1p1f1_185001-210012_yearmonths.nc').intpp
ukesm_85 = xr.open_dataset('ETOPO_intpp_Omon_UKESM1-0-LL_r1i1p1f2_185001-210012_yearmonths.nc').intpp

access_45 = xr.open_dataset('ETOPO_intpp_Omon_ACCESS-ESM1-5_ssp245_r1i1p1f1_185001-210012_yearmonths.nc').intpp
canesm_45 = xr.open_dataset('ETOPO_intpp_Omon_CanESM5_ssp245_r1i1p2f1_185001-210012_yearmonths.nc').intpp
cnrm_45 = xr.open_dataset('ETOPO_intpp_Omon_CNRM-ESM2-1_ssp245_r1i1p1f2_185001-210012_yearmonths.nc').intpp
gfdlcm4_45 = xr.open_dataset('ETOPO_intpp_Omon_GFDL-CM4_ssp245_r1i1p1f1_185001-210012_yearmonths.nc').intpp
gfdlesm_45 = xr.open_dataset('ETOPO_intpp_Omon_GFDL-ESM4_ssp245_r1i1p1f1_185001-210012_yearmonths.nc').intpp
ipsl_45 = xr.open_dataset('ETOPO_intpp_Omon_IPSL-CM6A-LR_ssp245_r1i1p1f1_185001-210012_yearmonths.nc').intpp
miroc_45 = xr.open_dataset('ETOPO_intpp_Omon_MIROC-ES2L_ssp245_r1i1p1f2_185001-210012_yearmonths.nc').intpp
mpiesm_45 = xr.open_dataset('ETOPO_intpp_Omon_MPI-ESM1-2-HR_ssp245_r1i1p1f1_185001-210012_yearmonths.nc').intpp
noresm_45 = xr.open_dataset('ETOPO_intpp_Omon_NorESM2-LM_ssp245_r1i1p1f1_185001-210012_yearmonths.nc').intpp
ukesm_45 = xr.open_dataset('ETOPO_intpp_Omon_UKESM1-0-LL_ssp245_r1i1p1f2_185001-210012_yearmonths.nc').intpp

years2100 = np.arange(1850.5,2100.6,1)
months = np.arange(1,13,1)

### assign coordinate to record dimension, redefine the months dimension, and chunk the dataset along months
chunky = {'record':251, 'time':1, 'ETOPO60Y':180, 'ETOPO60X':360}

access_85 = access_85.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
canesm_85 = canesm_85.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
cesm2_85 = cesm2_85.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
cnrm_85 = cnrm_85.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
gfdlcm4_85 = gfdlcm4_85.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
gfdlesm_85 = gfdlesm_85.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
ipsl_85 = ipsl_85.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
miroc_85 = miroc_85.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
mpiesm_85 = mpiesm_85.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
mriesm_85 = mriesm_85.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
noresm_85 = noresm_85.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
ukesm_85 = ukesm_85.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)

access_45 = access_45.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
canesm_45 = canesm_45.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
cnrm_45 = cnrm_45.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
gfdlcm4_45 = gfdlcm4_45.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
gfdlesm_45 = gfdlesm_45.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
ipsl_45 = ipsl_45.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
miroc_45 = miroc_45.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
mpiesm_45 = mpiesm_45.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
noresm_45 = noresm_45.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)
ukesm_45 = ukesm_45.assign_coords(record=('record', years2100), time=('time', months)).chunk(chunks=chunky)


# check size of the arrays in Gb and their chunks in Mb
print("Total size of array in Gb", ukesm_85.nbytes * 1e-9)
print("Total size of chunk in Mb", ukesm_85.sel(time=1).nbytes * 1e-6)


#%% calculate the multi-model mean and standard deviation

ProgressBar().register()

# check models first
    # concatenate them all into one xarray by adding model as a new dimension and chunck along that dimension
npp_85 = xr.concat((access_85,cesm2_85,canesm_85,cnrm_85,gfdlcm4_85,gfdlesm_85,ipsl_85,miroc_85,mpiesm_85,mriesm_85,noresm_85,ukesm_85), dim='model')
models = ['access', 'canesm', 'cesm2', 'cnrm', 'gfdl-cm4', 'gfdl-esm', 'ipsl', 'miroc', 'mpi-esm', 'mri-esm', 'noresm', 'ukesm']
chunky = {'model':1, 'record':251, 'time':1, 'ETOPO60Y':180, 'ETOPO60X':360}
npp_85 = npp_85.assign_coords(model=('model', models)).chunk(chunks=chunky)

    # concatenate them all into one xarray by adding model as a new dimension and chunck along that dimension
npp_45 = xr.concat((access_45,canesm_45,cnrm_45,gfdlcm4_45,gfdlesm_45,ipsl_45,miroc_45,mpiesm_45,noresm_45,ukesm_45), dim='model')
models = ['access', 'canesm', 'cnrm', 'gfdl-cm4', 'gfdl-esm', 'ipsl', 'miroc', 'mpi-esm', 'noresm', 'ukesm']
chunky = {'model':1, 'record':251, 'time':1, 'ETOPO60Y':180, 'ETOPO60X':360}
npp_45 = npp_45.assign_coords(model=('model', models)).chunk(chunks=chunky)

# take average of the models and the inter-model standard deviation and load into memory
npp_85_mean = npp_85.mean(dim='model').load().values
npp_85_std = npp_85.std(dim='model').load().values

# take average of the models and the inter-model standard deviation and load into memory
npp_45_mean = npp_45.mean(dim='model').load().values
npp_45_std = npp_45.std(dim='model').load().values


#%% correct array by setting upper and lower bounds and changing units (mol C m-2 s-1 --> g C m-2 yr-1)

npp_45_mean = npp_45_mean*86400*365*12
npp_45_mean[npp_45_mean < 1e-16] = np.nan
npp_45_mean[npp_45_mean > 1000] = np.nan
npp_45_std = npp_45_std*86400*365*12 
npp_45_std[np.isnan(npp_45_mean)] = np.nan

npp_85_mean = npp_85_mean*86400*365*12
npp_85_mean[npp_85_mean < 1e-16] = np.nan
npp_85_mean[npp_85_mean > 1000] = np.nan
npp_85_std = npp_85_std*86400*365*12 
npp_85_std[np.isnan(npp_85_mean)] = np.nan


#%% get signal to noise at each grid cell at each year and month   

@jit(nopython=True)
def lin_regression(arr1, arr2):
    m_x = np.mean(arr1)
    m_y = np.mean(arr2)
    pr = np.sum( (arr1 - m_x) * (arr2 - m_y) ) / ( np.sum( (arr1-m_x)**2 ) * np.sum( (arr2 - m_y)**2 ) )**0.5
    s_x = np.sqrt( np.sum( (arr1 - m_x)**2 ) / ( len(arr1) ) )
    s_y = np.sqrt( np.sum( (arr2 - m_y)**2 ) / ( len(arr2) ) )
    slope = pr * (s_y/s_x)
    inter = m_y - slope * m_x
    return slope, inter

@jit(nopython=True)
def get_signal_to_noise(ny, data):
    
    # initialise the array to save    
    tre_array = np.zeros(np.shape(data))
    std_array = np.zeros(np.shape(data))
    s2n_array = np.zeros(np.shape(data))
    
    for yr,year in enumerate(years2100):
        print(year)
        # define end year
        stop = yr+ny
        time = np.arange(1,ny+1,1)
    
        if (year > 2070):
            print("end")
            break
        else:
            for m in np.arange(len(data[0,:,0,0])):
                for y in np.arange(len(data[0,0,:,0])):
                    for x in np.arange(len(data[0,0,0,:])):
                        # find linear trend
                        slope, inter = lin_regression(time, data[yr:stop,m,y,x])
                        # find noise from the error around the linear trend
                        error = (slope * time + inter) - data[yr:stop,m,y,x]
                        # save arrays
                        tre_array[yr,m,y,x] = slope
                        std_array[yr,m,y,x] = np.std(error)
                        s2n_array[yr,m,y,x] = (ny*slope) / np.std(error)
                        
    return s2n_array, tre_array, std_array


#%% collect signal to noise array

ny = 30
data = npp_45_mean
(npp_45_mean_s2n, npp_45_mean_tre, npp_45_mean_std) = get_signal_to_noise(ny, data)

ny = 30
data = npp_85_mean
(npp_85_mean_s2n, npp_85_mean_tre, npp_85_mean_std) = get_signal_to_noise(ny, data)


#%% check the trends and S:N at a point in the north atlantic

fstic = 13
fslab = 15

fig = plt.figure(figsize=(10,9.5))
gs = GridSpec(4,1)

la = 138
lo = 310

ax1 = plt.subplot(gs[0])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.tick_params(labelsize=fstic)
plt.plot(np.arange(1850,2101,1),npp_85_mean[:,0,la,lo], color='k')
for y,year in enumerate(np.arange(1850,2101,10)):
    yr = y*10
    inter = np.mean(npp_85_mean[yr:yr+30,0,la,lo]) - npp_85_mean_tre[yr,0,la,lo] * 15.5
    plt.plot((year,year+30), (inter, inter + npp_85_mean_tre[yr,0,la,lo]*30), color='royalblue', alpha=0.5)
plt.ylabel('NPP', fontsize=fslab)
plt.xlim(1850,2100)
plt.ylim(25,60)

ax2 = plt.subplot(gs[1])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.tick_params(labelsize=fstic)
plt.plot(np.arange(1850,2101,1),npp_85_mean_tre[:,0,la,lo]*30, color='royalblue')
plt.plot((1850,2100),(0,0),'k--')
plt.ylabel('Trend*30', fontsize=fslab)
plt.xlim(1850,2100)

ax3 = plt.subplot(gs[2])
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.tick_params(labelsize=fstic)
plt.plot(np.arange(1850,2101,1),npp_85_mean_std[:,0,la,lo], color='royalblue')
plt.plot((1850,2100),(0,0),'k--')
plt.ylabel('Standard deviation', fontsize=fslab)
plt.xlim(1850,2100)

ax4 = plt.subplot(gs[3])
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.tick_params(labelsize=fstic)
plt.plot(np.arange(1850,2101,1),npp_85_mean_s2n[:,0,la,lo], color='royalblue')
plt.plot((1850,2100),(0,0),'k:')
plt.plot((1850,2100),(-1,-1),'k--')
plt.plot((1850,2100),(1,1),'k--')
plt.ylabel('Signal : Noise', fontsize=fslab)
plt.xlabel('year', fontsize=fslab)
plt.xlim(1850,2100)


xx = 0.025; yy=0.95
plt.text(xx,yy,'a', fontweight='bold', va='center', ha='center', transform=ax1.transAxes, fontsize=fslab+2)
plt.text(xx,yy,'b', fontweight='bold', va='center', ha='center', transform=ax2.transAxes, fontsize=fslab+2)
plt.text(xx,yy,'c', fontweight='bold', va='center', ha='center', transform=ax3.transAxes, fontsize=fslab+2)
plt.text(xx,yy,'d', fontweight='bold', va='center', ha='center', transform=ax4.transAxes, fontsize=fslab+2)


#%% save

#os.chdir('C://Users/pearseb/Dropbox/PostDoc/collaborations/CMIP6 hackathon 2021/figures')
os.chdir('/Users/pbuchanan/Dropbox/PostDoc/collaborations/CMIP6 hackathon 2021/figures')
fig.savefig('suppfig-signaltonoise.png', dpi=300, bbox_inches='tight')


#%% calculate the mean NPP trends and signal:noise weighted by tern foraging use area

npp_45_mean = np.ma.masked_where(np.isnan(npp_45_mean), npp_45_mean)
mask_45 = np.ma.getmask(npp_45_mean)
npp_45_std = np.ma.masked_where(mask_45, npp_45_std)

npp_85_mean = np.ma.masked_where(np.isnan(npp_85_mean), npp_85_mean)
mask_85 = np.ma.getmask(npp_85_mean)
npp_85_std = np.ma.masked_where(mask_85, npp_85_std)


# North Atlantic
la1 = 125; la2 = 155
lo1 = 250; lo2 = 340
print("longitude =",lon[0,lo1],lon[0,lo2])
print("latitude =",lat[la1,0],lat[la2,0])
mask_45_NA = mask_45[:,:,la1:la2,lo1:lo2]
density_45_NA = np.ma.masked_where(mask_45_NA[0,0,...], density[la1:la2,lo1:lo2])
npp_45_std_NA = np.ma.masked_where(mask_45_NA, npp_45_std[:,:,la1:la2,lo1:lo2])
npp_45_mean_NA = np.ma.masked_where(mask_45_NA, npp_45_mean[:,:,la1:la2,lo1:lo2])
npp_45_mean_s2n_NA = np.ma.masked_where(mask_45_NA, npp_45_mean_s2n[:,:,la1:la2,lo1:lo2])
mask_85_NA = mask_85[:,:,la1:la2,lo1:lo2]
density_85_NA = np.ma.masked_where(mask_85_NA[0,0,...], density[la1:la2,lo1:lo2])
npp_85_std_NA = np.ma.masked_where(mask_85_NA, npp_85_std[:,:,la1:la2,lo1:lo2])
npp_85_mean_NA = np.ma.masked_where(mask_85_NA, npp_85_mean[:,:,la1:la2,lo1:lo2])
npp_85_mean_s2n_NA = np.ma.masked_where(mask_85_NA, npp_85_mean_s2n[:,:,la1:la2,lo1:lo2])

# Benguela
la1 = 50; la2 = 90
lo1 = 330; lo2 = -1
print("longitude =",lon[0,lo1],lon[0,lo2])
print("latitude =",lat[la1,0],lat[la2,0])
mask_45_BE = mask_45[:,:,la1:la2,lo1:lo2]
density_45_BE = np.ma.masked_where(mask_45_BE[0,0,...], density[la1:la2,lo1:lo2])
npp_45_std_BE = np.ma.masked_where(mask_45_BE, npp_45_std[:,:,la1:la2,lo1:lo2])
npp_45_mean_BE = np.ma.masked_where(mask_45_BE, npp_45_mean[:,:,la1:la2,lo1:lo2])
npp_45_mean_s2n_BE = np.ma.masked_where(mask_45_BE, npp_45_mean_s2n[:,:,la1:la2,lo1:lo2])
mask_85_BE = mask_85[:,:,la1:la2,lo1:lo2]
density_85_BE = np.ma.masked_where(mask_85_BE[0,0,...], density[la1:la2,lo1:lo2])
npp_85_std_BE = np.ma.masked_where(mask_85_BE, npp_85_std[:,:,la1:la2,lo1:lo2])
npp_85_mean_BE = np.ma.masked_where(mask_85_BE, npp_85_mean[:,:,la1:la2,lo1:lo2])
npp_85_mean_s2n_BE = np.ma.masked_where(mask_85_BE, npp_85_mean_s2n[:,:,la1:la2,lo1:lo2])

# Amsterdam Island
la1 = 35; la2 = 70
lo1 = 30; lo2 = 80
print("longitude =",lon[0,lo1],lon[0,lo2])
print("latitude =",lat[la1,0],lat[la2,0])
mask_45_AI = mask_45[:,:,la1:la2,lo1:lo2]
density_45_AI = np.ma.masked_where(mask_45_AI[0,0,...], density[la1:la2,lo1:lo2])
npp_45_std_AI = np.ma.masked_where(mask_45_AI, npp_45_std[:,:,la1:la2,lo1:lo2])
npp_45_mean_AI = np.ma.masked_where(mask_45_AI, npp_45_mean[:,:,la1:la2,lo1:lo2])
npp_45_mean_s2n_AI = np.ma.masked_where(mask_45_AI, npp_45_mean_s2n[:,:,la1:la2,lo1:lo2])
mask_85_AI = mask_85[:,:,la1:la2,lo1:lo2]
density_85_AI = np.ma.masked_where(mask_85_AI[0,0,...], density[la1:la2,lo1:lo2])
npp_85_std_AI = np.ma.masked_where(mask_85_AI, npp_85_std[:,:,la1:la2,lo1:lo2])
npp_85_mean_AI = np.ma.masked_where(mask_85_AI, npp_85_mean[:,:,la1:la2,lo1:lo2])
npp_85_mean_s2n_AI = np.ma.masked_where(mask_85_AI, npp_85_mean_s2n[:,:,la1:la2,lo1:lo2])

# Southern Ocean
la1 = 0; la2 = 35
lo1 = 0; lo2 = 360
#print("longitude =",lon[0,lo1],lon[0,lo2])
print("latitude =",lat[la1,0],lat[la2,0])
mask_45_SO = mask_45[:,:,la1:la2,lo1:lo2]
density_45_SO = np.ma.masked_where(mask_45_SO[0,0,...], density[la1:la2,lo1:lo2])
npp_45_std_SO = np.ma.masked_where(mask_45_SO, npp_45_std[:,:,la1:la2,lo1:lo2])
npp_45_mean_SO = np.ma.masked_where(mask_45_SO, npp_45_mean[:,:,la1:la2,lo1:lo2])
npp_45_mean_s2n_SO = np.ma.masked_where(mask_45_SO, npp_45_mean_s2n[:,:,la1:la2,lo1:lo2])
mask_85_SO = mask_85[:,:,la1:la2,lo1:lo2]
density_85_SO = np.ma.masked_where(mask_85_SO[0,0,...], density[la1:la2,lo1:lo2])
npp_85_std_SO = np.ma.masked_where(mask_85_SO, npp_85_std[:,:,la1:la2,lo1:lo2])
npp_85_mean_SO = np.ma.masked_where(mask_85_SO, npp_85_mean[:,:,la1:la2,lo1:lo2])
npp_85_mean_s2n_SO = np.ma.masked_where(mask_85_SO, npp_85_mean_s2n[:,:,la1:la2,lo1:lo2])


plt.figure()
plt.pcolormesh(npp_85_mean_s2n_NA[0,0,...])
plt.figure()
plt.pcolormesh(npp_85_mean_s2n_BE[0,0,...])
plt.figure()
plt.pcolormesh(npp_85_mean_s2n_AI[0,0,...])
plt.figure()
plt.pcolormesh(npp_85_mean_s2n[0,3,...])


def weighting(data,dens):
    data = np.ma.masked_where(np.isnan(data), data)
    ww = dens / np.ma.sum(dens)
    tmp1 = data * ww
    data_w = np.ma.sum(np.ma.sum(tmp1, axis=3),axis=2)
    return data_w

npp_45_std_NAweighted = weighting(npp_45_std_NA, density_45_NA) 
npp_45_mean_NAweighted = weighting(npp_45_mean_NA, density_45_NA) 
npp_45_mean_s2n_NAweighted = weighting(npp_45_mean_s2n_NA, density_45_NA) 
npp_45_std_BEweighted = weighting(npp_45_std_BE, density_45_BE) 
npp_45_mean_BEweighted = weighting(npp_45_mean_BE, density_45_BE) 
npp_45_mean_s2n_BEweighted = weighting(npp_45_mean_s2n_BE, density_45_BE) 
npp_45_std_AIweighted = weighting(npp_45_std_AI, density_45_AI) 
npp_45_mean_AIweighted = weighting(npp_45_mean_AI, density_45_AI) 
npp_45_mean_s2n_AIweighted = weighting(npp_45_mean_s2n_AI, density_45_AI) 
npp_45_std_SOweighted = weighting(npp_45_std_SO, density_45_SO) 
npp_45_mean_SOweighted = weighting(npp_45_mean_SO, density_45_SO) 
npp_45_mean_s2n_SOweighted = weighting(npp_45_mean_s2n_SO, density_45_SO) 

npp_85_std_NAweighted = weighting(npp_85_std_NA, density_85_NA) 
npp_85_mean_NAweighted = weighting(npp_85_mean_NA, density_85_NA) 
npp_85_mean_s2n_NAweighted = weighting(npp_85_mean_s2n_NA, density_85_NA) 
npp_85_std_BEweighted = weighting(npp_85_std_BE, density_85_BE) 
npp_85_mean_BEweighted = weighting(npp_85_mean_BE, density_85_BE) 
npp_85_mean_s2n_BEweighted = weighting(npp_85_mean_s2n_BE, density_85_BE) 
npp_85_std_AIweighted = weighting(npp_85_std_AI, density_85_AI) 
npp_85_mean_AIweighted = weighting(npp_85_mean_AI, density_85_AI) 
npp_85_mean_s2n_AIweighted = weighting(npp_85_mean_s2n_AI, density_85_AI) 
npp_85_std_SOweighted = weighting(npp_85_std_SO, density_85_SO) 
npp_85_mean_SOweighted = weighting(npp_85_mean_SO, density_85_SO) 
npp_85_mean_s2n_SOweighted = weighting(npp_85_mean_s2n_SO, density_85_SO) 


#%% also compute the time of emergence based on no ability to adapt to conditions outside of that typical for 1850-1950

def time_of_emergence(data, window, poly):

    # remove mean of 1850-1950
    normalised = data - np.mean(data[0:100], axis=0)
    # find trends via decadal smoothing
    signal = np.zeros(np.shape(normalised))
    for m in np.arange(12):
        signal[:,m] = savitzky_golay(normalised[:,m], window, poly)
    # calculate standard deviation of 1850-1950 period
    std = np.std(normalised[0:100], axis=0)
    # signal:noise
    s2n = signal / std
    # set values > 1 when signal > noise, and vice versa
    emerge = np.float32(np.int64(np.abs(signal) > std))
    # eliminate transient emergences by decadal smoothing
    emerged = np.zeros(np.shape(normalised))
    for m in np.arange(12):
        emerged[:,m] = savitzky_golay(emerge[:,m], window, poly)
    return s2n, emerged


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial
    
    try:
        window_size = np.abs(np.int32(window_size))
        order = np.abs(np.int32(order))
    except(ValueError, msg):
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')


window = 11
poly = 1

npp_45_mean_toeS2N_NAweighted, npp_45_mean_toe_NAweighted = time_of_emergence(npp_45_mean_NAweighted, window, poly)
npp_45_mean_toeS2N_BEweighted, npp_45_mean_toe_BEweighted = time_of_emergence(npp_45_mean_BEweighted, window, poly)
npp_45_mean_toeS2N_AIweighted, npp_45_mean_toe_AIweighted = time_of_emergence(npp_45_mean_AIweighted, window, poly)
npp_45_mean_toeS2N_SOweighted, npp_45_mean_toe_SOweighted = time_of_emergence(npp_45_mean_SOweighted, window, poly)

npp_85_mean_toeS2N_NAweighted, npp_85_mean_toe_NAweighted = time_of_emergence(npp_85_mean_NAweighted, window, poly)
npp_85_mean_toeS2N_BEweighted, npp_85_mean_toe_BEweighted = time_of_emergence(npp_85_mean_BEweighted, window, poly)
npp_85_mean_toeS2N_AIweighted, npp_85_mean_toe_AIweighted = time_of_emergence(npp_85_mean_AIweighted, window, poly)
npp_85_mean_toeS2N_SOweighted, npp_85_mean_toe_SOweighted = time_of_emergence(npp_85_mean_SOweighted, window, poly)


#%% normalise the npp trends relative to the 1850-1950 conditions for plotting 

npp_45_norm_NAweighted = (npp_45_mean_NAweighted - np.average(npp_45_mean_NAweighted[0:100,:], axis=0)) / np.average(npp_45_mean_NAweighted[0:100,:], axis=0) * 100
npp_45_norm_BEweighted = (npp_45_mean_BEweighted - np.average(npp_45_mean_BEweighted[0:100,:], axis=0)) / np.average(npp_45_mean_BEweighted[0:100,:], axis=0) * 100
npp_45_norm_AIweighted = (npp_45_mean_AIweighted - np.average(npp_45_mean_AIweighted[0:100,:], axis=0)) / np.average(npp_45_mean_AIweighted[0:100,:], axis=0) * 100
npp_45_norm_SOweighted = (npp_45_mean_SOweighted - np.average(npp_45_mean_SOweighted[0:100,:], axis=0)) / np.average(npp_45_mean_SOweighted[0:100,:], axis=0) * 100

npp_85_norm_NAweighted = (npp_85_mean_NAweighted - np.average(npp_85_mean_NAweighted[0:100,:], axis=0)) / np.average(npp_85_mean_NAweighted[0:100,:], axis=0) * 100
npp_85_norm_BEweighted = (npp_85_mean_BEweighted - np.average(npp_85_mean_BEweighted[0:100,:], axis=0)) / np.average(npp_85_mean_BEweighted[0:100,:], axis=0) * 100
npp_85_norm_AIweighted = (npp_85_mean_AIweighted - np.average(npp_85_mean_AIweighted[0:100,:], axis=0)) / np.average(npp_85_mean_AIweighted[0:100,:], axis=0) * 100
npp_85_norm_SOweighted = (npp_85_mean_SOweighted - np.average(npp_85_mean_SOweighted[0:100,:], axis=0)) / np.average(npp_85_mean_SOweighted[0:100,:], axis=0) * 100


#%% two panel figure

fstic = 13
fslab = 15
alf1 = [0.7,0.7,0.7,0.7]
alf2 = [0.2,0.2,0.2,0.2]
ls = ['-','-','-','-']
lw = [1.5,1.5,1.5,1.5]
labs = ['North Atlantic (Jul-Aug)', 
        'Benguela Upwelling (Aug-Oct)',
        'Amsterdam Island (Nov-Dec)',
        'Southern Ocean (Jan-Mar)']
cols = ['k', 'royalblue', 'firebrick', 'goldenrod']

fig = plt.figure(figsize=(14,7))
gs = GridSpec(2,1)

ax1 = plt.subplot(gs[0])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.tick_params(labelsize=fstic)

plt.plot(years2100, np.average(npp_45_mean_toeS2N_NAweighted[:,6:8],axis=1), color=cols[0], linewidth=lw[0], alpha=alf1[0], linestyle='--')
plt.plot(years2100, np.average(npp_45_mean_toeS2N_BEweighted[:,7:10],axis=1), color=cols[1], linewidth=lw[1], alpha=alf1[1], linestyle='--')
plt.plot(years2100, np.average(npp_45_mean_toeS2N_AIweighted[:,10:12],axis=1), color=cols[2], linewidth=lw[2], alpha=alf1[2], linestyle='--')
plt.plot(years2100, np.average(npp_45_mean_toeS2N_SOweighted[:,0:4],axis=1), color=cols[3], linewidth=lw[3], alpha=alf1[3], linestyle='--')
plt.plot(years2100, np.average(npp_85_mean_toeS2N_NAweighted[:,6:8],axis=1), color=cols[0], linewidth=lw[0], alpha=alf1[0], linestyle=ls[0], label=labs[0])
plt.plot(years2100, np.average(npp_85_mean_toeS2N_BEweighted[:,7:10],axis=1), color=cols[1], linewidth=lw[1], alpha=alf1[1], linestyle=ls[1], label=labs[1])
plt.plot(years2100, np.average(npp_85_mean_toeS2N_AIweighted[:,10:12],axis=1), color=cols[2], linewidth=lw[2], alpha=alf1[2], linestyle=ls[2], label=labs[2])
plt.plot(years2100, np.average(npp_85_mean_toeS2N_SOweighted[:,0:4],axis=1), color=cols[3], linewidth=lw[3], alpha=alf1[3], linestyle=ls[3], label=labs[3])

plt.fill_between(years2100, -1, 1, facecolor='grey', alpha=0.25)
#plt.plot((1850,2100),(0,0),'k:', alpha=alf1[0])
#plt.plot((1850,2100),(-1,-1),'k--', alpha=alf1[0])
#plt.plot((1850,2100),(1,1),'k--', alpha=alf1[0])

plt.xlim(1850,2100)
plt.ylabel('Signal : Noise\n(no adaptation)', fontsize=fslab)


ax2 = plt.subplot(gs[1])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.tick_params(labelsize=fstic)

plt.plot(years2100[0:220], np.average(npp_45_mean_s2n_NAweighted[0:220,6:8],axis=1), color=cols[0], linewidth=lw[0], alpha=alf1[0], linestyle='--')
plt.plot(years2100[0:220], np.average(npp_45_mean_s2n_BEweighted[0:220,7:10],axis=1), color=cols[1], linewidth=lw[1], alpha=alf1[1], linestyle='--')
plt.plot(years2100[0:220], np.average(npp_45_mean_s2n_AIweighted[0:220,10:12],axis=1), color=cols[2], linewidth=lw[2], alpha=alf1[2], linestyle='--')
plt.plot(years2100[0:220], np.average(npp_45_mean_s2n_SOweighted[0:220,0:4],axis=1), color=cols[3], linewidth=lw[3], alpha=alf1[3], linestyle='--')
plt.plot(years2100[0:220], np.average(npp_85_mean_s2n_NAweighted[0:220,6:8],axis=1), color=cols[0], linewidth=lw[0], alpha=alf1[0], linestyle=ls[0])
plt.plot(years2100[0:220], np.average(npp_85_mean_s2n_BEweighted[0:220,7:10],axis=1), color=cols[1], linewidth=lw[1], alpha=alf1[1], linestyle=ls[1])
plt.plot(years2100[0:220], np.average(npp_85_mean_s2n_AIweighted[0:220,10:12],axis=1), color=cols[2], linewidth=lw[2], alpha=alf1[2], linestyle=ls[2])
plt.plot(years2100[0:220], np.average(npp_85_mean_s2n_SOweighted[0:220,0:4],axis=1), color=cols[3], linewidth=lw[3], alpha=alf1[3], linestyle=ls[3])


'''
window = 11; poly = 1
lwb = 2.0; alfb = 0.8
plt.plot(years2100, savitzky_golay(np.average(s2n_npp_NA_weighted[6:8,:],axis=0), window, order=poly), color=cols[0], label=labs[0], linewidth=lwb, alpha=alfb, linestyle=ls[0])
plt.plot(years2100, savitzky_golay(np.average(s2n_npp_BE_weighted[7:10,:],axis=0), window, order=poly), color=cols[1], label=labs[1], linewidth=lwb, alpha=alfb, linestyle=ls[1])
plt.plot(years2100, savitzky_golay(np.average(s2n_npp_AI_weighted[10:12,:],axis=0), window, order=poly), color=cols[2], label=labs[2], linewidth=lwb, alpha=alfb, linestyle=ls[2])
plt.plot(years2100, savitzky_golay(np.average(s2n_npp_SO_weighted[0:4,:],axis=0), window, order=poly), color=cols[3], label=labs[3], linewidth=lwb, alpha=alfb, linestyle=ls[3])
'''

plt.fill_between(years2100, -1, 1, facecolor='grey', alpha=0.25)
#plt.plot((1850,2100),(0,0),'k:', alpha=alf1[0])
#plt.plot((1850,2100),(-1,-1),'k--', alpha=alf1[0])
#plt.plot((1850,2100),(1,1),'k--', alpha=alf1[0])

plt.ylabel('Signal : Noise\n(adaptation)', fontsize=fslab)
plt.xlabel('year', fontsize=fslab)
plt.xlim(1850,2100)
plt.ylim(-5,5)


xx = 0.025; yy=0.95
plt.text(xx,yy,'a', fontweight='bold', va='center', ha='center', transform=ax1.transAxes, fontsize=fslab+2)
plt.text(xx,yy,'b', fontweight='bold', va='center', ha='center', transform=ax2.transAxes, fontsize=fslab+2)


# custom legend
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color=cols[0], lw=lw[0], alpha=alf1[0]),
                Line2D([0], [0], color=cols[1], lw=lw[1], alpha=alf1[1]),
                Line2D([0], [0], color=cols[2], lw=lw[2], alpha=alf1[2]),
                Line2D([0], [0], color=cols[3], lw=lw[3], alpha=alf1[3])]
fig.legend(custom_lines, labs, frameon=False, loc='upper center', ncol=2, bbox_to_anchor=(0.5,0.96), fontsize=fstic)

custom_lines = [Line2D([0], [0], color=cols[0], lw=lw[0], alpha=alf1[0], linestyle='--'),
                Line2D([0], [0], color=cols[0], lw=lw[0], alpha=alf1[0], linestyle='-')]
fig.legend(custom_lines, ['SSP2-4.5', 'SSP5-8.5'], frameon=False, loc='upper center', ncol=1, bbox_to_anchor=(0.25,0.7), fontsize=fstic)


print("North Atlantic mean change 2081-2100 (SSP2) =",np.average(npp_45_norm_NAweighted[231:251,6:8]))
print("Bengurla mean change 2081-2100 (SSP2) =",np.average(npp_45_norm_BEweighted[231:251,7:10]))
print("Amsterdam Island mean change 2081-2100 (SSP2) =",np.average(npp_45_norm_AIweighted[231:251,10:12]))
print("Southern Ocean mean change 2081-2100 (SSP2) =",np.average(npp_45_norm_SOweighted[231:251,0:4]))

print("North Atlantic mean change 2081-2100 (SSP5) =",np.average(npp_85_norm_NAweighted[231:251,6:8]))
print("Bengurla mean change 2081-2100 (SSP5) =",np.average(npp_85_norm_BEweighted[231:251,7:10]))
print("Amsterdam Island mean change 2081-2100 (SSP5) =",np.average(npp_85_norm_AIweighted[231:251,10:12]))
print("Southern Ocean mean change 2081-2100 (SSP5) =",np.average(npp_85_norm_SOweighted[231:251,0:4]))


### find years when the emergence occurs
print("ToE for NA (SSP2-4.5; non-conservative) =",years2100[np.abs(np.average(npp_45_mean_toeS2N_NAweighted[:,6:8],axis=1))>1.0][1])
print("ToE for BE (SSP2-4.5; non-conservative) =",years2100[np.abs(np.average(npp_45_mean_toeS2N_BEweighted[:,7:10],axis=1))>1.0][0])
print("ToE for AI (SSP2-4.5; non-conservative) =",years2100[np.abs(np.average(npp_45_mean_toeS2N_AIweighted[:,10:12],axis=1))>1.0][0])
print("ToE for SO (SSP2-4.5; non-conservative) =",years2100[np.abs(np.average(npp_45_mean_toeS2N_SOweighted[:,0:4],axis=1))>1.0][0])
print("ToE for NA (SSP5-8.5; non-conservative) =",years2100[np.abs(np.average(npp_85_mean_toeS2N_NAweighted[:,6:8],axis=1))>1.0][1])
print("ToE for BE (SSP5-8.5; non-conservative) =",years2100[np.abs(np.average(npp_85_mean_toeS2N_BEweighted[:,7:10],axis=1))>1.0][0])
print("ToE for AI (SSP5-8.5; non-conservative) =",years2100[np.abs(np.average(npp_85_mean_toeS2N_AIweighted[:,10:12],axis=1))>1.0][0])
print("ToE for SO (SSP5-8.5; non-conservative) =",years2100[np.abs(np.average(npp_85_mean_toeS2N_SOweighted[:,0:4],axis=1))>1.0][0])


print("ToE for NA (SSP2-4.5; conservative) =",years2100[0:220][np.abs(np.average(npp_45_mean_s2n_NAweighted[0:220,6:8],axis=1))>1.0][0])
print("ToE for NA (SSP5-8.5; conservative) =",years2100[0:220][np.abs(np.average(npp_85_mean_s2n_NAweighted[0:220,6:8],axis=1))>1.0][0])

# when signal is halved
print("ToE for NA (halved signal; SSP2-4.5; conservative) =",years2100[0:220][np.abs(np.average(npp_45_mean_s2n_NAweighted[0:220,6:8],axis=1))>2.0][0])
print("ToE for NA (halved signal; SSP5-8.5; conservative) =",years2100[0:220][np.abs(np.average(npp_85_mean_s2n_NAweighted[0:220,6:8],axis=1))>2.0][0])
print("ToE for NA (halved signal; SSP2-4.5; non-conservative) =",years2100[np.abs(np.average(npp_45_mean_toeS2N_NAweighted[:,6:8],axis=1))>2.0][0])
print("ToE for NA (halved signal; SSP5-8.5; non-conservative) =",years2100[np.abs(np.average(npp_85_mean_toeS2N_NAweighted[:,6:8],axis=1))>2.0][0])


#%%

#os.chdir('C://Users/pearseb/Dropbox/PostDoc/collaborations/CMIP6 hackathon 2021/figures')
os.chdir('/Users/pbuchanan/Dropbox/PostDoc/collaborations/CMIP6 hackathon 2021/figures')
fig.savefig('fig3_signaltonoise.png', dpi=300, bbox_inches='tight')
