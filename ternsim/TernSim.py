#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script tracks particles departing from a source location and flying with
an airspeed towards a destination location.
@author: Noam Vogt-Vincent
"""
import ternmethods as tm
import numpy as np
import os
from parcels import (FieldSet, ParticleSet, JITParticle, ErrorCode, Variable)
from netCDF4 import Dataset, num2date
from datetime import timedelta, datetime
from sys import argv

##############################################################################
# DIRECTORIES & PARAMETERS                                                   #
##############################################################################

param = {'model_name'        : [],
         'scenarios'         : ['SSP245',
                                'SSP585',
                                ],
         'mode'              : 'traj',

         'release_start_day' : 80,
         'release_end_day'   : 100,

         'number_of_releases': 3,
         'terns_per_release' : 2000,
         'release_lat'       : -70.,
         'release_lon_range' : [-60., -40.],       # [min, max]
         'target_lat'        : 65.,                # Where terns fly to
         'target_lon'        : -20.,               # Where terns fly to
         'sink_lat'          : 60.,                # Where terns are removed
         'airspeed'          : 10.,                # (m s-1)
         'fly_frac'          : 0.6,                # Fraction of day in flight
         'parcels_dt'        : timedelta(hours=1), # Parcels solver dt
         'out_dt'            : timedelta(hours=2),# Only used if mode == traj

         'var_name'          : ['uas', 'vas'],     # [zonal, meridional]
         'coordinate_name'   : ['lon', 'lat'],     # [lon, lat]

         'debug'             : False,               # Toggle to skip simulations
         'first_sim'         : 1                   # Only used if debug == True
         }

# Ask for model name if not given
try:
    param['model_name'] = argv[1]
except:
    param['model_name'] = input('Enter model name...')

dirs  = {'script': os.path.dirname(os.path.realpath(__file__)),
         'model': os.path.dirname(os.path.realpath(__file__)) + '/MODEL_DATA/' + param['model_name'] + '/',
         'traj': os.path.dirname(os.path.realpath(__file__)) + '/TRAJ_DATA/'}

# MODE NOTES:
# 'traj' : Full trajectory is recorded
# 'time' : Only time to reach destination is recorded

##############################################################################
# PROCESS INPUTS                                                             #
##############################################################################

print('Processing inputs...')
print('')

# CHECK MODE IS VALID
if not param['mode'] in ['traj', 'time']:
    raise NotImplementedError('Mode not understood!')

# FIND FILE NAMES
param['n_scen'] = len(param['scenarios'])

fh = {'u_hist' : (dirs['model'] + param['var_name'][0] + '_' +
                  param['model_name'] + '_HISTORICAL.nc'),
      'v_hist' : (dirs['model'] + param['var_name'][1] + '_' +
                  param['model_name'] + '_HISTORICAL.nc'),
      'u_scen' : [],
      'v_scen' : []}

for i in range(param['n_scen']):
    fh['u_scen'].append(dirs['model'] + param['var_name'][0] + '_' +
                        param['model_name'] + '_' + param['scenarios'][i] +
                        '.nc')
    fh['v_scen'].append(dirs['model'] + param['var_name'][1] + '_' +
                        param['model_name'] + '_' + param['scenarios'][i] +
                        '.nc')

fh['traj_hist'] = dirs['traj'] + param['model_name'] + '_HISTORICAL_TRAJ.nc'
fh['traj_scen'] = []

for i in range(param['n_scen']):
    fh['traj_scen'].append(dirs['traj'] + param['model_name'] + '_' +
                           param['scenarios'][i] + '_TRAJ.nc')


# GENERATE INITIAL POSITIONS
release = {'lon' : {'basis': np.linspace(param['release_lon_range'][0],
                                         param['release_lon_range'][1],
                                         num = param['terns_per_release'])},
           'lat' : {'basis': np.linspace(param['release_lat'],
                                         param['release_lat'],
                                         num = param['terns_per_release'])},}

# FIND START AND END TIMES
with Dataset(fh['u_hist'], mode='r') as nc:
    param['calendar'] = nc.variables['time'].calendar

    param['Ystart'] = {'hist' : num2date(nc.variables['time'][0],
                                         nc.variables['time'].units,
                                         calendar=param['calendar']).year}

    param['Yend']   = {'hist' : num2date(nc.variables['time'][-1],
                                         nc.variables['time'].units,
                                         calendar=param['calendar']).year}

    param['time_offset'] = (num2date(nc.variables['time'][0],
                                     nc.variables['time'].units) -
                            datetime(year  = param['Ystart']['hist'],
                                     month = 1,
                                     day   = 1)).total_seconds()

    param['endtime'] = {'hist' : (num2date(nc.variables['time'][-1],
                                           nc.variables['time'].units,
                                           calendar=param['calendar']) -
                                  num2date(nc.variables['time'][0],
                                           nc.variables['time'].units,
                                           calendar=param['calendar'])).total_seconds()}
    param['endtime']['hist'] -= 4*param['time_offset']

with Dataset(fh['u_scen'][0], mode='r') as nc:
    param['Ystart']['scen'] = num2date(nc.variables['time'][0],
                                       nc.variables['time'].units,
                                       calendar=param['calendar']).year

    param['Yend']['scen'] = num2date(nc.variables['time'][-1],
                                       nc.variables['time'].units,
                                       calendar=param['calendar']).year

    param['endtime']['scen'] = (num2date(nc.variables['time'][-1],
                                         nc.variables['time'].units,
                                         calendar=param['calendar']) -
                                num2date(nc.variables['time'][0],
                                         nc.variables['time'].units,
                                         calendar=param['calendar'])).total_seconds()
    param['endtime']['scen'] -= 4*param['time_offset']

release = tm.prepare_release(release, param)

# GENERATE THE FLYING VECTOR FIELD
fh['fly_field'] = dirs['model'] + 'fly_field.nc'
fly_field = tm.genTargetField(param['target_lon'],
                              param['target_lat'],
                              param['airspeed'],
                              fh)

##############################################################################
# TERN KERNELS                                                               #
##############################################################################

class ArcticTern_TRAJ(JITParticle):
    flight_time = Variable('flight_time', dtype=np.float32, initial=0.,
                           to_write=False)
    time_of_day = Variable('time_of_day', dtype=np.float32, initial=0.,
                           to_write=False)

class ArcticTern_TIME(JITParticle):
    flight_time = Variable('flight_time', dtype=np.float32, initial=0.,
                           to_write=True)
    time_of_day = Variable('time_of_day', dtype=np.float32, initial=0.,
                           to_write=False)
    # start_time  = Variable('start_time', dtype=np.int32, initial=0,
    #                        to_write=True)

def TimeOut(particle, fieldset, time):
    # Delete particle after excessive time (e.g. if particle gets stuck)
    if particle.flight_time >= 7776000: # 60 days
        particle.delete()

def DeleteParticle(particle, fieldset, time):
    #  Recovery kernel to delete a particle if it leaves the domain
    particle.delete()

def TernTools(particle, fieldset, time):
    # Remove the tern if they reach the end latitude
    if particle.lat > fieldset.sink_lat:
        particle.delete()

    # Update the particle age
    particle.flight_time += particle.dt
    particle.time_of_day += particle.dt
    if particle.time_of_day >= 86400.:
        particle.time_of_day = 0

def AdvectionRK4Tern(particle, fieldset, time):
    # Modified RK4 kernel to only activate during active hours
    # Modified from the parcels.kernels.advection module
    if particle.time_of_day < fieldset.night_time:
        (u1, v1) = fieldset.UV[particle]
        lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
        (u2, v2) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat1, lon1, particle]
        lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
        (u3, v3) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat2, lon2, particle]
        lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
        (u4, v4) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3, particle]
        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt


##############################################################################
# RELEASE THE TERNS                                                          #
##############################################################################

cs = {'time': ('time', 1),
      'lat' : ('lat', 64),
      'lon' : ('lon', 64)}

if not param['debug']:
    param['first_sim'] = 0  # Set first simulation to default of debug = False

for i in range(param['first_sim'], param['n_scen'] + 1):
    # Create fieldset
    if i == 0:
        print('Setting up historical simulation...')
        filenames = {'U': fh['u_hist'],
                     'V': fh['v_hist']}
    else:
        print('Setting up ' + param['scenarios'][i-1] + ' simulation...')
        filenames = {'U': fh['u_scen'][i-1],
                     'V': fh['v_scen'][i-1]}

    print('')

    variables = {'U': param['var_name'][0],
                 'V': param['var_name'][1]}

    dimensions = {'U': {'lon': param['coordinate_name'][0],
                        'lat': param['coordinate_name'][1],
                        'time': 'time'},
                  'V': {'lon': param['coordinate_name'][0],
                        'lat': param['coordinate_name'][1],
                        'time': 'time'}}

    wind_fieldset = FieldSet.from_netcdf(filenames,
                                         variables,
                                         dimensions)


    flight_fieldset = FieldSet.from_data(data = {'U': fly_field['u'],
                                                  'V': fly_field['v']},
                                          dimensions = {'lon': fly_field['lon'],
                                                        'lat': fly_field['lat']},
                                          allow_time_extrapolation=True)

    fieldset = FieldSet(U=wind_fieldset.U+flight_fieldset.U,
                        V=wind_fieldset.V+flight_fieldset.V)

    # Add constants
    fieldset.add_constant('sink_lat', param['sink_lat'])
    fieldset.add_constant('night_time', 86400.*param['fly_frac'])

    # Create particleset (dpeending on mode and hist/scen)
    if i == 0:
        if param['mode'] == 'traj':
            pset = ParticleSet.from_list(fieldset=fieldset,
                                         pclass=ArcticTern_TRAJ,
                                         lon = release['lon']['hist'],
                                         lat = release['lat']['hist'],
                                         time = release['time']['hist'])
        else:
            pset = ParticleSet.from_list(fieldset=fieldset,
                                         pclass=ArcticTern_TIME,
                                         lon = release['lon']['hist'],
                                         lat = release['lat']['hist'],
                                         time = release['time']['hist'])
        traj_fh = fh['traj_hist']
    else:
        if param['mode'] == 'traj':
            pset = ParticleSet.from_list(fieldset=fieldset,
                                         pclass=ArcticTern_TRAJ,
                                         lon = release['lon']['scen'],
                                         lat = release['lat']['scen'],
                                         time = release['time']['scen'])
        else:
            pset = ParticleSet.from_list(fieldset=fieldset,
                                         pclass=ArcticTern_TIME,
                                         lon = release['lon']['scen'],
                                         lat = release['lat']['scen'],
                                         time = release['time']['scen'])
        traj_fh = fh['traj_scen'][i-1]

    if param['mode'] == 'traj':
        traj_file = pset.ParticleFile(name=traj_fh,
                                      outputdt=param['out_dt'],
                                      )

    elif param['mode'] == 'time':
        traj_file = pset.ParticleFile(name=traj_fh,
                                      write_ondelete=True)
    else:
        raise NotImplementedError('Mode not understood!')

    # Run the simulation
    print('Set-up complete!')
    print('Releasing the birds!')

    Kernels = (pset.Kernel(AdvectionRK4Tern) +
               pset.Kernel(TernTools))

    if i == 0:
        pset.execute(Kernels,
                     endtime=param['endtime']['hist'],
                     dt = param['parcels_dt'],
                     recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
                     output_file=traj_file)
    else:
        pset.execute(Kernels,
                     endtime=param['endtime']['scen'],
                     dt = param['parcels_dt'],
                     recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
                     output_file=traj_file)

    print('')
    print('Simulation complete!')
    print('Exporting netcdf...')
    traj_file.export()

    print('Export complete!')
    print('')

    if i < param['n_scen']:
        print('Moving to the next simulation...')
        print('')
    else:
        print('Simulations complete!')
        print('The terns will miss you!')

