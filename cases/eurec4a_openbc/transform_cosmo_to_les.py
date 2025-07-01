from datetime import datetime, timedelta

import xarray as xr
import pandas as pd
import numpy as np

# Make sure `microhhpy` is in the Python path.
import microhhpy.thermo as thermo
import microhhpy.constants as cst

from global_settings import cosmo_path
from helpers import read_cosmo

"""
Settings
"""
lon_slice = slice(-60, -50)
lat_slice = slice(11, 18)

start = datetime(year=2020, month=2, day=1, hour=0)
end   = datetime(year=2020, month=2, day=12, hour=0)
dates = pd.date_range(start, end, freq='h')

# Read single file to get dimensions.
def base_name(date, type):
    if type == '2d':
        return f'{cosmo_path}/COSMO_CTRL_BC_2D/lffd{date.year:04d}{date.month:02d}{date.day:02d}{date.hour:02d}0000.nc'
    else:
        return f'{cosmo_path}/COSMO_CTRL_BC_3D/lffd{date.year:04d}{date.month:02d}{date.day:02d}{date.hour:02d}0000z.nc'

ds_2d = xr.open_dataset(base_name(start, '2d')).squeeze()
ds_3d = xr.open_dataset(base_name(start, '3d')).squeeze()

ds_2d = ds_2d.sel(rlon=lon_slice, rlat=lat_slice)
ds_3d = ds_3d.sel(rlon=lon_slice, rlat=lat_slice)

# Define new Dataset. The sliced area is sufficiently small to just pack everything
# (all time steps) into a single NetCDF file. Easy for lazily reading whatever you need.
dim3d = (dates.size, ds_3d.altitude.size, ds_2d.rlat.size, ds_2d.rlon.size)
dim2d = (dates.size, ds_2d.rlat.size, ds_2d.rlon.size)

ds = xr.Dataset(
    {
      'p': (['time', 'z', 'rlat', 'rlon'], np.empty(dim3d, dtype=np.float32)),
      'u': (['time', 'z', 'rlat', 'rlon'], np.empty(dim3d, dtype=np.float32)),
      'v': (['time', 'z', 'rlat', 'rlon'], np.empty(dim3d, dtype=np.float32)),
      'w': (['time', 'z', 'rlat', 'rlon'], np.empty(dim3d, dtype=np.float32)),
      'thl': (['time', 'z', 'rlat', 'rlon'], np.empty(dim3d, dtype=np.float32)),
      'qt': (['time', 'z', 'rlat', 'rlon'], np.empty(dim3d, dtype=np.float32)),
      'qr': (['time', 'z', 'rlat', 'rlon'], np.empty(dim3d, dtype=np.float32)),
      'p_bot': (['time', 'rlat', 'rlon'], np.empty(dim2d, dtype=np.float32)),
      'thl_sbot': (['time', 'rlat', 'rlon'], np.empty(dim2d, dtype=np.float32)),
      'qt_sbot':  (['time', 'rlat', 'rlon'], np.empty(dim2d, dtype=np.float32))
    },

    coords={
        'time': dates,
        'z': ds_3d.altitude.values,
        'rlat': ds_2d.rlat.values,
        'rlon': ds_2d.rlon.values
    }
)

"""
Parse all times.
"""
for t,date in enumerate(dates):
    print(f'Parsing {date}...')

    ds_2d = xr.open_dataset(base_name(date, '2d')).squeeze()
    ds_3d = xr.open_dataset(base_name(date, '3d')).squeeze()

    ds_2d = ds_2d.sel(rlon=lon_slice, rlat=lat_slice)
    ds_3d = ds_3d.sel(rlon=lon_slice, rlat=lat_slice)

    exner_2d = thermo.exner(ds_2d['PS'].values)
    exner_3d = thermo.exner(ds_3d['P'].values)

    # Short-cuts.
    ps = ds_2d['PS'].values
    Ts = ds_2d['T_S'].values
    u = ds_3d['U'].values
    v = ds_3d['V'].values
    w = ds_3d['W'].values
    p = ds_3d['P'].values
    T = ds_3d['T'].values
    th = T / exner_3d
    qv = ds_3d['QV'].values
    ql = ds_3d['QC'].values
    qi = ds_3d['QI'].values
    qr = ds_3d['QR'].values
    qs = ds_3d['QS'].values

    # 2D fields.
    ds['p_bot'][t,:,:] = ps
    ds['qt_sbot'][t,:,:] = thermo.qsat_liq(ps, Ts)
    ds['thl_sbot'][t,:,:] = Ts / exner_2d

    # 3D fields.
    ds['u'][t,:,:,:] = u
    ds['v'][t,:,:,:] = v
    ds['w'][t,:,:,:] = w
    ds['p'][t,:,:,:] = p
    ds['thl'][t,:,:,:] = th - cst.Lv / (cst.cp * exner_3d) * ql - cst.Ls / (cst.cp * exner_3d) * qi
    ds['qt'][t,:,:,:] = qv + ql + qi
    ds['qr'][t,:,:,:] = qr + qs

ds.to_netcdf(f'{cosmo_path}/COSMO_CTRL_BC_nD_LES_XL.nc')