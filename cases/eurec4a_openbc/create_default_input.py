#
# MicroHH
# Copyright (c) 2011-2023 Chiel van Heerwaarden
# Copyright (c) 2011-2023 Thijs Heus
# Copyright (c) 2014-2023 Bart van Stratum
#
# This file is part of MicroHH
#
# MicroHH is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MicroHH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
#

import argparse

import netCDF4 as nc4
import xarray as xr
import pandas as pd
import numpy as np

import ls2d

# Make sure `microhhpy` is in the Python path
import microhhpy.io as io
import microhhpy.thermo as thermo

# Custom scripts from same directory.
from global_settings import float_type, work_path, cosmo_path, start_date, end_date
import helpers as hlp

# Grid (horizontal/vertical) definition.
from domain_definition import vgrid, outer_dom, inner_dom, zstart_buffer


"""
Parse command line arguments.
"""
parser = argparse.ArgumentParser(description='EUREC4A openBC MicroHH setup')
parser.add_argument(
        '-d', '--domain', required=True,
        help='Name of domain ("inner" or "outer")')
args = parser.parse_args()

if args.domain not in ('inner', 'outer'):
    raise Exception('Invalid domain choice.')

dates = pd.date_range(start_date, end_date, freq='h')
time_sec = np.array((dates-dates[0]).total_seconds()).astype(np.int32)


"""
Read grid info, and calculate spatial interpolation factors.
"""
domain = inner_dom if args.domain == 'inner' else outer_dom
dim_xy = (domain.jtot, domain.itot)

# Slice out LES domain with safe margin.
lon_slice = slice(domain.proj_pad.lon.min()-0.5, domain.proj_pad.lon.max()+0.5)
lat_slice = slice(domain.proj_pad.lat.min()-0.5, domain.proj_pad.lat.max()+0.5)

ds_2d, _ = hlp.read_cosmo(start_date, cosmo_path, lon_slice, lat_slice, read_3d=False)

# Calculate horizontal interpolation factors.
if_s = hlp.Calc_xy_interpolation_factors(
        ds_2d.rlon.values,
        ds_2d.rlat.values,
        domain.proj.lon,
        domain.proj.lat,
        domain.itot,
        domain.jtot,
        float_type)


"""
Process hourly surface COSMO fields: SST (thl_bot) and qsat(SST) (qt_bot).
Unit conversions (T -> thl etc.) are done in `read_cosmo()`.
"""
print('Interpolating COSMO surface fields.')

ds_cosmo = xr.open_dataset(f'{cosmo_path}/COSMO_CTRL_BC_nD_LES.nc')
ds_cosmo = ds_cosmo.sel(time=slice(start_date, end_date))

fld_s = np.empty(dim_xy, float_type)

for t, time in enumerate(time_sec):

    hlp.interpolate_cosmo(fld_s, ds_cosmo.thl_sbot[t,:,:].values, if_s, None, float_type)
    fld_s.tofile(f'{work_path}/thl_bot_in.{time:07d}')

    hlp.interpolate_cosmo(fld_s, ds_cosmo.qt_sbot[t,:,:].values, if_s, None, float_type)
    fld_s.tofile(f'{work_path}/qt_bot_in.{time:07d}')


"""
Create `eurec4a_input.nc` file.
"""
settings = {
    'central_lon' : -57.7,
    'central_lat' : 13.3,
    'start_date'  : start_date,
    'end_date'    : end_date,
    'case_name'   : 'eurec4a_openbc',
    'era5_path'   : '/home/scratch1/bart/LS2D_ERA5/'}

era = ls2d.Read_era5(settings)
era.calculate_forcings(n_av=6, method='2nd')

ds_time = era.get_les_input(vgrid.z)
ds_time = ds_time.sel(lay=slice(0,135), lev=slice(0,136))
ds_mean = ds_time.mean(dim='time')

def add_var(nc_group, name, dims, values):
    var = nc_group.createVariable(name, float_type, dims)
    var[:] = values

nc_file = nc4.Dataset(f'{work_path}/eurec4a_input.nc', mode='w', datamodel='NETCDF4', clobber=True)
nc_file.createDimension('z', vgrid.ktot)
add_var(nc_file, 'z', ('z'), vgrid.z)

# Initial profiles. Not really used, but needed for `init` phase. The resulting
# initial 3D fields are overwritten by the interpolated fields from COSMO.
nc_init = nc_file.createGroup('init')
add_var(nc_init, 'thl', ('z'), ds_mean.thl.values)
add_var(nc_init, 'qt', ('z'), ds_mean.qt.values)
add_var(nc_init, 'u', ('z'), ds_mean.u.values)
add_var(nc_init, 'v', ('z'), ds_mean.v.values)

# Time dependent input. For now, spatially constant.
nc_tdep = nc_file.createGroup('timedep')
nc_tdep.createDimension('time_surface', time_sec.size)
nc_tdep.createDimension('time_ls', time_sec.size)

add_var(nc_tdep, 'time_surface', ('time_surface'), time_sec)
add_var(nc_tdep, 'time_ls', ('time_ls'), time_sec)

add_var(nc_tdep, 'p_sbot', ('time_surface'), ds_time.ps.values)
add_var(nc_tdep, 'u_geo', ('time_ls', 'z'), ds_time.ug.values)
add_var(nc_tdep, 'v_geo', ('time_ls', 'z'), ds_time.vg.values)

# Radiation.
# TODO...

nc_file.close()


"""
Create basestate.
"""
bs = thermo.calc_moist_basestate(
    ds_mean.thl.values,
    ds_mean.qt.values,
    float(ds_mean.ps),
    vgrid.z,
    vgrid.zsize,
    float_type)

# Only save the density part for the dynamic core.
# A `_0` is appended to the name; this way the density created by `microhh init` can be overwritten with this one.
thermo.save_basestate_density(bs['rho'], bs['rhoh'], f'{work_path}/rhoref_0.0000000')


"""
Create `eurec4a.ini` from eurec4a.ini.base`, with correct settings.
"""
ini = io.read_ini('eurec4a.ini.base')

ini['master']['npx'] = 2
ini['master']['npy'] = 4

ini['grid']['itot'] = domain.itot
ini['grid']['jtot'] = domain.jtot
ini['grid']['ktot'] = vgrid.ktot

ini['grid']['xsize'] = domain.xsize
ini['grid']['ysize'] = domain.ysize
ini['grid']['zsize'] = vgrid.zsize

ini['buffer']['zstart'] = zstart_buffer
ini['time']['endtime'] = time_sec[-1]
ini['force']['fc'] = ds_time.attrs['fc']
ini['boundary_lateral']['n_sponge'] = domain.n_sponge

ini['cross']['xz'] = domain.ysize/2
ini['cross']['yz'] = domain.xsize/2

io.check_ini(ini)
io.save_ini(ini, f'{work_path}/eurec4a.ini')
