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

from pathlib import Path
import argparse
import shutil
import glob
#import sys
import os

import netCDF4 as nc4
import xarray as xr
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

#sys.path.append('/home/bart/meteo/models/LS2D')
import ls2d

# Make sure `microhhpy` is in the Python path
import microhhpy.io as io
import microhhpy.thermo as thermo
from microhhpy.spatial import calc_vertical_grid_2nd
from microhhpy.openbc import create_era5_input
from microhhpy.interpolate import regrid_les

# Custom scripts from same directory.
from global_settings import env, domains, gd, zstart_buffer, start_date, end_date, float_type, host_model
import helpers as hlp


def interp(x_out, x_in, y_in):
    """
    Linear interpolation with extrapolation using Scipy's `interp1d`.
    """
    f = interp1d(x_in, y_in, fill_value='extrapolate')
    return f(x_out)


"""
Parse command line arguments and switch domain definition.
"""
parser = argparse.ArgumentParser(description='EUREC4A openBC MicroHH setup')
parser.add_argument('-d', '--domain', required=True, type=int, help='Domain number (0,1,N)')
args = parser.parse_args()

if args.domain > len(domains)-1:
    raise Exception('Domain number out of range.')


"""
Domain switches.
"""
domain = domains[args.domain]

parent = domain.parent
child = domain.child


"""
Read ERA5 and CAMS with (LS)2D, for backgroud radiation, and/or domain initialisation.
"""
settings = {
    'central_lon' : -57.7,
    'central_lat' : 13.3,
    'start_date'  : start_date,
    'end_date'    : end_date,
    'case_name'   : 'eurec4a_xl',
    'era5_path'   : env["ls2d_era5_path"],
    'cams_path'   : env["ls2d_cams_path"],
    }

dates = pd.date_range(start_date, end_date, freq='h')
time_sec = np.array((dates-dates[0]).total_seconds()).astype(np.int32)

era = ls2d.Read_era5(settings)
era.calculate_forcings(n_av=6, method='2nd')

# Interpolate to both outer and inner vertical grid.
era_time = era.get_les_input(gd['z'])
era_time = era_time.sel(lay=slice(0,135), lev=slice(0,136))
era_mean = era_time.mean(dim='time')

# Read CAMS data.
cams_vars = {
        'eac4_ml': [
            'dust_aerosol_0.03-0.55um_mixing_ratio',
            'dust_aerosol_0.55-0.9um_mixing_ratio',
            'dust_aerosol_0.9-20um_mixing_ratio',
            'hydrophilic_black_carbon_aerosol_mixing_ratio',
            'hydrophilic_organic_matter_aerosol_mixing_ratio',
            'hydrophobic_black_carbon_aerosol_mixing_ratio',
            'hydrophobic_organic_matter_aerosol_mixing_ratio',
            'sea_salt_aerosol_0.03-0.5um_mixing_ratio',
            'sea_salt_aerosol_0.5-5um_mixing_ratio',
            'sea_salt_aerosol_5-20um_mixing_ratio',
            'specific_humidity',
            'sulphate_aerosol_mixing_ratio',
            'temperature'],
        'eac4_sfc': [
            'surface_pressure'],
        }

cams = ls2d.Read_cams(settings, cams_vars)
cams_time = cams.get_les_input(gd['z'], 6)

# Interpolate CAMS `lay` variables onto ERA5 pressures.
tmp = cams_time.copy()
drop_vars = ['lay']
for name, da in tmp.data_vars.items():
    if 'lay' in da.dims:
        drop_vars.append(name)
tmp = tmp.drop_vars(drop_vars)

# Target pressure:
p_lay = era_time.p_lay

for name, da in cams_time.data_vars.items():
    if 'lay' in da.dims:
        out = np.zeros_like(p_lay)
        for t in range(out.shape[0]):
            out[t,:] = interp(p_lay[t], cams_time.p_lay[t,:], da[t,:])
        tmp[name] = (('time', 'lay'), out)

# Copy back to original name.
cams_time = tmp

# Limit aerosols to >= eps.
# CAMS sometimes has small negative values (-1e-16), which causes RRTMGP to fail.
eps = 1e-16
for i in range(1, 12):
    cams_time[f'aermr{i:02d}'] = np.maximum(cams_time[f'aermr{i:02d}'], eps)
    cams_time[f'aermr{i:02d}_lay'] = np.maximum(cams_time[f'aermr{i:02d}_lay'], eps)

cams_mean = cams_time.mean(dim='time')


"""
Process hourly surface COSMO/ERA fields: SST (thl_bot) and qsat(SST) (qt_bot).
Unit conversions for COSMO (T -> thl etc.) are done in `read_cosmo()`.
"""

fld_s = np.empty((domain.jtot, domain.itot), float_type)
tmp_s = np.empty((domain.jtot, domain.itot), float_type)

if host_model == 'COSMO':
    print('Interpolating COSMO surface fields.')

    ds_cosmo = xr.open_dataset(f'{env["cosmo_path"]}/COSMO_CTRL_BC_nD_LES_XL.nc')
    ds_cosmo = ds_cosmo.sel(time=slice(start_date, end_date))

    # Calculate horizontal interpolation factors.
    if_s = hlp.Calc_xy_interpolation_factors(
            ds_cosmo.rlon.values,
            ds_cosmo.rlat.values,
            domain.proj.lon,
            domain.proj.lat,
            domain.itot,
            domain.jtot,
            float_type)

    for t, time in enumerate(time_sec):

        hlp.interpolate_cosmo(fld_s, ds_cosmo.thl_sbot[t,:,:].values, if_s, None, float_type)
        fld_s.tofile(f'{domain.work_dir}/thl_bot_in.{time:07d}')

        hlp.interpolate_cosmo(fld_s, ds_cosmo.qt_sbot[t,:,:].values, if_s, None, float_type)
        fld_s.tofile(f'{domain.work_dir}/qt_bot_in.{time:07d}')

elif host_model == 'ERA5':
    print('Interpolating ERA5 surface fields.')

    # Calculate horizontal interpolation factors.
    if_s = hlp.Calc_xy_interpolation_factors(
            era.lons.data,
            era.lats.data,
            domain.proj.lon,
            domain.proj.lat,
            domain.itot,
            domain.jtot,
            float_type)

    for t, time in enumerate(time_sec):

        # Interpolate surface pressure to LES.
        hlp.interpolate_cosmo(tmp_s, era.ps[t,:,:], if_s, None, float_type)

        # Sea surface temperature.
        hlp.interpolate_cosmo(fld_s, era.sst[t,:,:].data, if_s, None, float_type)

        # For moisture, assume `q_bot = qsat(sst)`.
        qt_s = thermo.qsat(tmp_s, fld_s)
        qt_s.tofile(f'{domain.work_dir}/qt_bot_in.{time:07d}')

        # SST from absolute to potential temperature.
        fld_s /= thermo.exner(tmp_s)
        fld_s.tofile(f'{domain.work_dir}/thl_bot_in.{time:07d}')


"""
Create `eurec4a_input.nc` file.
"""
def add_nc_var(name, dims, nc, data):
    if name not in nc.variables:
        if dims is None:
            var = nc.createVariable(name, np.float64)
        else:
            var = nc.createVariable(name, np.float64, dims)
        var[:] = data

def add_nc_dim(name, size, nc):
    if name not in nc.dimensions:
        nc.createDimension(name, size)


nc_file = nc4.Dataset(f'{domain.work_dir}/eurec4a_input.nc', mode='w', datamodel='NETCDF4', clobber=True)
add_nc_dim('z', gd['ktot'], nc_file)
add_nc_var('z', ('z'), nc_file, gd['z'])

# Initial profiles. Not really used, but needed for `init` phase. The resulting
# initial 3D fields are overwritten by the interpolated fields from COSMO.
nc_init = nc_file.createGroup('init')
add_nc_var('thl', ('z'), nc_init, era_mean.thl.values)
add_nc_var('qt', ('z'), nc_init, era_mean.qt.values)
add_nc_var('u', ('z'), nc_init, era_mean.u.values)
add_nc_var('v', ('z'), nc_init, era_mean.v.values)

# Time dependent input. For now, spatially constant.
nc_tdep = nc_file.createGroup('timedep')
add_nc_dim('time_surface', time_sec.size, nc_tdep)
add_nc_dim('time_ls', time_sec.size, nc_tdep)

add_nc_var('time_surface', ('time_surface'), nc_tdep, time_sec)
add_nc_var('time_ls', ('time_ls'), nc_tdep, time_sec)

add_nc_var('p_sbot', ('time_surface'), nc_tdep, era_time.ps.values)
add_nc_var('u_geo', ('time_ls', 'z'), nc_tdep, era_time.ug.values)
add_nc_var('v_geo', ('time_ls', 'z'), nc_tdep, era_time.vg.values)

# Radiation.
nc_rad = nc_file.createGroup('radiation')
add_nc_dim('lay', era_mean.sizes['lay'], nc_rad)
add_nc_dim('lev', era_mean.sizes['lev'], nc_rad)

# Radiation variables on LES grid.
xm_air = 28.97; xm_h2o = 18.01528; eps = xm_h2o / xm_air
h2o = era_mean.qt / (eps - eps * era_mean.qt)
add_nc_var('h2o', ('z'), nc_init, h2o)
add_nc_var('o3',  ('z'), nc_init, era_mean.o3*1e-6)

# RFMIP background concentrations.
rfmip = {
    'co2': 0.00039754696655273437,
    'ch4': 1.8314709472656252e-06,
    'n2o': 3.269880065917969e-07,
    'n2': 0.781000018119812,
    'o2': 0.20900000631809235,
    'co': 1.199999957179898e-07,
    'ccl4': 8.306993103027344e-11,
    'cfc11': 2.330798645019531e-10,
    'cfc12': 5.205809936523438e-10,
    'hcfc22': 2.295420684814453e-10,
    'hfc143a': 1.525278091430664e-11,
    'hfc125': 1.5355008125305177e-11,
    'hfc23': 2.6890436172485352e-11,
    'hfc32': 8.336969375610351e-12,
    'hfc134a': 8.051573181152344e-11,
    'cf4': 8.109249114990234e-11,
    'no2': 0.0}

for group in (nc_init, nc_rad):
    for name, value in rfmip.items():
        add_nc_var(name, None, group, value)

# Radiation variables on radiation grid/levels:
add_nc_var('z_lay', ('lay'), nc_rad, era_mean.z_lay)
add_nc_var('z_lev', ('lev'), nc_rad, era_mean.z_lev)
add_nc_var('p_lay', ('lay'), nc_rad, era_mean.p_lay)
add_nc_var('p_lev', ('lev'), nc_rad, era_mean.p_lev)
add_nc_var('t_lay', ('lay'), nc_rad, era_mean.t_lay)
add_nc_var('t_lev', ('lev'), nc_rad, era_mean.t_lev)
add_nc_var('h2o',   ('lay'), nc_rad, era_mean.h2o_lay)
add_nc_var('o3',    ('lay'), nc_rad, era_mean.o3_lay*1e-6)

add_nc_var('aermr01', ('z'), nc_init, cams_mean.aermr01)
add_nc_var('aermr02', ('z'), nc_init, cams_mean.aermr02)
add_nc_var('aermr03', ('z'), nc_init, cams_mean.aermr03)
add_nc_var('aermr04', ('z'), nc_init, cams_mean.aermr04)
add_nc_var('aermr05', ('z'), nc_init, cams_mean.aermr05)
add_nc_var('aermr06', ('z'), nc_init, cams_mean.aermr06)
add_nc_var('aermr07', ('z'), nc_init, cams_mean.aermr07)
add_nc_var('aermr08', ('z'), nc_init, cams_mean.aermr08)
add_nc_var('aermr09', ('z'), nc_init, cams_mean.aermr09)
add_nc_var('aermr10', ('z'), nc_init, cams_mean.aermr10)
add_nc_var('aermr11', ('z'), nc_init, cams_mean.aermr11)

add_nc_var('aermr01', ('lay'), nc_rad, cams_mean.aermr01_lay)
add_nc_var('aermr02', ('lay'), nc_rad, cams_mean.aermr02_lay)
add_nc_var('aermr03', ('lay'), nc_rad, cams_mean.aermr03_lay)
add_nc_var('aermr04', ('lay'), nc_rad, cams_mean.aermr04_lay)
add_nc_var('aermr05', ('lay'), nc_rad, cams_mean.aermr05_lay)
add_nc_var('aermr06', ('lay'), nc_rad, cams_mean.aermr06_lay)
add_nc_var('aermr07', ('lay'), nc_rad, cams_mean.aermr07_lay)
add_nc_var('aermr08', ('lay'), nc_rad, cams_mean.aermr08_lay)
add_nc_var('aermr09', ('lay'), nc_rad, cams_mean.aermr09_lay)
add_nc_var('aermr10', ('lay'), nc_rad, cams_mean.aermr10_lay)
add_nc_var('aermr11', ('lay'), nc_rad, cams_mean.aermr11_lay)

# Time dep background radiation and aerosols.
add_nc_dim('time_rad', time_sec.size, nc_tdep)
add_nc_dim('lay', era_mean.z_lay.size, nc_tdep)
add_nc_dim('lev', era_mean.z_lev.size, nc_tdep)

add_nc_var('time_rad', ('time_rad'), nc_tdep, time_sec)

add_nc_var('z_lay',  ('time_rad', 'lay'), nc_tdep, era_time.z_lay[:,:])
add_nc_var('z_lev',  ('time_rad', 'lev'), nc_tdep, era_time.z_lev[:,:])
add_nc_var('p_lay',  ('time_rad', 'lay'), nc_tdep, era_time.p_lay[:,:])
add_nc_var('p_lev',  ('time_rad', 'lev'), nc_tdep, era_time.p_lev[:,:])
add_nc_var('t_lay',  ('time_rad', 'lay'), nc_tdep, era_time.t_lay[:,:])
add_nc_var('t_lev',  ('time_rad', 'lev'), nc_tdep, era_time.t_lev[:,:])
add_nc_var('o3_bg',  ('time_rad', 'lay'), nc_tdep, era_time.o3_lay[:,:]*1e-6)
add_nc_var('o3',     ('time_rad', 'z'  ), nc_tdep, era_time.o3[:,:]*1e-6)
add_nc_var('h2o_bg', ('time_rad', 'lay'), nc_tdep, era_time.h2o_lay[:,:])

add_nc_var('aermr01',  ('time_rad', 'z'),   nc_tdep, cams_time.aermr01[:,:])
add_nc_var('aermr02',  ('time_rad', 'z'),   nc_tdep, cams_time.aermr02[:,:])
add_nc_var('aermr03',  ('time_rad', 'z'),   nc_tdep, cams_time.aermr03[:,:])
add_nc_var('aermr04',  ('time_rad', 'z'),   nc_tdep, cams_time.aermr04[:,:])
add_nc_var('aermr05',  ('time_rad', 'z'),   nc_tdep, cams_time.aermr05[:,:])
add_nc_var('aermr06',  ('time_rad', 'z'),   nc_tdep, cams_time.aermr06[:,:])
add_nc_var('aermr07',  ('time_rad', 'z'),   nc_tdep, cams_time.aermr07[:,:])
add_nc_var('aermr08',  ('time_rad', 'z'),   nc_tdep, cams_time.aermr08[:,:])
add_nc_var('aermr09',  ('time_rad', 'z'),   nc_tdep, cams_time.aermr09[:,:])
add_nc_var('aermr10',  ('time_rad', 'z'),   nc_tdep, cams_time.aermr10[:,:])
add_nc_var('aermr11',  ('time_rad', 'z'),   nc_tdep, cams_time.aermr11[:,:])

add_nc_var('aermr01_bg', ('time_rad', 'lay'), nc_tdep, cams_time.aermr01_lay[:,:])
add_nc_var('aermr02_bg', ('time_rad', 'lay'), nc_tdep, cams_time.aermr02_lay[:,:])
add_nc_var('aermr03_bg', ('time_rad', 'lay'), nc_tdep, cams_time.aermr03_lay[:,:])
add_nc_var('aermr04_bg', ('time_rad', 'lay'), nc_tdep, cams_time.aermr04_lay[:,:])
add_nc_var('aermr05_bg', ('time_rad', 'lay'), nc_tdep, cams_time.aermr05_lay[:,:])
add_nc_var('aermr06_bg', ('time_rad', 'lay'), nc_tdep, cams_time.aermr06_lay[:,:])
add_nc_var('aermr07_bg', ('time_rad', 'lay'), nc_tdep, cams_time.aermr07_lay[:,:])
add_nc_var('aermr08_bg', ('time_rad', 'lay'), nc_tdep, cams_time.aermr08_lay[:,:])
add_nc_var('aermr09_bg', ('time_rad', 'lay'), nc_tdep, cams_time.aermr09_lay[:,:])
add_nc_var('aermr10_bg', ('time_rad', 'lay'), nc_tdep, cams_time.aermr10_lay[:,:])
add_nc_var('aermr11_bg', ('time_rad', 'lay'), nc_tdep, cams_time.aermr11_lay[:,:])

nc_file.close()


"""
Create basestate.
"""
bs = thermo.calc_moist_basestate(
    era_mean.thl.values,
    era_mean.qt.values,
    float(era_mean.ps),
    gd['z'],
    gd['zsize'],
    float_type)

# Only save the density part for the dynamic core.
# `_ext` is appended to the name; this way the density created by `microhh init` can be overwritten with this one.
thermo.save_basestate_density(bs['rho'], bs['rhoh'], f'{domain.work_dir}/rhoref_ext.0000000')


"""
Create 2D field with coriolis frequency.
"""
omega = 2*np.pi/86400
fc = 2 * omega * np.sin(np.deg2rad(domain.proj.lat))
fc.astype(float_type).tofile(f'{domain.work_dir}/fc.0000000')


"""
Create `eurec4a.ini` from eurec4a.ini.base`, with correct settings.
"""
ini = io.read_ini('eurec4a.ini.base')

# Check turbulence recycling
imax = domain.itot / domain.npx
jmax = domain.jtot / domain.npy

if ini['boundary_lateral']['recycle_offset'] + domain.n_sponge > imax + domain.n_ghost:
    raise Exception(f'Recycle offset x too large; max value = {imax + domain.n_ghost - domain.n_sponge}')
if ini['boundary_lateral']['recycle_offset'] + domain.n_sponge > jmax + domain.n_ghost:
    raise Exception(f'Recycle offset y too large; max value = {imax + domain.n_ghost - domain.n_sponge}')

ini['master']['npx'] = domain.npx
ini['master']['npy'] = domain.npy

ini['grid']['itot'] = domain.itot
ini['grid']['jtot'] = domain.jtot
ini['grid']['ktot'] = gd['ktot']

ini['grid']['xsize'] = domain.xsize
ini['grid']['ysize'] = domain.ysize
ini['grid']['zsize'] = gd['zsize']

ini['buffer']['zstart'] = zstart_buffer
ini['buffer']['loadfreq'] = 3600 if parent is None else 600

ini['time']['endtime'] = time_sec[-1]
ini['time']['datetime_utc'] = start_date.strftime('%Y-%m-%d %H:%M:%S')

ini['boundary_lateral']['loadfreq'] = domain.lbc_freq
ini['boundary_lateral']['n_sponge'] = domain.n_sponge
ini['boundary_lateral']['tau_sponge'] = domain.n_sponge * domain.dx / (10 * 5)

if parent is None:
    ini['boundary_lateral']['slist'] = ['thl', 'qt', 'qr']
    ini['boundary_lateral']['sw_recycle[east]'] = True
    ini['boundary_lateral']['sw_recycle[north]'] = True
    ini['radiation']['dt_rad'] = 300
else:
    ini['boundary_lateral']['slist'] = ['thl', 'qt', 'qr', 'nr']
    ini['radiation']['dt_rad'] = 60

if child is not None:
    ini['subdomain']['sw_subdomain'] = True

    ini['subdomain']['xstart'] = child.xstart_in_parent
    ini['subdomain']['ystart'] = child.ystart_in_parent
    ini['subdomain']['xend'] = child.xstart_in_parent + child.xsize
    ini['subdomain']['yend'] = child.ystart_in_parent + child.ysize

    ini['subdomain']['grid_ratio_ij'] = int(domain.dx / child.dx)
    ini['subdomain']['grid_ratio_k'] = 1
    ini['subdomain']['n_ghost'] = child.n_ghost
    ini['subdomain']['n_sponge'] = child.n_sponge

    ini['subdomain']['sw_save_wtop'] = True
    ini['subdomain']['sw_save_buffer'] = True

    ini['subdomain']['savetime_bcs'] = child.lbc_freq
    ini['subdomain']['savetime_buffer'] = 600
    ini['subdomain']['zstart_buffer'] = zstart_buffer

else:
    ini['subdomain']['sw_subdomain'] = False

# NOTE: 10, 100 not part of MIP. 10 = needed for 10 m wind, temp, etc.
ini['cross']['xy'] = [10, 400, 1000, 1500, 3000, 5000, gd['zsize']]
#ini['cross']['xz'] = domain.ysize / 2
#ini['cross']['yz'] = domain.xsize / 2

# Column stats for BCO.
x,y = domain.proj.to_xy(lon=-59.432, lat=13.165)
if x < 0 or x > domain.xsize or x < 0 or y > domain.ysize:
    print('Warning: BCO coordinates outside domain. Setting column to domain center')
    x = domain.xsize / 2
    y = domain.ysize / 2

ini['column']['coordinates[x]'] = x
ini['column']['coordinates[y]'] = y

# 3D dump, local near domain center. Only for 100 m domain.
if args.domain == 1:
    ini['dump']['swdump_sub'] = True
    ini['dump']['mpicoordx'] = [domain.npx//2-1, domain.npx//2]
    ini['dump']['mpicoordy'] = [domain.npy//2-1, domain.npy//2, domain.npy//2+1]
else:
    ini['dump']['swdump_sub'] = False


io.check_ini(ini)
io.save_ini(ini, f'{domain.work_dir}/eurec4a.ini')


"""
Link/copy required RRTMGP lookup tables.
"""
rrtmgp_path = f'{env["microhh_path"]}/rte-rrtmgp-cpp/'
rrtmgp_data_path = f'{env["microhh_path"]}/rte-rrtmgp-cpp/rrtmgp-data'

to_copy = [
        (f'{env["gpt_veerman_path"]}/rrtmgp-gas-lw-g056-cf2.nc', 'coefficients_lw.nc'),
        (f'{env["gpt_veerman_path"]}/rrtmgp-gas-sw-g049-cf2.nc', 'coefficients_sw.nc'),
        (f'{rrtmgp_data_path}/rrtmgp-clouds-lw.nc', 'cloud_coefficients_lw.nc'),
        (f'{rrtmgp_data_path}/rrtmgp-clouds-sw.nc', 'cloud_coefficients_sw.nc'),
        (f'{rrtmgp_path}/data/aerosol_optics.nc', 'aerosol_optics.nc')]

for f in to_copy:
    target = f'{domain.work_dir}/{f[1]}'
    if not os.path.exists(target):
        shutil.copy(f[0], target)


if parent is None:
    """
    Create input from COSMO 3D files.
    """

    if host_model == 'COSMO':

        cosmo = xr.open_dataset(f'{env["cosmo_path"]}/COSMO_CTRL_BC_nD_LES_XL.nc')
        cosmo = cosmo.sel(time=slice(start_date, end_date))

        fields_ls = {
            'u': cosmo.u[:,:,:,:],
            'v': cosmo.v[:,:,:,:],
            'w': cosmo.w[:,:,:,:],
            'thl': cosmo.thl[:,:,:,:],
            'qt': cosmo.qt[:,:,:,:],
            'qr': cosmo.qr[:,:,:,:],
        }

        lon_ls = cosmo.rlon.data
        lat_ls = cosmo.rlat.data

        p_ls = cosmo.p[:,:,:,:]
        z_ls = np.broadcast_to(cosmo.z.values[None,:,None,None], cosmo.thl.shape)
        time_ls = cosmo.time - cosmo.time[0]
        time_ls = (time_ls.values.astype(np.float32)/1e9).astype(np.int32)

    elif host_model == 'ERA5':

        fields_ls = {
            'u': era.u[:,:,:,:],
            'v': era.v[:,:,:,:],
            'w': era.wls[:,:,:,:],
            'thl': era.thl[:,:,:,:],
            'qt': era.qt[:,:,:,:],
            'qr': era.qr[:,:,:,:],
        }

        lon_ls = era.lons.data
        lat_ls = era.lats.data

        p_ls = era.p[:,:,:,:]
        z_ls = era.z[:,:,:,:]
        time_ls = era.time_sec


    # Standard dev. of Gaussian filter applied to interpolated fields (m).
    sigma_h = 1_000

    create_era5_input(
        fields_ls,
        lon_ls,
        lat_ls,
        z_ls,
        p_ls,
        time_ls,
        gd['z'],
        gd['zsize'],
        zstart_buffer,
        bs['rho'],
        bs['rhoh'],
        domain,
        sigma_h,
        perturb_size=3,
        perturb_amplitude={'thl': 0.1, 'qt': 0.1e-3},
        perturb_max_height=1000,
        clip_at_zero=['qt', 'qr'],
        name_suffix='ext',
        output_dir=domain.work_dir,
        ntasks=8,
        dtype=float_type)


else:
    """
    Create/copy input from parent domain.
    """
    # Regrid t=0 restart files.
    fields_3d = {
            'u': 0,
            'v': 0,
            'w': 0,
            'thl': 0,
            'qt': 0,
            'qr': 0,
            'nr': 0}

    # Regrid all 2D pressure @ TOD files.
    fields_2d = {
            'phydro_tod': '*'}

    regrid_les(
            fields_3d,
            fields_2d,
            parent.xsize,
            parent.ysize,
            gd['z'],
            gd['zh'],
            parent.itot,
            parent.jtot,
            domain.xsize,
            domain.ysize,
            gd['z'],
            gd['zh'],
            domain.itot,
            domain.jtot,
            domain.xstart_in_parent,
            domain.ystart_in_parent,
            parent.work_dir,
            domain.work_dir,
            float_type,
            name_suffix='ext')


    # Link boundary conditions from parent to child domain.
    # Symlinks can be created before parent simulation is finished,
    # but the parent domain must stay ahead of child domain
    # to ensure boundary condition files are available when needed.
    buffer_fields = ['u', 'v', 'w', 'thl', 'qt']
    lbc_fields = ['u', 'v', 'w', 'thl', 'qt', 'qr', 'nr']

    buffer_freq = 600

    for field in buffer_fields:
        for time in range(0, time_sec[-1], buffer_freq):
            src = Path(f'{parent.work_dir}/{field}_buffer_out.{time:07d}').resolve()
            dst = Path(f'{domain.work_dir}/{field}_buffer.{time:07d}').resolve()
            os.symlink(src, dst)
            
    for field in lbc_fields:
        for edge in ('west', 'north', 'east', 'south'):
            for time in range(0, time_sec[-1], domain.lbc_freq):
                src = Path(f'{parent.work_dir}/lbc_{field}_{edge}_out.{time:07d}').resolve()
                dst = Path(f'{domain.work_dir}/lbc_{field}_{edge}.{time:07d}').resolve()
                os.symlink(src, dst)

    for time in range(0, time_sec[-1], domain.lbc_freq):
        src = Path(f'{parent.work_dir}/w_top_out.{time:07d}').resolve()
        dst = Path(f'{domain.work_dir}/w_top.{time:07d}').resolve()
        os.symlink(src, dst)



    # Link boundary conditions from parent to child domain.
    # Only link the `_out` files, without `_out` they are LBCs used as input for the parent domain.
    #def link_files(src_pattern):
    #    print(f'Linking files {src_pattern} to {domain.work_dir}.')

    #    files = glob.glob(src_pattern)

    #    for f in files:
    #        src = Path(f).resolve()
    #        name = src.name.replace('_out', '')
    #        dst = Path(domain.work_dir) / name
    #        dst.symlink_to(src)

    #link_files(f'{parent.work_dir}/lbc_*_out.*')
    #link_files(f'{parent.work_dir}/w_top_out.*')
    #link_files(f'{parent.work_dir}/*_buffer_out.*')
