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
import os

import netCDF4 as nc4
import xarray as xr
import pandas as pd
import numpy as np

import ls2d

# Make sure `microhhpy` is in the Python path
import microhhpy.io as io
import microhhpy.thermo as thermo
from microhhpy.spatial import calc_vertical_grid_2nd
from microhhpy.openbc import create_era5_input
from microhhpy.interpolate import regrid_les

# Custom scripts from same directory.
from global_settings import outer_dom, inner_dom, start_date, end_date, float_type
from global_settings import work_path, cosmo_path, microhh_path, gpt_veerman_path, ls2d_era5_path
import helpers as hlp

# Vertical grid definition. Switch between domains is set in `global_settings.py`.
from domain_definition import vgrid, zstart_buffer


"""
Parse command line arguments and switch domain definition.
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


domain = inner_dom if args.domain == 'inner' else outer_dom
child = domain.child
parent = domain.parent
dim_xy = (domain.jtot, domain.itot)

exp_path = f'{work_path}/{args.domain}'


"""
Calculate exact vertical grid LES. The grid definition from
(LS)2D is slightly different from the one in MicroHH!
"""
gd = calc_vertical_grid_2nd(vgrid.z, vgrid.zsize)


"""
Read grid info, and calculate spatial interpolation factors.
"""
ds_cosmo = xr.open_dataset(f'{cosmo_path}/COSMO_CTRL_BC_nD_LES_XL.nc')
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


"""
Process hourly surface COSMO fields: SST (thl_bot) and qsat(SST) (qt_bot).
Unit conversions (T -> thl etc.) are done in `read_cosmo()`.
"""
print('Interpolating COSMO surface fields.')

fld_s = np.empty(dim_xy, float_type)

for t, time in enumerate(time_sec):

    hlp.interpolate_cosmo(fld_s, ds_cosmo.thl_sbot[t,:,:].values, if_s, None, float_type)
    fld_s.tofile(f'{exp_path}/thl_bot_in.{time:07d}')

    hlp.interpolate_cosmo(fld_s, ds_cosmo.qt_sbot[t,:,:].values, if_s, None, float_type)
    fld_s.tofile(f'{exp_path}/qt_bot_in.{time:07d}')


"""
Create `eurec4a_input.nc` file.
"""
settings = {
    'central_lon' : -57.7,
    'central_lat' : 13.3,
    'start_date'  : start_date,
    'end_date'    : end_date,
    'case_name'   : 'eurec4a_openbc',
    'era5_path'   : ls2d_era5_path}

era = ls2d.Read_era5(settings)
era.calculate_forcings(n_av=6, method='2nd')

ds_time = era.get_les_input(vgrid.z)
ds_time = ds_time.sel(lay=slice(0,135), lev=slice(0,136))
ds_mean = ds_time.mean(dim='time')

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


nc_file = nc4.Dataset(f'{exp_path}/eurec4a_input.nc', mode='w', datamodel='NETCDF4', clobber=True)
add_nc_dim('z', vgrid.kmax, nc_file)
add_nc_var('z', ('z'), nc_file, vgrid.z)

# Initial profiles. Not really used, but needed for `init` phase. The resulting
# initial 3D fields are overwritten by the interpolated fields from COSMO.
nc_init = nc_file.createGroup('init')
add_nc_var('thl', ('z'), nc_init, ds_mean.thl.values)
add_nc_var('qt', ('z'), nc_init, ds_mean.qt.values)
add_nc_var('u', ('z'), nc_init, ds_mean.u.values)
add_nc_var('v', ('z'), nc_init, ds_mean.v.values)

# Time dependent input. For now, spatially constant.
nc_tdep = nc_file.createGroup('timedep')
add_nc_dim('time_surface', time_sec.size, nc_tdep)
add_nc_dim('time_ls', time_sec.size, nc_tdep)

add_nc_var('time_surface', ('time_surface'), nc_tdep, time_sec)
add_nc_var('time_ls', ('time_ls'), nc_tdep, time_sec)

add_nc_var('p_sbot', ('time_surface'), nc_tdep, ds_time.ps.values)
add_nc_var('u_geo', ('time_ls', 'z'), nc_tdep, ds_time.ug.values)
add_nc_var('v_geo', ('time_ls', 'z'), nc_tdep, ds_time.vg.values)

# Radiation.
nc_rad = nc_file.createGroup('radiation')
add_nc_dim('lay', ds_mean.sizes['lay'], nc_rad)
add_nc_dim('lev', ds_mean.sizes['lev'], nc_rad)

# Radiation variables on LES grid.
xm_air = 28.97; xm_h2o = 18.01528; eps = xm_h2o / xm_air
h2o = ds_mean.qt / (eps - eps * ds_mean.qt)
add_nc_var('h2o', ('z'), nc_init, h2o)
add_nc_var('o3',  ('z'), nc_init, ds_mean.o3*1e-6)

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
add_nc_var('z_lay', ('lay'), nc_rad, ds_mean.z_lay)
add_nc_var('z_lev', ('lev'), nc_rad, ds_mean.z_lev)
add_nc_var('p_lay', ('lay'), nc_rad, ds_mean.p_lay)
add_nc_var('p_lev', ('lev'), nc_rad, ds_mean.p_lev)
add_nc_var('t_lay', ('lay'), nc_rad, ds_mean.t_lay)
add_nc_var('t_lev', ('lev'), nc_rad, ds_mean.t_lev)
add_nc_var('h2o',   ('lay'), nc_rad, ds_mean.h2o_lay)
add_nc_var('o3',    ('lay'), nc_rad, ds_mean.o3_lay*1e-6)

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
# `_ext` is appended to the name; this way the density created by `microhh init` can be overwritten with this one.
thermo.save_basestate_density(bs['rho'], bs['rhoh'], f'{exp_path}/rhoref_ext.0000000')


"""
Create 2D field with coriolis frequency.
"""
omega = 2*np.pi/86400
fc = 2 * omega * np.sin(np.deg2rad(domain.proj.lat))
fc.astype(float_type).tofile(f'{exp_path}/fc.0000000')


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
ini['grid']['ktot'] = vgrid.kmax

ini['grid']['xsize'] = domain.xsize
ini['grid']['ysize'] = domain.ysize
ini['grid']['zsize'] = vgrid.zsize

ini['buffer']['zstart'] = zstart_buffer
ini['buffer']['loadfreq'] = 3600 if args.domain == 'outer' else 600

ini['time']['endtime'] = time_sec[-1]
ini['time']['datetime_utc'] = start_date.strftime('%Y-%m-%d %H:%M:%S')

ini['boundary_lateral']['loadfreq'] = domain.lbc_freq
ini['boundary_lateral']['n_sponge'] = domain.n_sponge
ini['boundary_lateral']['tau_sponge'] = domain.n_sponge * domain.dx / (10 * 5)

if args.domain == 'outer':
    ini['boundary_lateral']['slist'] = ['thl', 'qt', 'qr']
    ini['boundary_lateral']['sw_recycle[east]'] = True
    ini['boundary_lateral']['sw_recycle[north]'] = True
else:
    ini['boundary_lateral']['slist'] = ['thl', 'qt', 'qr', 'nr']

if child is not None:
    ini['subdomain']['sw_subdomain'] = True

    ini['subdomain']['xstart'] = child.xstart_in_parent
    ini['subdomain']['ystart'] = child.ystart_in_parent
    ini['subdomain']['xend'] = child.xstart_in_parent + child.xsize
    ini['subdomain']['yend'] = child.ystart_in_parent + child.ysize

    ini['subdomain']['grid_ratio'] = int(domain.dx / child.dx)
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
ini['cross']['xy'] = [10, 400, 1000, 1500, 3000, 5000, vgrid.zsize]
ini['cross']['xz'] = domain.ysize/2
ini['cross']['yz'] = domain.xsize/2

# Column stats for BCO.
x,y = domain.proj.to_xy(lon=-59.432, lat=13.165)
if x < 0 or x > domain.xsize or x < 0 or y > domain.ysize:
    print('Warning: BCO coordinates outside domain. Setting column to domain center')
    x = domain.xsize / 2
    y = domain.ysize / 2

ini['column']['coordinates[x]'] = x
ini['column']['coordinates[y]'] = y

io.check_ini(ini)
io.save_ini(ini, f'{exp_path}/eurec4a.ini')


"""
Link/copy required RRTMGP lookup tables.
"""
rrtmgp_path = f'{microhh_path}/rte-rrtmgp-cpp/'
rrtmgp_data_path = f'{microhh_path}/rte-rrtmgp-cpp/rrtmgp-data'

to_copy = [
        (f'{gpt_veerman_path}/rrtmgp-gas-lw-g056-cf2.nc', 'coefficients_lw.nc'),
        (f'{gpt_veerman_path}/rrtmgp-gas-sw-g049-cf2.nc', 'coefficients_sw.nc'),
        (f'{rrtmgp_data_path}/rrtmgp-clouds-lw.nc', 'cloud_coefficients_lw.nc'),
        (f'{rrtmgp_data_path}/rrtmgp-clouds-sw.nc', 'cloud_coefficients_sw.nc'),
        (f'{rrtmgp_path}/data/aerosol_optics.nc', 'aerosol_optics.nc')]

for f in to_copy:
    target = f'{exp_path}/{f[1]}'
    if not os.path.exists(target):
        shutil.copy(f[0], target)


if args.domain == 'outer':
    """
    Create input from COSMO 3D files.
    """

    cosmo = xr.open_dataset(f'{cosmo_path}/COSMO_CTRL_BC_nD_LES_XL.nc')
    cosmo = cosmo.sel(time=slice(start_date, end_date))

    fields_cosmo = {
        'u': cosmo.u[:,:,:,:],
        'v': cosmo.v[:,:,:,:],
        'w': cosmo.w[:,:,:,:],
        'thl': cosmo.thl[:,:,:,:],
        'qt': cosmo.qt[:,:,:,:],
        'qr': cosmo.qr[:,:,:,:],
    }

    p_cosmo = cosmo.p[:,:,:,:]
    z_cosmo = np.broadcast_to(cosmo.z.values[None,:,None,None], cosmo.thl.shape)
    time_cosmo = cosmo.time - cosmo.time[0]
    time_cosmo = (time_cosmo.values.astype(np.float32)/1e9).astype(np.int32)

    # Standard dev. of Gaussian filter applied to interpolated fields (m).
    #sigma_h = 10_000
    sigma_h = 0.

    # MicroHHpy does not sprechen Deutsch. Just pretent it's ERA5...
    create_era5_input(
        fields_cosmo,
        cosmo.rlon.data,
        cosmo.rlat.data,
        z_cosmo,
        p_cosmo,
        time_cosmo,
        vgrid.z,
        vgrid.zsize,
        zstart_buffer,
        bs['rho'],
        bs['rhoh'],
        outer_dom,
        sigma_h,
        perturb_size=3,
        perturb_amplitude={'thl': 0.1, 'qt': 0.1e-3},
        perturb_max_height=1000,
        clip_at_zero=['qt', 'qr'],
        name_suffix='ext',
        output_dir=exp_path,
        ntasks=8,
        dtype=float_type)


else:
    """
    Create/copy input from parent domain.
    """
    parent_exp_path = f'{work_path}/outer'

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
            outer_dom.xsize,
            outer_dom.ysize,
            gd['z'],
            gd['zh'],
            outer_dom.itot,
            outer_dom.jtot,
            inner_dom.xsize,
            inner_dom.ysize,
            gd['z'],
            gd['zh'],
            inner_dom.itot,
            inner_dom.jtot,
            inner_dom.xstart_in_parent,
            inner_dom.ystart_in_parent,
            parent_exp_path,
            exp_path,
            float_type,
            method='nearest',
            name_suffix='ext')


    # Link boundary conditions from parent to child domain.
    # Only link the `_out` files, without `_out` they are LBCs used as input for the parent domain.
    def link_files(src_pattern):

        files = glob.glob(src_pattern)

        for f in files:
            src = Path(f).resolve()
            name = src.name.replace('_out', '')
            dst = Path(exp_path) / name
            dst.symlink_to(src)

    link_files(f'{parent_exp_path}/lbc_*_out.*')
    link_files(f'{parent_exp_path}/w_top_out.*')
    link_files(f'{parent_exp_path}/*_buffer_out.*')
