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

from datetime import datetime, timedelta

import matplotlib.pyplot as pl
import xarray as xr
import numpy as np

# Custom modules, from `microhh/python/`.
import microhh_tools as mht
import microhh_lbc_tools as mlt

# Custom scripts from same directory.
import helpers as hlp

pl.close('all')

dtype = np.float32

# Data paths
cosmo_path = '/home/scratch2/bart/eurec4a_cosmo/'
work_path = 'test_case'

#import grid_definition as gd  # Full domain.
import grid_definition_dev as gd  # Develop domain.

start = datetime(year=2020, month=2, day=1, hour=0)
end   = datetime(year=2020, month=2, day=2, hour=0)
nt    = int((end - start).total_seconds()/3600)+1

# Short-cuts.
hgrid = gd.hgrid_outer_pad
vgrid = gd.vgrid

dim_xy = (hgrid.jtot, hgrid.itot)
dim_xyz = (vgrid.ktot, hgrid.jtot, hgrid.itot)

ngc = gd.n_ghost
nbuf = gd.n_buffer
nlbc = ngc + nbuf

lon_slice = slice(hgrid.lon.min()-0.5, hgrid.lon.max()+0.5)
lat_slice = slice(hgrid.lat.min()-0.5, hgrid.lat.max()+0.5)


"""
Read grid info, and calculate spatial interpolation factors.
"""
ds_2d, ds_3d = hlp.read_cosmo(start, cosmo_path, lon_slice, lat_slice)

# Calculate interpolation factors at different locations staggered LES grid.
if_u = hlp.Calc_xy_interpolation_factors(
        ds_2d.rlon.values, ds_2d.rlat.values,
        hgrid.lon_u, hgrid.lat_u,
        hgrid.itot, hgrid.jtot,
        dtype)

if_v = hlp.Calc_xy_interpolation_factors(
        ds_2d.rlon.values, ds_2d.rlat.values,
        hgrid.lon_v, hgrid.lat_v,
        hgrid.itot, hgrid.jtot,
        dtype)

if_s = hlp.Calc_xy_interpolation_factors(
        ds_2d.rlon.values, ds_2d.rlat.values,
        hgrid.lon, hgrid.lat,
        hgrid.itot, hgrid.jtot,
        dtype)

# Calculate vertical interpolation factors.
if_z = hlp.Calc_z_interpolation_factors(ds_3d.altitude.values, vgrid.z, dtype)


"""
Interpolated fields contain ghost cells. Define Numpy slices
to obtain the inner domain, LBCs, et cetera.
"""
s_inner = np.s_[:, +ngc:-ngc, +ngc:-ngc]
s_inner_2d = np.s_[+ngc:-ngc, +ngc:-ngc]

ss_west = np.s_[:, :, :nlbc]
ss_east = np.s_[:, :, -nlbc:]
ss_south = np.s_[:, :nlbc, :]
ss_north = np.s_[:, -nlbc:, :]

su_west = np.s_[:, :, :nlbc+1]
su_east = np.s_[:, :, -nlbc:]
su_south = np.s_[:, :nlbc, :]
su_north = np.s_[:, -nlbc:, :]

sv_west = np.s_[:, :, :nlbc]
sv_east = np.s_[:, :, -nlbc:]
sv_south = np.s_[:, :nlbc+1, :]
sv_north = np.s_[:, -nlbc:, :]


"""
Process hourly COSMO data.
"""
date = start
while date <= end:
    print(f'Processing {date}')
    time = int((date - start).total_seconds())

    ds_2d, ds_3d = hlp.read_cosmo(date, cosmo_path, lon_slice, lat_slice)

    """
    Process surface boundary conditions.
    """
    thl_s = np.empty(dim_xy, dtype)
    hlp.interpolate_cosmo(thl_s, ds_2d.thl_s.values, if_s, None, dtype)
    thl_s[s_inner_2d].tofile(f'{work_path}/thl_sbot_in.{time:07d}')

    qt_s = np.empty(dim_xy, dtype)
    hlp.interpolate_cosmo(qt_s, ds_2d.qsat_s.values, if_s, None, dtype)
    qt_s[s_inner_2d].tofile(f'{work_path}/qt_sbot_in.{time:07d}')


    """
    Process atmospheric fields and LBCs.
    """
    tmp = np.empty(dim_xyz, dtype)
    for fld in ['thl', 'qt', 'qr']:
        hlp.interpolate_cosmo(tmp, ds_3d[fld].values, if_s, if_z, dtype)
        if date == start:
            tmp[s_inner].tofile(f'{work_path}/{fld}_0.{time:07d}')
    del tmp

    u = np.empty(dim_xyz, dtype)
    v = np.empty(dim_xyz, dtype)

    hlp.interpolate_cosmo(u, ds_3d['u'].values, if_u, if_z, dtype)
    hlp.interpolate_cosmo(v, ds_3d['v'].values, if_v, if_z, dtype)

    if date == start:
        u[s_inner].tofile(f'{work_path}/u_0.{time:07d}')
        v[s_inner].tofile(f'{work_path}/v_0.{time:07d}')

    # Calculate vertical velocity that results in `dui/dxi=0`.
    # TODO..................









    date += timedelta(hours=1)
