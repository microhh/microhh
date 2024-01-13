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
import grid_definition as gd
import helpers as hlp

pl.close('all')

def read_cosmo(date, cosmo_path):
    """
    Read COSMO for single time step (hour).
    """
    base_name = f'lffd{date.year:04d}{date.month:02d}{date.day:02d}{date.hour:02d}0000'
    ds_2d = xr.open_dataset(f'{cosmo_path}/COSMO_CTRL_BC_2D/{base_name}.nc').squeeze()
    ds_3d = xr.open_dataset(f'{cosmo_path}/COSMO_CTRL_BC_3D/{base_name}z.nc').squeeze()
    return ds_2d, ds_3d


if __name__ == '__main__':

    dtype = np.float32

    # Data paths
    cosmo_path = '/home/scratch2/bart/eurec4a_cosmo/'

    start = datetime(year=2020, month=2, day=1, hour=0)
    end   = datetime(year=2020, month=2, day=2, hour=0)

    # Short-cuts.
    hgrid = gd.hgrid_outer
    vgrid = gd.vgrid

    # Read first time step to get lat/lon COSMO grid.
    ds_2d, ds_3d = read_cosmo(start, cosmo_path)

    # Calculate interpolation factors at different locations staggered LES grid.
    # Easy, since COSMO is on regular lat/lon grid!
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


    u_les = np.empty((vgrid.kmax, hgrid.jtot, hgrid.itot), dtype)
    hlp.interpolate_cosmo(
            u_les, ds_3d.U.values,
            if_u.il, if_u.jl, if_z.kl,
            if_u.fx, if_u.fy, if_z.fz,
            dtype)

    vmin = -10
    vmax = 10

    pl.figure()
    pl.pcolormesh(ds_3d.rlon, ds_3d.rlat, ds_3d.U[0,:,:], vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r)
    pl.pcolormesh(hgrid.lon, hgrid.lat, u_les[0,:,:], vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r)
    pl.plot(hgrid.bbox_lon, hgrid.bbox_lat, 'w:', linewidth=1)
    pl.legend()
