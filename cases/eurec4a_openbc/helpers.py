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

import xarray as xr
import numpy as np
from numba import jit, prange
from datetime import timedelta

import constants


@jit(nopython=True, nogil=True, fastmath=True)
def index_left(array, value, size):
    """
    Find index left of `value` in `array`. If out of bounds,
    an exception is raised.
    """

    for i in range(size-1):
        if array[i] <= value and array[i+1] > value:
            return i
    raise Exception("Interpolation out of bounds!")


@jit(nopython=True, nogil=True, fastmath=True)
def index_left_oob(array, value, size):
    """
    Find index left of `values` in `array`. If out of bounds,
    the first or second to last index is returned.
    """

    # Find index in domain.
    for i in range(size-1):
        if array[i] <= value and array[i+1] > value:
            return i

    # Out-of-bounds check.
    if array[0] > value:
        return 0
    elif array[-1] <= value:
        return size-2


@jit(nopython=True, nogil=True, fastmath=True)
def nearest_index(x_in, y_in, x_val, y_val, itot, jtot):
    """
    Find nearest (x_val, y_val) in `(x_in, y_in)` array.
    """
    i_min = 0
    j_min = 0
    d_min = 1e9

    for i in range(itot):
        for j in range(jtot):
            d = np.sqrt( (x_in[i]-x_val)**2 + (y_in[j]-y_val)**2 )
            if d < d_min:
                i_min = i
                j_min = j
                d_min = d

    return i_min, j_min


@jit(nopython=True, nogil=True, fastmath=True, parallel=True)
def calc_xy_interpolation_factors(
        il, jl,
        fx, fy,
        lon_cosmo, lat_cosmo,
        lon_les, lat_les,
        TF=np.float64):
    """
    Calculation of horizontal interpolation indexes and factors.
    """
    jtot_les, itot_les = lon_les.shape
    itot_cosmo = lon_cosmo.size
    jtot_cosmo = lat_cosmo.size

    for j in prange(jtot_les):
        for i in prange(itot_les):

            # Find index in ERA5 field west/south of LES grid point.
            il[j,i] = index_left(lon_cosmo, lon_les[j,i], itot_cosmo)
            jl[j,i] = index_left(lat_cosmo, lat_les[j,i], jtot_cosmo)

            # Interpolation factors.
            fx[j,i] = TF(1) - ((lon_les[j,i] - lon_cosmo[il[j,i]]) / (lon_cosmo[il[j,i]+1] - lon_cosmo[il[j,i]]))
            fy[j,i] = TF(1) - ((lat_les[j,i] - lat_cosmo[jl[j,i]]) / (lat_cosmo[jl[j,i]+1] - lat_cosmo[jl[j,i]]))


class Calc_xy_interpolation_factors:
    """
    Help class to calculate the horizontal interpolation settings
    for each (scalar, u, v,) LES grid location.
    """
    def __init__(
            self,
            lon_cosmo, lat_cosmo,
            lon_les, lat_les,
            itot, jtot,
            dtype):

        self.il = np.zeros((jtot, itot), dtype=np.uint16)
        self.jl = np.zeros((jtot, itot), dtype=np.uint16)

        # Interpolation factors:
        self.fx = np.zeros((jtot, itot), dtype=dtype)
        self.fy = np.zeros((jtot, itot), dtype=dtype)

        calc_xy_interpolation_factors(
                self.il, self.jl,
                self.fx, self.fy,
                lon_cosmo, lat_cosmo,
                lon_les, lat_les)


class Calc_z_interpolation_factors:
    def __init__(self, z_cosmo, z_les, dtype):
        """
        Help class to calculate the vertical interpolation factors.
        Easy, since COSMO is on a fixed height grid, so `z_cosmo` is 1D.
        """

        self.kl = np.zeros_like(z_les, np.uint16)
        self.fz = np.zeros_like(z_les, dtype)

        for k in range(z_les.size):
            self.kl[k] = index_left_oob(z_cosmo, z_les[k], z_cosmo.size)
            self.fz[k] = 1.- ((z_les[k] - z_cosmo[self.kl[k]]) / (z_cosmo[self.kl[k]+1] - z_cosmo[self.kl[k]]))


@jit(nopython=True, nogil=True, fastmath=True, parallel=True)
def interpolate_cosmo_3d(
        fld_out, fld_cosmo,
        il, jl, kl,
        fx, fy, fz,
        TF=np.float64):
    """
    Fast tri-linear interpolation of ERA5 onto LES grid.
    """
    ktot_les, jtot_les, itot_les = fld_out.shape
    ktot_cosmo, jtot_cosmo, itot_cosmo = fld_cosmo.shape

    for j in prange(jtot_les):
        for i in prange(itot_les):
            for k in prange(ktot_les):

                # Short cuts
                ill = il[j,i]
                jll = jl[j,i]
                kll = kl[k]

                fxl = fx[j,i]
                fxr = TF(1) - fxl

                fyl = fy[j,i]
                fyr = TF(1) - fyl

                fzl = fz[k]
                fzr = TF(1) - fzl

                fld_out[k,j,i] =  \
                        fxl * fyl * fzl * fld_cosmo[kll,   jll,   ill  ] + \
                        fxr * fyl * fzl * fld_cosmo[kll,   jll,   ill+1] + \
                        fxl * fyr * fzl * fld_cosmo[kll,   jll+1, ill  ] + \
                        fxr * fyr * fzl * fld_cosmo[kll,   jll+1, ill+1] + \
                        fxl * fyl * fzr * fld_cosmo[kll+1, jll,   ill  ] + \
                        fxr * fyl * fzr * fld_cosmo[kll+1, jll,   ill+1] + \
                        fxl * fyr * fzr * fld_cosmo[kll+1, jll+1, ill  ] + \
                        fxr * fyr * fzr * fld_cosmo[kll+1, jll+1, ill+1]


@jit(nopython=True, nogil=True, fastmath=True, parallel=True)
def interpolate_cosmo_2d(
        fld_out, fld_cosmo,
        il, jl,
        fx, fy,
        TF=np.float64):
    """
    Fast bi-linear interpolation of ERA5 onto LES grid.
    """
    jtot_les, itot_les = fld_out.shape
    jtot_cosmo, itot_cosmo = fld_cosmo.shape

    for j in prange(jtot_les):
        for i in prange(itot_les):

            # Short cuts
            ill = il[j,i]
            jll = jl[j,i]

            fxl = fx[j,i]
            fxr = TF(1) - fxl

            fyl = fy[j,i]
            fyr = TF(1) - fyl

            fld_out[j,i] =  \
                    fxl * fyl * fld_cosmo[jll,   ill  ] + \
                    fxr * fyl * fld_cosmo[jll,   ill+1] + \
                    fxl * fyr * fld_cosmo[jll+1, ill  ] + \
                    fxr * fyr * fld_cosmo[jll+1, ill+1]


def interpolate_cosmo(fld_les, fld_cosmo, if_xy, if_z=None, dtype=np.float64):
    """
    Interpolate 2D or 3D field from COSMO to LES grid.
    """
    if fld_les.ndim == 2:
        interpolate_cosmo_2d(
            fld_les, fld_cosmo,
            if_xy.il, if_xy.jl,
            if_xy.fx, if_xy.fy,
            dtype)

    elif fld_les.ndim == 3:
        interpolate_cosmo_3d(
            fld_les, fld_cosmo,
            if_xy.il, if_xy.jl, if_z.kl,
            if_xy.fx, if_xy.fy, if_z.fz,
            dtype)

    else:
        raise Exception('Can only interpolate 2D (xy) or 3D (xyz) fields!')


def check_grid_decomposition(itot, jtot, ktot, npx, npy):
    """
    Check whether grid / MPI decomposition is valid
    """

    err = False
    if itot%npx != 0:
        print('ERROR in grid: itot%npx != 0')
        err = True

    if itot%npy != 0:
        print('ERROR in grid: itot%npy != 0')
        err = True

    if jtot%npx != 0 and npy > 1:
        print('ERROR in grid: jtot%npx != 0')
        err = True

    if jtot%npy != 0:
        print('ERROR in grid: jtot%npy != 0')
        err = True

    if ktot%npx != 0:
        print('ERROt in grid: ktot%npx != 0')
        err = True

    if err:
        print('Grid: itot={}, jtot={}, ktot={}, npx={}, npy={}'.format(
            itot, jtot, ktot, npx, npy))
        raise Exception('Invalid grid configuration!')
    else:
        print('Grid: itot={}, jtot={}, ktot={}, npx={}, npy={}: OKAY!'.format(
            itot, jtot, ktot, npx, npy))


class Grid_stretched_manual:
    def __init__(self, ktot, dz0, heights, factors):
        self.ktot = ktot
        self.dz0  = dz0

        self.z = np.zeros(ktot)
        self.zh = np.zeros(ktot+1)
        self.dz = np.zeros(ktot)
        self.zsize = None

        self.z[0]  = dz0/2.
        self.dz[0] = dz0

        def index(z, goal):
            return np.where(z-goal>0)[0][0]-1

        for k in range(1, ktot):
            self.dz[k] = self.dz[k-1] * factors[index(heights, self.z[k-1])]
            self.z[k] = self.z[k-1] + self.dz[k]

        self.zsize = self.z[ktot-1] + 0.5*self.dz[ktot-1]

        self.zh[1:-1] = self.z[1:] - self.z[:-1]
        self.zh[-1] = self.zsize


def esat_liq(T):
    """
    Calculate saturation vapor pressure at absolute temperature T
    """
    x = np.maximum(-75., T-constants.T0);
    return constants.c00+x*(constants.c10+x*(constants.c20+x*(constants.c30+x*(constants.c40+x*(constants.c50+x*(constants.c60+x*(constants.c70+x*(constants.c80+x*(constants.c90+x*constants.c100)))))))))


def qsat_liq(T, p):
    """
    Calculate saturation specific humidity from absolute temperature and pressure.
    """
    return constants.ep*esat_liq(T)/(p-(1.-constants.ep)*esat_liq(T))


def exner(p):
    """
    Calculate exner function as function of absolute pressure p.
    """
    return (p/constants.p0)**(constants.Rd/constants.cp)


def read_cosmo(date, cosmo_path, lon_slice, lat_slice):
    """
    Read COSMO for single time step (hour).
    """
    base_name = f'lffd{date.year:04d}{date.month:02d}{date.day:02d}{date.hour:02d}0000'

    ds_2d = xr.open_dataset(f'{cosmo_path}/COSMO_CTRL_BC_2D/{base_name}.nc').squeeze()
    ds_3d = xr.open_dataset(f'{cosmo_path}/COSMO_CTRL_BC_3D/{base_name}z.nc').squeeze()

    ds_2d = ds_2d.sel(rlon=lon_slice, rlat=lat_slice)
    d3_2d = ds_3d.sel(rlon=lon_slice, rlat=lat_slice)

    # Calculate derived properties.
    dims_2d = ('rlat', 'rlon')
    dims_3d = ('altitude', 'rlat', 'rlon')

    exner_2d = exner(ds_2d['PS'].values)
    exner_3d = exner(ds_3d['P'].values)

    ds_2d['qsat_s'] = (dims_2d, qsat_liq(ds_2d['T_S'].values, ds_2d['PS'].values))
    ds_2d['thl_s'] = (dims_2d, ds_2d['T_S'].values / exner_2d)

    ds_3d['qt'] = (dims_3d, ds_3d['QV'].values + ds_3d['QC'].values + ds_3d['QI'].values)
    ds_3d['qr'] = (dims_3d, ds_3d['QR'].values + ds_3d['QS'].values)
    ds_3d['th'] = (dims_3d, ds_3d['T'].values / exner_3d)

    thl = ds_3d['th'].values - \
            constants.Lv * ds_3d['QC'].values / (constants.cp * exner_3d) - \
            constants.Ls * ds_3d['QI'].values / (constants.cp * exner_3d)

    ds_3d['thl'] = (dims_3d, thl)

    return ds_2d, ds_3d
