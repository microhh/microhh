#
#  MicroHH
#  Copyright (c) 2011-2024 Chiel van Heerwaarden
#  Copyright (c) 2011-2024 Thijs Heus
#  Copyright (c) 2014-2024 Bart van Stratum
#
#  This file is part of MicroHH
#
#  MicroHH is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MicroHH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
#

# Standard library

# Third-party.
import numpy as np
from numba import jit, prange

# Local library
from puhhpy.logger import logger


@jit(nopython=True, nogil=True, fastmath=True)
def _index_left(array, value, size):
    """
    Find index left of `value` in `array`. If out of bounds, an exception is raised.
    """

    for i in range(size-1):
        if array[i] <= value and array[i+1] > value:
            return i

    # Don't use `logging` in Numba kernels.
    raise Exception("interpolation out of bounds!")


@jit(nopython=True, nogil=True, fastmath=True)
def _index_left_oob(array, value, size):
    """
    Find index left of `value` in `array`.
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
def _nearest_index(x_in, y_in, x_val, y_val, itot, jtot):
    """
    Find nearest `i,j` index given two 1D arrays and `x,y` values.
    Written out using Numba, to prevent having to do a `np.meshgrid` on potentially large 1D arrays.
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
def _calc_rect_to_curv_interpolation_factors(
        il, jl,
        im, jm,
        fx, fy,
        lon_era, lat_era,
        lon_les, lat_les,
        TF=np.float64):
    """
    Numba kernel to calculate the horizontal interpolation indexes and factors.
    """
    jtot_les, itot_les = lon_les.shape
    itot_era = lon_era.size
    jtot_era = lat_era.size

    for j in prange(jtot_les):
        for i in prange(itot_les):

            # Find index in ERA5 field west/south of LES grid point.
            il[j,i] = _index_left(lon_era, lon_les[j,i], itot_era)
            jl[j,i] = _index_left(lat_era, lat_les[j,i], jtot_era)

            # Interpolation factors.
            fx[j,i] = TF(1) - ((lon_les[j,i] - lon_era[il[j,i]]) / (lon_era[il[j,i]+1] - lon_era[il[j,i]]))
            fy[j,i] = TF(1) - ((lat_les[j,i] - lat_era[jl[j,i]]) / (lat_era[jl[j,i]+1] - lat_era[jl[j,i]]))

            # Find nearest grid point in ERA.
            im[j,i], jm[j,i] = _nearest_index(lon_era, lat_era, lon_les[j,i], lat_les[j,i], itot_era, jtot_era)


class Rect_to_curv_interpolation_factors:
    """
    Class to calculate horizontal interpolation indexes and factors, for interpolating from a
    rectilinear input grid (1D arrays for lat/lon) onto a non-rectilinear grid (2D arrays for lat/lon).

    Arguments:
    ---------
    lon_in : np.ndarray, shape (1,)
        Input longitude on rectilinear grid.
    lat_in : np.ndarray, shape (1,)
        Input latitude on rectilinear grid.
    lon_out : np.ndarray, shape (2,)
        Output longitude on non-rectilinear grid.
    lat_out : np.ndarray, shape (2,)
        Output latitude on non-rectilinear grid.
    dtype : Numpy float type
        Numpy floating point datatype.

    Returns:
    -------
    None
    """
    def __init__(
            self,
            lon_in, lat_in,
            lon_out, lat_out,
            dtype):

        jtot, itot = lon_out.shape

        # Index left.
        self.il = np.zeros((jtot, itot), dtype=np.uint32)
        self.jl = np.zeros((jtot, itot), dtype=np.uint32)

        # Nearest grid point.
        self.im = np.zeros((jtot, itot), dtype=np.uint32)
        self.jm = np.zeros((jtot, itot), dtype=np.uint32)

        # Interpolation factors.
        self.fx = np.zeros((jtot, itot), dtype=dtype)
        self.fy = np.zeros((jtot, itot), dtype=dtype)

        _calc_rect_to_curv_interpolation_factors(
                self.il, self.jl,
                self.im, self.jm,
                self.fx, self.fy,
                lon_in, lat_in,
                lon_out, lat_out,
                TF=dtype)