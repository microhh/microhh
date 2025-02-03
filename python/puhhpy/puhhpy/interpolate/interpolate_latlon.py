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

# Local library
from puhhpy.logger import logger
from .interpolate_kernels import Rect_to_curv_interpolation_factors

def interp_rect_to_curv_latlon(
        proj_pad,
        z_out,
        zh_out,
        fields_in,
        lon_in,
        lat_in,
        z_in,
        zh_in,
        rho_in,
        rhoh_in,
        output_dir='.',
        dtype=np.float64):
        """
        Interpolate the fields in the `fields_in` dictionary to the output lon/lat grid defined by the `proj` projection instance..

        Input requirements:
        - All 3D input fields must have dimensions `(level, lat, lon)`.
        - The input `lon_in` and `lat_in` must be one-dimensional, as only regulare lat/lon grids (e.g. ERA5) are supported.

        Output:
        - Each interpolated field is stored in binary format at `f'{output_dir}/{name}.0000000'`.

        Special handling for velocity fields:
        If `u`, `v`, and `w` are present in `fields_in`, they are interpolated and corrected 
        to be divergence-free. This correction uses the base state density and vertical grid definition
        from MicroHH, and ensures that the mean output vertical velocity matches the mean input vertical
        velocity, by applying small corrections to the domain-mean horizontal divergence.

        IMPORTANT:  
        For the divergence correction to work correctly, the input base state density and vertical grid definition  
        must precisely match the MicroHH definition. Use `puhhpy.spatial.Vertical_grid_2nd` to calculate the vertical grid  
        information correctly.

        Arguments:
        ----------
        proj_pad : `puhhpy.spatial.Projection` instance.
            Output projection definition, with ghost cells.
        z_out : np.ndarray, shape (1,)
            Output full levels (m).
        zh_out : np.ndarray, shape (1,)
            Output half levels (m).
        fields_in : dict
            Dictionary with field name : input array pairs.
            Input arrays should have shape (3,).
        lon_in : np.ndarray, shape (1,)
            Input longitude (degrees).
        lat_in : np.ndarray, shape (1,)
            Input latitude (degrees).
        z_in : np.ndarray, shape (3,)
            Input full level height (m).
        zh_in : np.ndarray, shape (3,)
            Input half level height (m).
        rhoref : np.ndarray, shape (1,)
            Input full level base state density (kg m-3).
        rhorefh : np.ndarray, shape (1,)
            Input half level base state density (kg m-3).
        output_dir : str
            Path where the binaries are saved.

        Returns:
        --------
        None
        """

        logger.info('interpolating t=0 from rectilinear to LES grid.')

        # Checks.
        if lon_in.ndim != 1 or lat_in.ndim != 1:
            logger.critical('input lat/lon has to be a 1D array!')
        if z_in.ndim != 3 or zh_in.ndim != 3:
            logger.critical('input height has to be a 3D array!')
        if z_out.ndim != 1 or zh_out.ndim != 1:
            logger.critical('output height has to be a 1D array!')

        # Strip off the mask (if present):
        if isinstance(lon_in, np.ma.masked_array):
             lon_in = lon_in.data
        if isinstance(lat_in, np.ma.masked_array):
             lat_in = lat_in.data

        # Interpolation indexes/factors.
        if 'u' in fields_in.keys():
            if_u = Rect_to_curv_interpolation_factors(
                lon_in, lat_in, proj_pad.lon_u, proj_pad.lat_u, dtype)

        if 'v' in fields_in.keys():
            if_v = Rect_to_curv_interpolation_factors(
                lon_in, lat_in, proj_pad.lon_v, proj_pad.lat_v, dtype)

        if_s = Rect_to_curv_interpolation_factors(
            lon_in, lat_in, proj_pad.lon, proj_pad.lat, dtype)

        """
        If all momentum fields are present, correct horizontal velocities.
        """
        processed_uvw = False
        if 'u' in fields_in.keys() and 'v' in fields_in.keys() and 'w' in fields_in.keys():
            logger.debug('all momentum fields present, correction horizontal divergence.')

            processed_uvw = True