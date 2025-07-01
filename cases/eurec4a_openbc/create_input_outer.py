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

import xarray as xr
import pandas as pd
import numpy as np

# Make sure `microhhpy` is in the Python path.
import microhhpy.thermo as thermo
from microhhpy.openbc import create_era5_input

# Custom modules, from `microhh/python/`.
#import microhh_tools as mht
#import microhh_lbc_tools as mlt

# Custom scripts from same directory.
import helpers as hlp
import constants
from global_settings import float_type, work_path, cosmo_path, start_date, end_date

# Grid (horizontal/vertical) definition.
from domain_definition import vgrid, outer_dom, zstart_buffer

"""
Read basestate
"""
rho, rhoh = thermo.read_basestate_density(f'{work_path}/rhoref_0.0000000', dtype=float_type)

"""
Parse COSMO data.
"""
cosmo = xr.open_dataset(f'{cosmo_path}/COSMO_CTRL_BC_nD_LES.nc')
cosmo = cosmo.sel(time=slice(start_date, end_date))

fields_era = {
    'u': cosmo.u[:,:,:,:],
    'v': cosmo.v[:,:,:,:],
    'w': cosmo.w[:,:,:,:],
    'thl': cosmo.thl[:,:,:,:],
    'qt': cosmo.qt[:,:,:,:],
    'qr': cosmo.qr[:,:,:,:],
}

p_era = cosmo.p[:,:,:,:]
z_era = np.broadcast_to(cosmo.z.values[None,:,None,None], cosmo.thl.shape)
time_era = cosmo.time - cosmo.time[0]
time_era = (time_era.values.astype(np.float32)/1e9).astype(np.int32)

# Standard dev. of Gaussian filter applied to interpolated fields (m).
#sigma_h = 10_000
sigma_h = 0.

# MicroHHpy does not sprechen Deutsch. Just pretent it's ERA5...
create_era5_input(
    fields_era,
    cosmo.rlon.data,
    cosmo.rlat.data,
    z_era,
    p_era,
    time_era,
    vgrid.z,
    vgrid.zsize,
    zstart_buffer,
    rho,
    rhoh,
    outer_dom,
    sigma_h,
    perturb_size=3,
    perturb_amplitude={'thl': 0.1, 'qt': 0.1e-3},
    perturb_max_height=1000,
    name_suffix='0',
    output_dir=work_path,
    ntasks=12,
    dtype=float_type)
