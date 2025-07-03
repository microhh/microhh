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

import matplotlib.pyplot as pl
import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from scipy import interpolate

import ls2d

# Make sure `microhhpy` is in the Python path.
import microhhpy.thermo as thermo

# Custom scripts from current directory.
import helpers as hlp
from global_settings import cosmo_path, float_type
from domain_definition import vgrid

pl.close('all')

# NOTE: create background for full time period! This is quite an expensive operation,
# so don't repeat this for every test, change of time period, etc.
start = datetime(year=2020, month=2, day=1, hour=0)
end   = datetime(year=2020, month=2, day=12, hour=0)
ntime = int((end-start).total_seconds()/3600+1)

# Interpolate COSMO straight onto desired LES grid.
z_out = vgrid.z


"""
Read ERA5 for background profiles radiation.
"""
settings = {
    'central_lon' : -57.7,
    'central_lat' : 13.3,
    'start_date'  : start,
    'end_date'    : end,
    'case_name'   : 'eurec4a_openbc',
    'era5_path'   : '/home/scratch1/bart/LS2D_ERA5/'}

if 'era' not in locals():

    """
    Read ERA5 with (LS)2D.
    """
    era = ls2d.Read_era5(settings)
    era.calculate_forcings(n_av=6, method='2nd')
    
    les_input = era.get_les_input(z_out)
    les_input = les_input.sel(lay=slice(0,135), lev=slice(0,136))

    # Drop variables we don't need or don't overwrite with COSMO data.
    drop = [
        'zs', 'wls', 'dtthl_advec', 'dtqt_advec', 'dtu_advec', 'dtv_advec',
        't_soil', 'theta_soil' ,'type_soil', 'type_low_veg', 'type_high_veg', 'root_frac_low_veg',
        'root_frac_high_veg', 'lai_low_veg', 'lai_high_veg', 'c_low_veg', 'c_high_veg',
        'z0m', 'z0h', 'sst', 'ts', 'wth', 'wq']

    les_input = les_input.drop_vars(drop)


if 'cosmo' not in locals():

    """
    Read pre-processed COSMO dataset.
    """
    cosmo = xr.open_dataset(f'{cosmo_path}/COSMO_CTRL_BC_nD_LES.nc')

    # Location of inner 500 x 300 km2 domain.
    lon_slice = slice(-60, -55.4)
    lat_slice = slice(11.9, 14.6)

    cosmo = cosmo.sel(rlon=lon_slice, rlat=lat_slice)
    cosmo_m = cosmo.mean(dim=('rlat', 'rlon'))

    cosmo_m['T'] = cosmo_m['thl'] * thermo.exner(cosmo_m['p'].values)

    # Interpolate onto LES grid.
    cosmo_z = cosmo_m.interp(z=z_out, kwargs={'fill_value': 'extrapolate'})



"""
Blend COSMO (below ~20 km) and ERA5 (20 km -> TOA) profiles.
"""
for t, date in enumerate(cosmo.time):
    print(f'Parsing COSMO for {date.values}')

    def interp(x, xp, fp):
        """
        Interpolate `fp(xp)` to `fp(x)`. Use Scipy for interpolating
        pressure values (decreasing, so `np.interp` won't work) and
        (linear) extrapolating near surface.
        """
        return interpolate.interp1d(xp, fp, assume_sorted=False, fill_value='extrapolate')(x)

    # COSMO TOD is ~20km. Blend in ERA5 profiles between 50 < p < 100 hPa,
    # for RRTMGP background profiles temperature, h2o, and ozone.
    p_low = 10000
    p_high = 5000
    
    # Blend factor: 0 where COSMO is used, 1 where ERA is used, blend in between.
    fac_lay = np.maximum(0, np.minimum(1, (les_input.p_lay[t,:] - p_low) / (p_high - p_low)))
    fac_lev = np.maximum(0, np.minimum(1, (les_input.p_lev[t,:] - p_low) / (p_high - p_low)))

    # Interpolate COSMO to LES pressure grid.
    t_lay_cosmo = interp(les_input.p_lay[t,:],  cosmo_m.p[t,:], cosmo_m.T[t,:])
    t_lev_cosmo = interp(les_input.p_lev[t,:],  cosmo_m.p[t,:], cosmo_m.T[t,:])
    qt_lay_cosmo = interp(les_input.p_lay[t,:], cosmo_m.p[t,:], cosmo_m.qt[t,:])

    xm_air = 28.97; xm_h2o = 18.01528; eps = xm_h2o / xm_air
    h2o_lay_cosmo = qt_lay_cosmo / (eps - eps * qt_lay_cosmo)

    # Blend profiles.
    t_lay = fac_lay * les_input.t_lay[t,:] + (1-fac_lay) * t_lay_cosmo
    t_lev = fac_lev * les_input.t_lev[t,:] + (1-fac_lev) * t_lev_cosmo
    h2o_lay = fac_lay * les_input.h2o_lay[t,:] + (1-fac_lay) * h2o_lay_cosmo

    # Overwrite blended profiles. NOTE: no need to overwrite presure
    # or height values, as COSMO is interpolated to the ERA5 pressures.
    les_input['t_lay'][t,:] = t_lay
    les_input['t_lev'][t,:] = t_lev
    les_input['h2o_lay'][t,:] = h2o_lay

    # TODO: aerosols Mirjam.
    # .................

#"""
#Save in NetCDF, as input for other scripts.
#"""
## Cast float arrays to correct dtype.
#if dtype == np.float32:
#    for v in les_input.variables:
#        if les_input[v].dtype == np.float64:
#            les_input[v] = les_input[v].astype(dtype)
#    
#les_input.to_netcdf('eurec4a_mean_profiles.nc')