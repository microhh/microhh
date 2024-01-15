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

# Custom scripts from current directory.
import helpers as hlp
from grid_definition import vgrid

pl.close('all')

cosmo_path = '/home/scratch2/bart/eurec4a_cosmo/'
dtype = np.float32

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
    'era5_path'   : '/home/scratch1/bart/LS2D/'}


if 'era' not in locals():
    era = ls2d.Read_era5(settings)
    era.calculate_forcings(n_av=8, method='2nd')
    
    les_input = era.get_les_input(z_out)
    les_input = les_input.sel(lay=slice(0,135), lev=slice(0,136))

    # Drop variables we don't need or don't overwrite with COSMO data.
    drop = [
        'zs', 'wls', 'dtthl_advec', 'dtqt_advec', 'dtu_advec', 'dtv_advec',
        't_soil', 'theta_soil' ,'type_soil', 'type_low_veg', 'type_high_veg', 'root_frac_low_veg',
        'root_frac_high_veg', 'lai_low_veg', 'lai_high_veg', 'c_low_veg', 'c_high_veg',
        'z0m', 'z0h', 'sst', 'ts', 'wth', 'wq']

    les_input = les_input.drop_vars(drop)


"""
Read COSMO data for domain averaged profiles, and blending of COSMO and ERA for radiation.
"""
# Location of inner 500 x 300 km2 domain.
lon_slice = slice(-60, -55.4)
lat_slice = slice(11.9, 14.6)

dates = pd.date_range(start, end, freq='1h')
for t,date in enumerate(dates):
    print(f'Parsing COSMO for {date}')

    ds_2d, ds_3d = hlp.read_cosmo(date, cosmo_path, lon_slice, lat_slice)

    # Average COSMO over 500x300 km2 domain.
    p = ds_3d.P.mean(axis=(1,2))
    u = ds_3d.U.mean(axis=(1,2))
    v = ds_3d.V.mean(axis=(1,2))
    T = ds_3d.T.mean(axis=(1,2))
    thl = ds_3d.thl.mean(axis=(1,2))
    qt = ds_3d.qt.mean(axis=(1,2))
    ps = ds_2d.PS.mean()
    z_cosmo = ds_3d.altitude.values

    def interp(x, xp, fp):
        """
        Interpolate `fp(xp)` to `fp(x)`. Use Scipy for interpolating
        pressure values (decreasing, so `np.interp` won't work) and
        (linear) extrapolating near surface.
        """
        return interpolate.interp1d(xp, fp, assume_sorted=False, fill_value='extrapolate')(x)

    # Interpolate to LES grid, and overwrite mean ERA5 profiles.
    les_input['thl'][t,:] = interp(z_out, z_cosmo, thl)
    les_input['qt'][t,:] = interp(z_out, z_cosmo, qt)
    les_input['u'][t,:] = interp(z_out, z_cosmo, u)
    les_input['v'][t,:] = interp(z_out, z_cosmo, v)
    les_input['p'][t,:] = interp(z_out, z_cosmo, p)
    les_input['ps'][t] = ps

    # COSMO TOD is ~20km. Blend in ERA5 profiles between 50 < p < 100 hPa,
    # for RRTMGP background profiles temperature, h2o, and ozone.
    p_low = 10000
    p_high = 5000
    
    # Blend factor: 0 where COSMO is used, 1 where ERA is used, blend in between.
    fac_lay = np.maximum(0, np.minimum(1, (les_input.p_lay[t,:] - p_low) / (p_high - p_low)))
    fac_lev = np.maximum(0, np.minimum(1, (les_input.p_lev[t,:] - p_low) / (p_high - p_low)))

    # Interpolate COSMO to LES pressure grid.
    t_lay_cosmo = interp(les_input.p_lay[t,:], p, T)
    t_lev_cosmo = interp(les_input.p_lev[t,:], p, T)
    qt_lay_cosmo = interp(les_input.p_lay[t,:], p, qt)

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

"""
Save in NetCDF, as input for other scripts.
"""
# Cast float arrays to correct dtype.
if dtype == np.float32:
    for v in les_input.variables:
        if les_input[v].dtype == np.float64:
            les_input[v] = les_input[v].astype(dtype)
    
les_input.to_netcdf('eurec4a_mean_profiles.nc')
