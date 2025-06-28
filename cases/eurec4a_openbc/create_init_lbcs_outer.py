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

# Custom modules, from `microhh/python/`.
#import microhh_tools as mht
#import microhh_lbc_tools as mlt

# Custom scripts from same directory.
import helpers as hlp
import constants
from global_settings import float_type, work_path, cosmo_path, start_date, end_date

# Grid (horizontal/vertical) definition.
from domain_definition import outer_dom as domain

# More help scripts, from (LS)2D repo.
#from microhh_grid import Vertical_grid_2nd
#from microhh_thermo import Basestate_moist

"""
Settings.
"""
perturb_size = 3    # Perturb in nxnxn blocks
perturb_ampl = dict(thl = 0.1, qt =  0.1e-3)


"""
Read basestate
"""
bs = thermo.read_moist_basestate(f'{work_path}/thermo_basestate_0.0000000')


"""
Parse COSMO data.

"""








#fields_era = {
#    'u': era5.u[:,:,:,:],
#    'v': era5.v[:,:,:,:],
#    'w': era5.wls[:,:,:,:],
#    'thl': era5.thl[:,:,:,:],
#    'qt': era5.qt[:,:,:,:],
#}
#
#p_era = era5.p[:,:,:,:]
#z_era = era5.z[:,:,:,:]
#time_era = era5.time_sec
#
## Standard dev. of Gaussian filter applied to interpolated fields (m).
#sigma_h = 10_000
#

# MicroHH does not sprechen Deutsch. Just pretent it's ERA5...
#create_era5_input(
#    fields_era,
#    era5.lons.data,   # Strip off array masks.
#    era5.lats.data,
#    z_era,
#    p_era,
#    time_era,
#    gd['z'],
#    gd['zsize'],
#    bs['rho'],
#    bs['rhoh'],
#    dom0,
#    sigma_h,
#    perturb_size=4,
#    perturb_amplitude={'thl': 0.1, 'qt': 0.1e-3},
#    name_suffix='era5',
#    output_dir=work_dir,
#    ntasks=16,
#    dtype=TF)


















#"""
#Short-cuts.
#"""
#hgrid = gd.hgrid_outer_pad
#vgrid = gd.vgrid
#ktot = vgrid.ktot
#
#dim_xy = (hgrid.jtot, hgrid.itot)
#dim_xyz = (vgrid.ktot, hgrid.jtot, hgrid.itot)
#
#ngc = gd.n_ghost
#nbuf = gd.n_buffer
#nlbc = ngc + nbuf
#
#lon_slice = slice(hgrid.lon.min()-0.5, hgrid.lon.max()+0.5)
#lat_slice = slice(hgrid.lat.min()-0.5, hgrid.lat.max()+0.5)
#
## Vertical grid MicroHH.
#uhh_gd = Vertical_grid_2nd(vgrid.z, vgrid.zsize)
#
## Time (not really needed...)
#dates = pd.date_range(start_date, end_date, freq='h')
#time_sec = np.array((dates-dates[0]).total_seconds()).astype(np.int32)
#
#
#"""
#Read base state density created by `create_basestate.py`.
#"""
#bs = np.fromfile(f'{work_path}/rhoref_0.0000000', dtype=dtype)
#rho = bs[:ktot]
#rhoh = bs[ktot:]
#
#
#"""
#Read grid info, and calculate spatial interpolation factors.
#"""
#ds_2d, ds_3d = hlp.read_cosmo(start_date, cosmo_path, lon_slice, lat_slice)
#
## Calculate interpolation factors at different locations staggered LES grid.
#if_u = hlp.Calc_xy_interpolation_factors(
#        ds_2d.rlon.values, ds_2d.rlat.values,
#        hgrid.lon_u, hgrid.lat_u,
#        hgrid.itot, hgrid.jtot,
#        dtype)
#
#if_v = hlp.Calc_xy_interpolation_factors(
#        ds_2d.rlon.values, ds_2d.rlat.values,
#        hgrid.lon_v, hgrid.lat_v,
#        hgrid.itot, hgrid.jtot,
#        dtype)
#
#if_s = hlp.Calc_xy_interpolation_factors(
#        ds_2d.rlon.values, ds_2d.rlat.values,
#        hgrid.lon, hgrid.lat,
#        hgrid.itot, hgrid.jtot,
#        dtype)
#
## Calculate vertical interpolation factors.
#if_z = hlp.Calc_z_interpolation_factors(ds_3d.altitude.values, vgrid.z, dtype)
#
#
#"""
#Interpolated fields contain ghost cells. Define Numpy slices
#to obtain the inner domain, LBCs (ghost + sponge cells), et cetera.
#"""
#s_inner = np.s_[:, +ngc:-ngc, +ngc:-ngc]
#s_inner_2d = np.s_[+ngc:-ngc, +ngc:-ngc]
#
#ss_west = np.s_[:, :, :nlbc]
#ss_east = np.s_[:, :, -nlbc:]
#ss_south = np.s_[:, :nlbc, :]
#ss_north = np.s_[:, -nlbc:, :]
#
#su_west = np.s_[:, :, :nlbc+1]
#su_east = np.s_[:, :, -nlbc:]
#su_south = np.s_[:, :nlbc, :]
#su_north = np.s_[:, -nlbc:, :]
#
#sv_west = np.s_[:, :, :nlbc]
#sv_east = np.s_[:, :, -nlbc:]
#sv_south = np.s_[:, :nlbc+1, :]
#sv_north = np.s_[:, -nlbc:, :]
#
#
#"""
#Create Xarray dataset with correct dimensions/coordinates/.. for LBCs
#"""
#lbc_ds = mlt.get_lbc_xr_dataset(
#        ('u', 'v', 'w', 'thl', 'qt', 'qr', 'nr'),
#        gd.hgrid_outer.xsize,
#        gd.hgrid_outer.ysize,
#        gd.hgrid_outer.itot,
#        gd.hgrid_outer.jtot,
#        vgrid.z, vgrid.zh,
#        time_sec,
#        gd.n_ghost,
#        gd.n_buffer,
#        dtype)
#
#
#"""
#Process hourly COSMO data.
#"""
#for t, date in enumerate(dates):
#    timer = hlp.Timer()
#
#    print(f'Processing {date}, ', end='')
#
#    ds_2d, ds_3d = hlp.read_cosmo(date, cosmo_path, lon_slice, lat_slice)
#
#    """
#    Process atmospheric fields and LBCs.
#    """
#    tmp = np.empty(dim_xyz, dtype)
#
#    for fld in ['thl', 'qt', 'qr']:
#        hlp.interpolate_cosmo(tmp, ds_3d[fld].values, if_s, if_z, dtype)
#
#        # Add random perturbations to entire 3D field, and make sure that fld >= 0.
#        if fld in perturb_ampl.keys():
#            hlp.block_perturb_fld(tmp, perturb_size, perturb_ampl[fld])
#            tmp[tmp<0] = 0.
#
#        if date == start_date:
#            tmp[s_inner].tofile(f'{work_path}/{fld}_0.0000000')
#
#        # Store LBCs.
#        lbc_ds[f'{fld}_west'][t,:,:,:] = tmp[ss_west]
#        lbc_ds[f'{fld}_east'][t,:,:,:] = tmp[ss_east]
#        lbc_ds[f'{fld}_south'][t,:,:,:] = tmp[ss_south]
#        lbc_ds[f'{fld}_north'][t,:,:,:] = tmp[ss_north]
#
#    del tmp
#
#    u = np.empty(dim_xyz, dtype)
#    v = np.empty(dim_xyz, dtype)
#
#    hlp.interpolate_cosmo(u, ds_3d['U'].values, if_u, if_z, dtype)
#    hlp.interpolate_cosmo(v, ds_3d['V'].values, if_v, if_z, dtype)
#
#    if date == start_date:
#        u[s_inner].tofile(f'{work_path}/u_0.0000000')
#        v[s_inner].tofile(f'{work_path}/v_0.0000000')
#
#    # Store LBCs.
#    lbc_ds.u_west[t,:,:,:] = u[su_west]
#    lbc_ds.u_east[t,:,:,:] = u[su_east]
#    lbc_ds.u_south[t,:,:,:] = u[su_south]
#    lbc_ds.u_north[t,:,:,:] = u[su_north]
#
#    lbc_ds.v_west[t,:,:,:] = v[sv_west]
#    lbc_ds.v_east[t,:,:,:] = v[sv_east]
#    lbc_ds.v_south[t,:,:,:] = v[sv_south]
#    lbc_ds.v_north[t,:,:,:] = v[sv_north]
#
#    # Calculate vertical velocity field that satisfies `dui/dxi==0`.
#    w = np.zeros((vgrid.ktot+1, hgrid.jtot, hgrid.itot), dtype)
#
#    hlp.calc_w_from_uv(
#            w, u, v,
#            rho, rhoh,
#            uhh_gd.dz[uhh_gd.kstart:uhh_gd.kend],
#            1./hgrid.dx,
#            1./hgrid.dy,
#            ngc, hgrid.itot-ngc,
#            ngc, hgrid.jtot-ngc,
#            vgrid.ktot)
#
#    w_top = w[-1, ngc:-ngc, ngc:-ngc]
#    w = w[:-1, +ngc:-ngc, +ngc:-ngc]
#
#    if date == start_date:
#        w.tofile(f'{work_path}/w_0.0000000'.format(work_path))
#
#    # Save w_top as boundary conditions for MicroHH.
#    w_top.tofile(f'{work_path}/w_top.{time_sec[t]:07d}')
#    print(f'<w_top> = {(w_top.mean()*100):+.9f} cm/s, ', end='')
#
#    del u,v,w
#
#    # Interpolate `w` from COSMO to get the LBC values.
#    # NOTE: this might result in a small inconsistency between
#    #       `w` used for the initial field and the LBCs,
#    #       since the initial fields calculated to guarantee `dui/dxi==0`.
#    w = np.empty(dim_xyz, dtype)
#    hlp.interpolate_cosmo(w, ds_3d['W'].values, if_s, if_z, dtype)
#
#    # Store LBCs.
#    lbc_ds.w_west[t,:-1,:,:] = w[ss_west]
#    lbc_ds.w_east[t,:-1,:,:] = w[ss_east]
#    lbc_ds.w_south[t,:-1,:,:] = w[ss_south]
#    lbc_ds.w_north[t,:-1,:,:] = w[ss_north]
#
#    # Calculate geostrophic wind components.
#    #fc = 2*(2*np.pi/86400.) * np.sin(np.deg2rad(hgrid.lat))
#
#    #if date == start:
#    #    fc[s_inner_2d].tofile(f'{work_path}/fc.0000000')
#
#    #p = np.empty(dim_xyz, dtype)
#    #hlp.interpolate_cosmo(p, ds_3d['P'].values, if_s, if_z, dtype)
#
#    #vg =  1/(uhh_bs.rho[1:-1,None,None] * fc[None,:,:]) * np.gradient(p, axis=2) / hgrid.dx
#    #ug = -1/(uhh_bs.rho[1:-1,None,None] * fc[None,:,:]) * np.gradient(p, axis=1) / hgrid.dy
#
#    #pl.figure()
#    #pl.plot(ug.mean(axis=(1,2)), vgrid.z, 'b-o', label='ug, COSMO')
#    #pl.plot(vg.mean(axis=(1,2)), vgrid.z, 'r-o', label='vg, COSMO')
#    #pl.plot(ds_bg.ug[0,:], ds_bg.z, 'b--', label='ug, ERA')
#    #pl.plot(ds_bg.vg[0,:], ds_bg.z, 'r--', label='vg, ERA')
#    #pl.legend()
#
#    elapsed = timer.stop()
#    print(f'processing took {elapsed}')
#
#
#"""
#Save LBCs in binary format.
#"""
#mlt.write_lbcs_as_binaries(
#        lbc_ds, dtype, work_path)
#
## Bonus (aka debug): save in NetCDF.
#lbc_ds.to_netcdf(f'{work_path}/lbc_input.nc')
#