#
# This file is part of LS2D.
#
# Copyright (c) 2017-2025 Wageningen University & Research
# Author: Bart van Stratum (WUR)
#
# LS2D is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LS2D is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with LS2D.  If not, see <http://www.gnu.org/licenses/>.
#

from datetime import timedelta
import argparse
import shutil
import sys
import os

import pandas as pd
import numpy as np
import xarray as xr

# pip install ls2d
import ls2d

# pip install microhhpy
from microhhpy.real import create_input_from_regular_latlon
from microhhpy.real import create_sst_from_regular_latlon
from microhhpy.real import regrid_les, link_bcs_from_parent, link_buffer_from_parent
from microhhpy.land import create_land_surface_input, Land_surface_input
from microhhpy.thermo import calc_moist_basestate, save_basestate_density, read_basestate_density
from microhhpy.io import read_ini, check_ini, save_ini, save_case_input
from microhhpy.chem import get_rfmip_species, fit_gaussian_curve
from microhhpy.spatial import calc_vertical_grid_2nd
from microhhpy.chem import calc_tuv_photolysis
from microhhpy.utils import get_data_file
from microhhpy.logger import logger
from microhhpy.constants import xm_cams
from microhhpy.chem.emission_input import Emission_input

# Local settings and scripts.
from global_settings import sw_chemistry
from global_settings import float_type, ls2d_settings, env, outer_dom, vgrid
from global_settings import cams_eac4_variables, chemical_species, lumping_species
from corso_emissions import Corso_emissions


def ensure_binary_precision(path, shape, dtype):
    """
    Cast a raw MicroHH binary file to the requested dtype.
    """
    dtype = np.dtype(dtype)
    path = str(path)
    count = int(np.prod(shape))
    size = os.path.getsize(path)

    if size == count * dtype.itemsize:
        return

    if size == count * np.dtype(np.float32).itemsize:
        source_dtype = np.float32
    elif size == count * np.dtype(np.float64).itemsize:
        source_dtype = np.float64
    else:
        raise ValueError(
            f'Unexpected size for {path}: {size} bytes for shape {shape}')

    data = np.fromfile(path, dtype=source_dtype)
    data.astype(dtype).tofile(path)
    logger.info(f'Converted {path} from {np.dtype(source_dtype)} to {dtype}')


def read_era5_cams(ls2d_settings, start_date, end_date, cams_variables, vgrid):
    """
    Read / process ERA5 and CAMS data using (LS)2D.
    Returns the 3D fields (needed for open boundaries),
    and vertical profiles (needed for e.g. basestate).
    """
    logger.info('Reading ERA5 and CAMS using (LS)2D')

    ls2d_settings['start_date'] = start_date
    ls2d_settings['end_date'] = end_date

    era5 = ls2d.Read_era5(ls2d_settings)
    era5.calculate_forcings(n_av=0, method='2nd')
    era5_les = era5.get_les_input(vgrid.z)

    # Remove top ERA5 level, to stay within minium reference pressure RRTMGP.
    era5_les = era5_les.sel(lay=slice(0,135), lev=slice(0,136))

    # Mean profiles, only used for base-state density and dummy input model.
    era5_mean = era5_les.mean(dim='time')

    # Read CAMS data.
    cams = ls2d.Read_cams(ls2d_settings, cams_variables)
    cams.ds_ml = cams.ds_ml.rename({'go3': 'o3'})

    cams_les = cams.get_les_input(vgrid.z)

    return era5, era5_les, era5_mean, cams, cams_les


def create_nc_input(era5_1d, era5_1d_mean, cams_1d, df_tuv, domain, case_name):
    """
    Create `case_input.nc` file.
    """
    logger.info(f'Creating {case_name}_input.nc')

    # RFMIP concentrations as background species for RTE+RRTMGP.
    lon = domain.proj.central_lon
    lat = domain.proj.central_lat
    rfmip = get_rfmip_species(lat, lon, exp=0)

    # Scalar fields from CAMS. Input CAMS = mass mixing ratio (kg/kg), convert to volume mixing ratio.
    species_cams = {}
    for specie in chemical_species:
        if specie == 'co2':
            logger.warning('CO2 not available, setting to zero...')
            species_cams['co2'] = np.zeros_like(cams_1d['no2'])
        elif specie not in lumping_species:
            species_cams[specie] = cams_1d[specie].values * xm_cams['air'] / xm_cams[specie]

    # Sum lumped species as sum of converted volume mixing ratios.
    for output_specie, sub_species in lumping_species.items():
        species_cams[output_specie] = np.zeros_like(cams_1d['no2'])
        for sub_specie in sub_species:
            species_cams[output_specie] += cams_1d[sub_specie].values * xm_cams['air'] / xm_cams[sub_specie]

    eps = xm_cams['h2o'] / xm_cams['air']
    h2o = era5_1d['qt'][0,:] / (eps - eps * era5_1d['qt'][0,:])
    nudgefac = np.ones(vgrid.kmax) / 10800      # s-1

    init_profiles = {
            'z': vgrid.z,
            'thl': era5_1d['thl'][0,:],
            'qt': era5_1d['qt'][0,:],
            'u': era5_1d['u'][0,:],
            'v': era5_1d['v'][0,:],
            'o3': era5_1d['o3'][0,:]*1e-6,
            'h2o': h2o,
            'nudgefac': nudgefac}

    for name, conc in species_cams.items():
        init_profiles[name] = conc[0,:]

    radiation  = {
        'z_lay': era5_1d_mean['z_lay'  ],
        'z_lev': era5_1d_mean['z_lev'  ],
        'p_lay': era5_1d_mean['p_lay'  ],
        'p_lev': era5_1d_mean['p_lev'  ],
        't_lay': era5_1d_mean['t_lay'  ],
        't_lev': era5_1d_mean['t_lev'  ],
        'o3':    era5_1d_mean['o3_lay' ]*1e-6,
        'h2o':   era5_1d_mean['h2o_lay']}

    # NOTE: not used with heterogeneous surface, but still required by MicroHH.
    soil_index = int(era5_1d.type_soil-1)  # -1 = Fortran -> C indexing
    soil = {
            'z': era5_1d.zs[::-1],
            'theta_soil': era5_1d.theta_soil[0,::-1],
            't_soil': era5_1d.t_soil[0,::-1],
            'index_soil': np.ones(4) * soil_index,
            'root_frac': era5_1d.root_frac_low_veg[::-1]}

    for name, conc in rfmip.items():
        if name not in species_cams:
            init_profiles[name] = conc
        radiation[name] = conc

    # Photolysis rates.
    if sw_chemistry:
        time_chem = np.array(
            (df_tuv.index - df_tuv.index[0]).total_seconds(),
            dtype=float_type
        )
        emi_isop = np.zeros_like(time_chem)
        emi_no = np.zeros_like(time_chem)

        tdep_chem = {
            'time_chem': time_chem,
            'jo31d': df_tuv.jo31d,
            'jh2o2': df_tuv.jh2o2,
            'jno2': df_tuv.jno2,
            'jno3': df_tuv.jno3,
            'jn2o5': df_tuv.jn2o5,
            'jch2or': df_tuv.jch2or,
            'jch2om': df_tuv.jch2om,
            'jch3o2h': df_tuv.jch3o2h,
            'emi_isop': emi_isop,
            'emi_no': emi_no,
        }

        for name, value in tdep_chem.items():
            value = np.asarray(value, dtype=float_type)
            if np.any(~np.isfinite(value)):
                logger.warning(f'Chemistry timedep {name} contains NaNs/Infs; replacing with zero.')
            tdep_chem[name] = np.nan_to_num(value, nan=0.0, posinf=0.0, neginf=0.0)
    else:
        tdep_chem = None

    # Time dependent emissions niet nu
    tdep_source = None

    # Large-scale forcings and inflow.
    tdep_ls = {
        'time_ls': era5_1d.time_sec,
        'u_geo' : era5_1d.ug,
        'v_geo' : era5_1d.vg,
        'thl_nudge' : era5_1d.thl,
        'qt_nudge' : era5_1d.qt,
        'u_nudge' : era5_1d.u,
        'v_nudge' : era5_1d.v,
        'thl_ls' : era5_1d.dtthl_advec,
        'qt_ls' : era5_1d.dtqt_advec,
        'u_ls' : era5_1d.dtu_advec,
        'v_ls' : era5_1d.dtv_advec,
        'w_ls' : era5_1d.wls}

    for name, conc in species_cams.items():
        tdep_ls[f'{name}_nudge'] = conc
        tdep_ls[f'{name}_inflow'] = conc

    # Save in NetCDF format.
    save_case_input(
        case_name = case_name,
        init_profiles = init_profiles,
        radiation = radiation,
        soil = soil,
        tdep_ls = tdep_ls,
        tdep_chem = tdep_chem,
        tdep_source = tdep_source,
        output_dir = domain.work_dir)


def create_ini(domain, era5_1d, species, fields_emis, source_ktot, case_name):
    """
    Read base .ini file and fill in details.
    """
    logger.info(f'Creating {case_name}.ini')

    child = domain.child
    parent = domain.parent

    ini = read_ini(f'{case_name}.ini.base')

    ini['master']['npx'] = domain.npx
    ini['master']['npy'] = domain.npy

    ini['grid']['itot'] = domain.itot
    ini['grid']['jtot'] = domain.jtot
    ini['grid']['ktot'] = vgrid.kmax

    ini['grid']['xsize'] = domain.xsize
    ini['grid']['ysize'] = domain.ysize
    ini['grid']['zsize'] = vgrid.zsize

    ini['grid']['lat'] = domain.proj.central_lat
    ini['grid']['lon'] = domain.proj.central_lon

    ini['buffer']['zstart'] = 0.75 * vgrid.zsize

    ini['boundary']['scalar_outflow'] = species

    ini['chemistry']['swchemistry'] = sw_chemistry
    ini['deposition']['swdeposition'] = sw_chemistry

    ini['force']['fc'] = era5_1d.fc
    ini['force']['nudgelist'] = ['thl', 'qt', 'u', 'v'] + species
    ini['force']['timedeplist_nudge'] = ['thl', 'qt', 'u', 'v'] + species

    ini['fields']['slist'] = species
    ini['advec']['fluxlimit_list'] = ['qt'] + species
    ini['limiter']['limitlist'] = ['qt'] + species

    ini['time']['endtime'] = (domain.end_date - domain.start_date).total_seconds()
    ini['time']['datetime_utc'] = domain.start_date.strftime('%Y-%m-%d %H:%M:%S')

    ini['source']['swsource'] = '3d'
    ini['source']['ktot'] = int(source_ktot)
    ini['source']['swtimedep'] = False
    ini['source']['swheat'] = False
    ini['source']['sourcelist'] = fields_emis

    for k in [
        'source_x0', 'source_y0', 'source_z0',
        'sigma_x', 'sigma_y', 'sigma_z',
        'strength', 'swvmr',
        'line_x', 'line_y', 'line_z'
    ]:
        if k in ini['source']:
            del ini['source'][k]

    chem_vars = ['no', 'no2', 'o3', 'co', 'co2']
    path_vars = ['no_path', 'no2_path', 'o3_path', 'co_path', 'co2_path']
    ini['cross']['crosslist'] += chem_vars + path_vars

    ini['cross']['xz'] = domain.ysize/2
    ini['cross']['yz'] = domain.xsize/2

    lat = 52.48
    lon = 4.60
    x, y = domain.proj.to_xy(lon, lat)
    ini['column']['coordinates[x]'] = x
    ini['column']['coordinates[y]'] = y


    # Check if all None values are set.
    check_ini(ini)

    # Write to output .ini file.
    save_ini(ini, f'{domain.work_dir}/{case_name}.ini')


def create_surface_input(era5, era5_mean, domain, env):
    """
    Create land-surface (vegetation) and sea (SST) input.
    """
    logger.info(f'Creating spatial (land-) surface input')

    # Default soil depths IFS.
    z_soil = np.array([-0.035, -0.175, -0.64 , -1.945])[::-1]

    # Land-surface / vegetation properties from global LCC dataset (100 m resolution).
    lu_lcc = create_land_surface_input(
        domain.proj.lon,
        domain.proj.lat,
        z_soil,
        land_use_source='lcc_100m',
        land_use_tiff=env['lcc_path'],
        save_binaries=True,
        output_dir=domain.work_dir,
        save_netcdf=True,
        netcdf_file='lsm_input.nc')

    land_surface_2d_fields = [
        'alb_dif', 'alb_dir',
        'c_veg', 'cs_veg', 'gD',
        'index_veg', 'lai',
        'lambda_stable', 'lambda_unstable',
        'rs_soil_min', 'rs_veg_min',
        'water_mask', 'z0h', 'z0m',
    ]

    for field in land_surface_2d_fields:
        ensure_binary_precision(
            f'{domain.work_dir}/{field}.0000000',
            (domain.jtot, domain.itot),
            float_type)

    ensure_binary_precision(
        f'{domain.work_dir}/root_frac.0000000',
        (len(z_soil), domain.jtot, domain.itot),
        float_type)

    # TODO: Init soil from HiHydroSoil. For now spatially homogeneous.
    soil = Land_surface_input(
        domain.itot,
        domain.jtot,
        4,
        exclude_veg=True,
        debug=True,
        float_type=float_type
    )

    soil.theta_soil[:,:,:] = era5_mean.theta_soil.values[::-1, None, None]
    soil.t_soil[:,:,:] = era5_mean.t_soil.values[::-1, None, None]
    soil.index_soil[:,:,:] = int(era5_mean.type_soil) - 1  # FORTRAN -> C
    soil.to_binaries(path=domain.work_dir, allow_overwrite=True)

    # Create SSTs from ERA5.
    sst_les = create_sst_from_regular_latlon(
        era5.sst[0],
        era5.lons,
        era5.lats,
        domain.proj.lon,
        domain.proj.lat,
        float_type=float_type)

    if np.any(np.isnan(sst_les)):
        logger.warning('SSTs contain NaNs! Setting to 290K...')
        sst_les[:,:] = 290.

    # We don't know water temperatures over land, and the extrapolated SSTs are of course not a very accurate estimation...
    # TODO: get inland water mask from Corine/LCC and let user define water/lake temperatures?
    sst_les[sst_les < 280] = 280
    sst_les.tofile(f'{domain.work_dir}/t_bot_water.0000000')


def calc_photolysis_rates(env, domain):
    """
    Calculate photolysis rates using TUV wrapper.
    """
    name = 'microh'   # Must be exactly 6 characters!

    # Default input file. Only start/end date and lat/lon location are updated.
    input_file = get_data_file('microhh_tuv.base')

    tuv_df = calc_tuv_photolysis(
            input_file,
            env['tuv_path'],
            name,
            domain.start_date,
            domain.end_date,
            domain.proj.central_lon,
            domain.proj.central_lat,
            suppress_stdout=True)

    return tuv_df


### deze functie is flink aangepast, nu voor 3d emissies. Hopelijk werkt dit

def create_3d_emissions(chemical_species, domain, vgrid, no_no2_ratio=0.95,
                        n_layers_area=1, sigma_x_pt=100, sigma_y_pt=100, sigma_z_pt=100):
    np_float = float_type

    area_nc_file = "NOX+CO2_area_latlon_approx200m.nc"
    var_nox = "NOX_kg_m2_s_approx200m"
    var_co2 = "CO2_kg_m2_s_approx200m"
    stacks_excel = "puntbronnen_latlon_all.xlsx"

    dx = domain.xsize / domain.itot
    dy = domain.ysize / domain.jtot

    x = (np.arange(domain.itot) + 0.5) * dx
    y = (np.arange(domain.jtot) + 0.5) * dy
    z = vgrid.z

    if hasattr(vgrid, "dz"):
        dz_arr = np.asarray(vgrid.dz, dtype=np_float)
    else:
        z_edges = np.zeros(vgrid.kmax + 1, dtype=np_float)
        z_edges[1:-1] = 0.5 * (z[:-1] + z[1:])
        z_edges[-1] = vgrid.zsize
        dz_arr = np.diff(z_edges)

    rho_ref = np.ones(vgrid.kmax, dtype=np_float)
    fields_emis = [s for s in ["co2", "no", "no2"] if s in chemical_species]

    E3 = Emission_input(
        fields=fields_emis,
        times=[0],
        x=x, y=y, z=z,
        dz=dz_arr,
        rho_ref=rho_ref,
        float_type=np_float,
    )

    ds = xr.open_dataset(area_nc_file)

    target_lon = xr.DataArray(domain.proj.lon, dims=("y", "x"))
    target_lat = xr.DataArray(domain.proj.lat, dims=("y", "x"))

    F_nox = ds[var_nox].interp(lat=target_lat, lon=target_lon, method="nearest").values
    F_co2 = ds[var_co2].interp(lat=target_lat, lon=target_lon, method="nearest").values

    F_nox = np.nan_to_num(F_nox, nan=0.0)
    F_co2 = np.nan_to_num(F_co2, nan=0.0)

    no_kmol_m2_s  = F_nox * no_no2_ratio / xm_cams["no2"]
    no2_kmol_m2_s = F_nox * (1 - no_no2_ratio) / xm_cams["no2"]
    co2_kmol_m2_s = F_co2 / xm_cams["co2"]

    nl = int(max(1, min(n_layers_area, vgrid.kmax)))

    for kk in range(nl):
        if "no" in E3.data:
            E3.data["no"][0, kk, :, :] += (no_kmol_m2_s * xm_cams["air"] / dz_arr[kk]).astype(np_float)
        if "no2" in E3.data:
            E3.data["no2"][0, kk, :, :] += (no2_kmol_m2_s * xm_cams["air"] / dz_arr[kk]).astype(np_float)
        if "co2" in E3.data:
            E3.data["co2"][0, kk, :, :] += (co2_kmol_m2_s * xm_cams["air"] / dz_arr[kk]).astype(np_float)

    df = pd.read_excel(stacks_excel)
    df.columns = df.columns.str.strip()

    for _, row in df.iterrows():
        lat = pd.to_numeric(row["lat"], errors="coerce")
        lon = pd.to_numeric(row["lon"], errors="coerce")
        z0  = pd.to_numeric(row["height_m"], errors="coerce")
        nox_kty = pd.to_numeric(row["nox_kty"], errors="coerce")
        co2_kty = pd.to_numeric(row["co2_kty"], errors="coerce")

        if np.isnan(lat) or np.isnan(lon) or np.isnan(z0):
            continue

        x0, y0 = domain.proj.to_xy(lon, lat)

        if not (0.0 <= x0 <= domain.xsize and 0.0 <= y0 <= domain.ysize):
            continue

        if np.isfinite(nox_kty):
            nox_kg_s = nox_kty * 1e6 / (365.25 * 24 * 3600)
            no_kmol_s  = nox_kg_s * no_no2_ratio / xm_cams["no2"]
            no2_kmol_s = nox_kg_s * (1 - no_no2_ratio) / xm_cams["no2"]

            if "no" in E3.data and no_kmol_s > 0:
                E3.add_gaussian("no", no_kmol_s, 0, x0, y0, z0,
                                sigma_x_pt, sigma_y_pt, sigma_z_pt, True)
            if "no2" in E3.data and no2_kmol_s > 0:
                E3.add_gaussian("no2", no2_kmol_s, 0, x0, y0, z0,
                                sigma_x_pt, sigma_y_pt, sigma_z_pt, True)

        if np.isfinite(co2_kty) and "co2" in E3.data:
            co2_kg_s = co2_kty * 1e6 / (365.25 * 24 * 3600)
            co2_kmol_s = co2_kg_s / xm_cams["co2"]

            if co2_kmol_s > 0:
                E3.add_gaussian("co2", co2_kmol_s, 0, x0, y0, z0,
                                sigma_x_pt, sigma_y_pt, sigma_z_pt, True)

    E3.clip()
    source_ktot = int(E3.kmax)
    E3.to_binary(path=domain.work_dir)

    logger.info(f"Wrote 3D emission files with sourcelist = {fields_emis} and ktot = {source_ktot}")
    return fields_emis, source_ktot


def copy_lookup_tables(env, domain):
    """
    Copy required land-surface and radiation lookup tables.
    """
    logger.info(f'Copying lookup tables')

    microhh_path = env['microhh_path']
    gpt_path = env['gpt_path']

    rrtmgp_path = f'{microhh_path}/rte-rrtmgp-cpp/'
    rrtmgp_data_path = f'{microhh_path}/rte-rrtmgp-cpp/rrtmgp-data'

    to_copy = [
            (f'{gpt_path}/rrtmgp-gas-lw-g056-cf2.nc', 'coefficients_lw.nc'),
            (f'{gpt_path}/rrtmgp-gas-sw-g049-cf2.nc', 'coefficients_sw.nc'),
            (f'{rrtmgp_data_path}/rrtmgp-clouds-lw.nc', 'cloud_coefficients_lw.nc'),
            (f'{rrtmgp_data_path}/rrtmgp-clouds-sw.nc', 'cloud_coefficients_sw.nc'),
            (f'{rrtmgp_path}/data/aerosol_optics.nc', 'aerosol_optics.nc'),
            (f'{microhh_path}/misc/van_genuchten_parameters.nc', 'van_genuchten_parameters.nc')]

    for f in to_copy:
        target = f'{domain.work_dir}/{f[1]}'
        if not os.path.exists(target):
            shutil.copy(f[0], target)



#def main():
if True:

    # Create work directory.
    if not os.path.exists(outer_dom.work_dir):
        os.makedirs(outer_dom.work_dir)

    # Initial / boundary conditions from ERA5 / CAMS using (LS)2D.
    era5_3d, era5_1d, era5_1d_mean, cams_3d, cams_1d = read_era5_cams(
        ls2d_settings, outer_dom.start_date, outer_dom.end_date, cams_eac4_variables, vgrid)

    # Setup emissions from corso_ps_catalogue_v2.0 and corso_ps_* time/vertical profiles.
    # emissions = create_emissions(chemical_species, outer_dom, env, sigma_x=100, sigma_y=100, no_no2_ratio=0.95)

    # ---- DEBUG geostrophic wind ----
    print("fc =", era5_3d.fc if hasattr(era5_3d, "fc") else "no fc attr")
    print("ug_mean first time:", era5_3d.ug_mean[0, :20] if hasattr(era5_3d, "ug_mean") else "no ug_mean")
    print("vg_mean first time:", era5_3d.vg_mean[0, :20] if hasattr(era5_3d, "vg_mean") else "no vg_mean")

    print("era5_1d ug first time:", era5_1d["ug"][0, :20].values if "ug" in era5_1d else "no ug in era5_1d")
    print("era5_1d vg first time:", era5_1d["vg"][0, :20].values if "vg" in era5_1d else "no vg in era5_1d")

    # 3d emissions
    fields_emis, source_ktot = create_3d_emissions(
        chemical_species, outer_dom, vgrid, no_no2_ratio=0.95, n_layers_area=1, sigma_x_pt=100, sigma_y_pt=100, sigma_z_pt=100
    )

    # Calculate photolysis rates for KPP
    if sw_chemistry:
        df_tuv = calc_photolysis_rates(env, outer_dom)
    else:
        df_tuv = None

    # Create `case_input.nc` NetCDF file.
    create_nc_input(era5_1d, era5_1d_mean, cams_1d, df_tuv, outer_dom, ls2d_settings['case_name'])

    # Create `case.ini` from `case.ini.base`, filling in details.
    create_ini(outer_dom, era5_1d, chemical_species, fields_emis, source_ktot, ls2d_settings['case_name'])

    # Create land-surface (vegetation) and sea (SST) input.
    create_surface_input(era5_3d, era5_1d_mean, outer_dom, env)

    # Copy surface and radiation lookup tables.
    copy_lookup_tables(env, outer_dom)

    # Convert `case_input.nc` to float32 
    # ds = xr.open_dataset("tata_steel_nl_input.nc")

    # ds = ds.astype(np.float32)

    # ds.to_netcdf("tata_steel_nl_input.nc")

    # print("Done: converted to float32")


#if __name__ == '__main__':
#    main()
