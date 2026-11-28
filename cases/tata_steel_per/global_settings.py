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

import sys

from datetime import datetime
import numpy as np

# pip install ls2d
import ls2d

# pip install microhhpy
from microhhpy.spatial import Domain, plot_domains
from microhhpy.utils import check_domain_decomposition

from corso_emissions import Corso_emissions


"""
Case settings.
"""
float_type = np.float32       # KPP does not support float32.

sw_chemistry = True   # Use KPP chemistry (TODO).
sw_debug = False      # Debug mini domain.


"""
Environment settings.
"""
# env_eddy = {
#     'era5_path': '/home/scratch1/bart/LS2D_ERA5',
#     'cams_path': '/home/scratch1/bart/LS2D_CAMS',
#     'lcc_path': '/home/scratch1/bart/LCC/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif',
#     'microhh_path': '/home/bart/meteo/models/microhh',
#     'tuv_path': '/home/bart/meteo/models/microhhpy/external/TUV/V5.4',
#     'gpt_path': '/home/bart/meteo/models/coefficients_veerman',
#     'cdsapirc': '/home/bart/.cdsapirc',
#     'corso_path': '/home/scratch1/bart/emissions/corso',
#     'work_path': 'test'
# }

env_eddy = {
    'era5_path': '/Users/xiaoxing/Documents/mhh/mhh_chem/microhh/LS2D_ERA5/',
    'cams_path': '/Users/xiaoxing/Documents/mhh/mhh_chem/microhh/LS2D_CAMS/',
    'lcc_path': '/Users/xiaoxing/Documents/mhh/mhh_chem/microhh/LCC/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif',
    'microhh_path': '/Users/xiaoxing/Documents/mhh/mhh_chem/microhh',
    'tuv_path': '/Users/xiaoxing/Documents/mhh/mhh_chem/microhh//TUV/V5.4',
    'gpt_path': '/Users/xiaoxing/Documents/mhh/mhh_chem/microhh/Coefficients_Veerman/',
    'cdsapirc': '/Users/xiaoxing/.cdsapirc',
    'corso_path': '/Users/xiaoxing/Documents/mhh/mhh_chem/LS2D_microhh_cases-main/corso/corso/',
    'work_path': '.'
}


env_snellius = {
    'era5_path': '/gpfs/home4/jinma/tata/LS2D_ERA5',
    'cams_path': '/gpfs/home4/jinma/tata/LS2D_CAMS',
    'lcc_path': '/gpfs/home4/jinma/tata/LCC/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif',
    'microhh_path': '/gpfs/scratch1/shared/jinma/maarten/microhh_gpu/microhh',
    'tuv_path': '/gpfs/home4/jinma/tata/TUV/V5.4',
    'gpt_path': '/gpfs/home4/jinma/tata/Coefficients_Veerman',
    'cdsapirc': '/home/bstratum/.cdsapirc',
    'corso_path': '/gpfs/home4/jinma/tata/corso',
    'work_path': '.'
}

env = env_snellius


"""
LS2D settings.
"""
ls2d_settings = {
    'case_name'   : 'tata_steel_nl',
    'central_lat' : 52.48,
    'central_lon' : 4.60,
    # This start/end date is only used for the ERA/CAMS download.
    # The start/end date of each domain is set in the `Domain()` definitions below.
    'start_date'  : datetime(year=2023, month=7, day=8, hour=6),
    'end_date'    : datetime(year=2023, month=7, day=8, hour=18),
    # 'start_date'  : datetime(year=2021, month=2, day=23, hour=12),
    # 'end_date'    : datetime(year=2021, month=2, day=25, hour=12),
    'area_size'   : 3,  # Download +/- area_size degree of ERA5/CAMS data.
    'era5_path'   : env['era5_path'],
    'cams_path'   : env['cams_path'],
    'cdsapirc'    : env['cdsapirc'],
    'write_log'   : False,
    'data_source' : 'CDS',
    'era5_expver' : 1,  # 1=normal, 5=near-realtime ERA5.
    }


"""
Setup domains with projection and nesting.
"""
lon = ls2d_settings['central_lon']
lat = ls2d_settings['central_lat']
proj_str = f'+proj=lcc +lat_1={lat-1} +lat_2={lat+1} +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

if sw_debug:
    print('Using debug domain...')

    outer_dom = Domain(
        xsize = 19_200,
        ysize = 19_200,
        itot = 96,
        jtot = 96,
        lon = 4.60,
        lat = 52.48,
        anchor = 'center',
        start_date = datetime(year=2023, month=7, day=8, hour=6),
        end_date = datetime(year=2023, month=7, day=8, hour=18),
        proj_str = proj_str,
        work_dir = env['work_path']
        )

    # Cheating
    outer_dom.npx = 12
    outer_dom.npy = 16

    domains = [outer_dom]

    # Vertical grid.
    vgrid = ls2d.grid.Grid_linear_stretched(kmax=96, dz0=20, alpha=0.01)
    zstart_buffer = 0.75 * vgrid.zsize
    #vgrid.plot()

else:
    print('Using production domain...')
    outer_dom = Domain(
        xsize=67_200,
        ysize=48_000,
        itot=336,
        jtot=240,
        lon=4.60,
        lat=52.48,
        anchor='center',
        start_date = datetime(year=2023, month=7, day=8, hour=6),
        end_date = datetime(year=2023, month=7, day=8, hour=18),
        proj_str=proj_str,
        work_dir=env['work_path']
        )

    # Cheating
    outer_dom.npx = 16
    outer_dom.npy = 12

    domains = [outer_dom]

    # Vertical grid.
    vgrid = ls2d.grid.Grid_linear_stretched(kmax=96, dz0=20, alpha=0.015)
    zstart_buffer = 0.75 * vgrid.zsize


"""
Variables to download and read from CAMS.
"""
cams_eac4_variables = {
    'eac4_ml': [
         'carbon_monoxide',
         'nitrogen_monoxide',
         'nitrogen_dioxide',
         'nitric_acid',
         'hydrogen_peroxide',
         'ozone',
         'formaldehyde',
         'hydroperoxy_radical',
         'hydroxyl_radical',
         'nitrate_radical',
         'dinitrogen_pentoxide',
         'paraffins',
         'ethene',
         'olefins',
         'isoprene',
         'methanol',
         'ethanol',
         'propane',
         'propene',
         'terpenes',
         'peroxides',
         'methyl_peroxide',
         'methylperoxy_radical',
         'peroxy_acetyl_radical',
         'acetone_product',
         'temperature',
         'specific_humidity'],
    'eac4_sfc': [
        'surface_pressure'],
    }

# NOTE to self: only available <2021.
cams_egg4_variables = {
    'egg4_ml': [
        'carbon_dioxide',
        'methane',
        'temperature',
        'specific_humidity'],
    'egg4_sl': [
        'logarithm_of_surface_pressure']
    }


if sw_chemistry:
    #chemical_species = ['co', 'no', 'no2', 'hno3', 'h2o2', 'o3', 'hcho', 'ho2', 'oh', 'no3', 'n2o5', 'rooh', 'c3h6', 'ro2', 'co2']
    chemical_species = ['co', 'no', 'no2', 'o3', 'hcho', 'c3h6', 'co2']

    lumping_species = {
            'c3h6': ['par', 'c2h4', 'ole', 'c5h8', 'ch3oh', 'c2h5oh', 'c3h8', 'c3h6', 'c10h16'],
            'rooh': ['rooh', 'ch3ooh'], 
            #'ro2':  ['ch3o2', 'c2o3', 'aco2', 'ic3h7o2', 'hypropo2']}   # Last two species not available in ADS.
            'ro2':  ['ch3o2', 'c2o3', 'aco2']}
else:
    chemical_species = ['co2']
    lumping_species = {}


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    # Plot horizontal domain with emissions.
    margin = 0.01    # ~10 km

    lon_min = outer_dom.proj.lon.min() + margin
    lon_max = outer_dom.proj.lon.max() - margin

    lat_min = outer_dom.proj.lat.min() + margin
    lat_max = outer_dom.proj.lat.max() - margin

    emiss = Corso_emissions(env['corso_path'])
    emiss.filter_emissions(lon_min, lon_max, lat_min, lat_max)

    # Plot domains and emissions.
    fig, ax = plot_domains(
        domains,
        scatter_lon=emiss.df_emiss.longitude,
        scatter_lat=emiss.df_emiss.latitude,
        use_projection=True,
        osm_background=True,
        zoom_level=9,
        labels=['115.2 x 115.2 km @ 100 m'])

    #ax.scatter(86.1996, 22.7886, marker='x', s=100, label='Tata steel', transform=ccrs.PlateCarree())
    #plt.legend(loc='lower right')

    # Plot vertical grid.
    vgrid.plot()

    # Check domain decomposition.
    for d in domains:
        check_domain_decomposition(d.itot, d.jtot, vgrid.kmax, d.npx, d.npy)
