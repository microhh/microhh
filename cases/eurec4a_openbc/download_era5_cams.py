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

from datetime import datetime
import sys

sys.path.append('/home/bart/meteo/models/LS2D')
import ls2d

settings = {
    'central_lon' : -57.7,
    'central_lat' : 13.3,
    'start_date'  : datetime(year=2020, month=2, day=1, hour=0),
    'end_date'    :datetime(year=2020, month=2, day=12, hour=0),
    'area_size'   : 6,
    'case_name'   : 'eurec4a_xl',
    'era5_path'   : '/home/scratch1/bart/LS2D_ERA5/',
    'cams_path'   : '/home/scratch1/bart/LS2D_CAMS/',
    'era5_expver' : 1,
    'write_log'   : False,
    'data_source' : 'CDS',
    'ntasks'      : 1,
    'cdsapirc'    : '/home/bart/.cdsapirc'
    }

# Download ERA5 meteo.
#ls2d.download_era5(settings)

# Download CAMS aerosols.
cams_vars = {
        'eac4_ml': [
            'dust_aerosol_0.03-0.55um_mixing_ratio',
            'dust_aerosol_0.55-0.9um_mixing_ratio',
            'dust_aerosol_0.9-20um_mixing_ratio',
            'hydrophilic_black_carbon_aerosol_mixing_ratio',
            'hydrophilic_organic_matter_aerosol_mixing_ratio',
            'hydrophobic_black_carbon_aerosol_mixing_ratio',
            'hydrophobic_organic_matter_aerosol_mixing_ratio',
            'sea_salt_aerosol_0.03-0.5um_mixing_ratio',
            'sea_salt_aerosol_0.5-5um_mixing_ratio',
            'sea_salt_aerosol_5-20um_mixing_ratio',
            'specific_humidity',
            'sulphate_aerosol_mixing_ratio',
            'temperature'],
        'eac4_sfc': [
            'surface_pressure'],
        }

settings['area_size'] = 2
ls2d.download_cams(settings, variables=cams_vars, grid=0.25)
