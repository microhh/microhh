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

import ls2d

settings = {
    'central_lon' : -57.7,
    'central_lat' : 13.3,
    'start_date'  : datetime(year=2020, month=2, day=1, hour=0),
    'end_date'    :datetime(year=2020, month=2, day=12, hour=0),
    'area_size'   : 3,
    'case_name'   : 'eurec4a_openbc',
    'era5_path'   : '/home/scratch1/bart/LS2D/',
    'era5_expver' : 1,
    'write_log'   : False,
    'data_source' : 'CDS',
    'ntasks'      : 1
    }

ls2d.download_era5(settings)
