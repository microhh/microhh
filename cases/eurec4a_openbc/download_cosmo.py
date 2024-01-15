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
import urllib.request
import os

# The simulation period is from  00 UTC on February 1st 2020 to 00 UTC on February 12th 2020.
start = datetime(year=2020, month=2, day=1,  hour=0)
end   = datetime(year=2020, month=2, day=12, hour=0)

remote_url = 'https://swift.dkrz.de/v1/dkrz_0913c8f3-e7b6-4f94-9221-06880d4ccfea'
output_dir = '/home/scratch2/bart/eurec4a_cosmo/'

date = start
while date <= end:
    print(f'Downloading {date}')

    filename_3d = f'lffd{date.year:04d}{date.month:02d}{date.day:02d}{date.hour:02d}0000z.nc'
    filename_2d = f'lffd{date.year:04d}{date.month:02d}{date.day:02d}{date.hour:02d}0000.nc'

    remote_3d = f'{remote_url}/COSMO_CTRL_BC_3D/{filename_3d}'
    remote_2d = f'{remote_url}/COSMO_CTRL_BC_2D/{filename_2d}'

    local_3d = f'{output_dir}/COSMO_CTRL_BC_3D/{filename_3d}'
    local_2d = f'{output_dir}/COSMO_CTRL_BC_2D/{filename_2d}'

    if not os.path.exists(local_3d):
        urllib.request.urlretrieve(remote_3d, local_3d) 

    if not os.path.exists(local_2d):
        urllib.request.urlretrieve(remote_2d, local_2d) 

    date += timedelta(hours=1)

