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
import numpy as np

import helpers as hlp

float_type = np.float32

# Data paths.
syst = 'eddy'
#syst = 'snellius'

if syst == 'eddy':
    cosmo_path = '/home/scratch2/bart/eurec4a_cosmo/'
    work_path = 'develop_case'
    microhh_path = '/home/bart/meteo/models/microhh/'
    gpt_veerman_path = '/home/bart/meteo/models/coefficients_veerman'
    ls2d_era5_path = '/home/scratch1/bart/LS2D_ERA5'

elif  syst == 'snellius':
    cosmo_path = '/gpfs/work3/0/lesmodels/eurec4a'
    work_path = '/scratch-shared/stratum2/eurec4a_test/'
    microhh_path = '/home/stratum2/meteo/models/microhh'
    gpt_veerman_path = '/home/stratum2/meteo/models/coefficients_veerman'
    ls2d_era5_path = '/gpfs/work3/0/lesmodels/ls2d_era5'

start_date = datetime(year=2020, month=2, day=1, hour=0)
end_date   = datetime(year=2020, month=2, day=1, hour=4)

case = 'develop'
#case = 'test'
