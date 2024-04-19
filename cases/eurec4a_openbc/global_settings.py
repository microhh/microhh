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

dtype = np.float64

cosmo_path = '/home/scratch2/bart/eurec4a_cosmo/'
work_path = 'outer'

start_date = datetime(year=2020, month=2, day=1, hour=0)
end_date   = datetime(year=2020, month=2, day=2, hour=0)

# Switch between domain settings in `grid_definition.py`.
#domain_name = 'develop'
domain_name = 'test'
#domain_name = 'mip'
