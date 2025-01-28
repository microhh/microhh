#
#  MicroHH
#  Copyright (c) 2011-2024 Chiel van Heerwaarden
#  Copyright (c) 2011-2024 Thijs Heus
#  Copyright (c) 2014-2024 Bart van Stratum
#
#  This file is part of MicroHH
#
#  MicroHH is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MicroHH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
#

# Standard library
from datetime import datetime
import sys

# Third-party.
import numpy as np
sys.path.append('/home/bart/meteo/models/LS2D')  # TMP
import ls2d

# Local library (puhhpy).
from puhhpy.spatial import Vertical_grid_2nd
from puhhpy.thermo import Basestate_moist

# Local library.
from domain_definition import domains


"""
Settings
"""
start_date = datetime(year=2022, month=5, day=1, hour=12)
end_date   = datetime(year=2022, month=5, day=1, hour=18)
work_dir = 'test/'
dtype = np.float64


"""
Grid definition.
Spatial projection is defined in `domain_definition.py`.

NOTE: vertical grid definition in (LS)2D is not identical to MicroHH's grid.
      Use `puhhpy.spatial.Vertical_grid_2nd()` to make sure that the grid 
      matches the one from MicroHH, otherwise the intial fields won't
      be divergence free.
"""
hgrid = domains[0]

_ = ls2d.grid.Grid_linear_stretched(kmax=128, dz0=20, alpha=0.01)
vgrid = Vertical_grid_2nd(_.z, _.zsize, remove_ghost=True)


"""
Read ERA5, for now with (LS)2D.
"""
settings = {
    'central_lat' : hgrid.proj.central_lat,
    'central_lon' : hgrid.proj.central_lon,
    'case_name'   : 'slocs_rf',
    'era5_path'   : '/home/scratch1/bart/LS2D_ERA5/',
    'start_date'  : start_date,
    'end_date'    : end_date,
    }

# 3D fields ERA5.
era_3d = ls2d.Read_era5(settings)

# Forcings are not used, only mean profiles to calculate base state.
era_3d.calculate_forcings(n_av=1, method='2nd')
era_1d = era_3d.get_les_input(vgrid.z)


"""
Calculate and save base state density.
NOTE: `puhhpy.thermo.Basestate_moist` includes one vertical ghost cell.
"""
bs = Basestate_moist(
    era_1d['thl'][0,:].values,
    era_1d['qt'][0,:].values,
    era_1d['ps'][0].values,
    vgrid.z,
    vgrid.zsize,
    dtype)

bs.to_binary(f'{work_dir}/basestate.0000000')