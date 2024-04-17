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

import xarray as xr
import numpy as np

## More help scripts, from (LS)2D repo.
from microhh_grid import Vertical_grid_2nd
from microhh_thermo import Basestate_moist
from global_settings import dtype, work_path

# Import correct grid settings.
import grid_definition as gd

# Calculate base state density from time averaged vertical profiles.
ds_bg = xr.open_dataset('eurec4a_mean_profiles.nc')
ds_bg = ds_bg.mean(dim='time')

# Re-produce exact vertical grid that MicroHH used internally.
uhh_gd = Vertical_grid_2nd(gd.vgrid.z, gd.vgrid.zsize)

# Calculate base state following identical procedurea as MicroHH uses.
bs = Basestate_moist(uhh_gd, ds_bg.thl.values, ds_bg.qt.values, float(ds_bg.ps))

# Save as binary file.
rho = bs.rho[uhh_gd.kstart:uhh_gd.kend]
rhoh = bs.rhoh[uhh_gd.kstart:uhh_gd.kend+1]
        
bs = np.concatenate((rho, rhoh)).astype(dtype)
bs.tofile(f'{work_path}/rhoref_0.0000000')
