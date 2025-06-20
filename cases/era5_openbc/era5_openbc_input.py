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
import ls2d

sys.path.append('/home/bart/meteo/models/microhhpy')
from microhhpy.spatial import Domain, plot_domains, Vertical_grid_2nd
from microhhpy.thermo import Basestate_moist
from microhhpy.openbc import create_era5_input

from microhhpy.logger import logger
logger.setLevel('DEBUG')


"""
Settings
"""
start_date = datetime(year=2020, month=2, day=5, hour=12)
end_date   = datetime(year=2020, month=2, day=5, hour=16)
work_dir = 'test/'
TF = np.float64

n_ghost = 3     # Ghost cells used by MicroHH (advection dependent!)
n_sponge = 5    # Sponge cells used at the lateral boundaries.


"""
Read ERA5 data, for now (for simplicity) with (LS)2D.
"""
settings = {
    'start_date'  : start_date,
    'end_date'    : end_date,
    'central_lon' : -59.4,
    'central_lat' : 13.16,
    'area_size'   : 5,      # +/-, in degrees
    'case_name'   : 'barbados',
    'era5_path'   : '/home/scratch1/bart/LS2D_ERA5/',
    'cdsapirc'    : '/home/bart/.cdsapirc_ads',
    'era5_expver' : 1,
    'data_source' : 'CDS',
    'write_log'   : False
    }

ls2d.download_era5(settings)
era5 = ls2d.Read_era5(settings)
era5.calculate_forcings(n_av=3, method='2nd')



"""
Vertical grid definition.

NOTE: vertical grid definition in (LS)2D is not identical to MicroHH's grid.
      Use `microhhpy.spatial.Vertical_grid_2nd()` to make sure that the grid 
      matches the one from MicroHH, otherwise the initial fields won't
      be divergence free.
"""
_g = ls2d.grid.Grid_linear_stretched(kmax=128, dz0=20, alpha=0.01)
vgrid = Vertical_grid_2nd(_g.z, _g.zsize, remove_ghost=True)


"""
Define projection used for LES coordinates (m) to real world (lat/lon) transforms.
"""
dom0 = Domain(
    xsize=25_600,
    ysize=25_600,
    itot=64,
    jtot=64,
    n_ghost=3,
    n_sponge=5,
    lon=settings['central_lon'],
    lat=settings['central_lat'],
    anchor='center',
    proj_str='+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs +type=crs'
    )

plot_domains([dom0], use_projection=True)


"""
Calculate and save base state density.
"""
era5_1d = era5.get_les_input(vgrid.z)

bs = Basestate_moist(
    era5_1d['thl'][0,:].values,
    era5_1d['qt'][0,:].values,
    era5_1d['ps'][0].values,
    vgrid.z,
    vgrid.zsize,
    remove_ghost=True,
    dtype=TF)
    
bs.to_binary(f'{work_dir}/basestate.0000000')


"""
Create initial fields and boundary conditions, tri-linearly interpolated from ERA5.
The horizontal velocity fields are corrected to match the horizontal divergence between ERA5 and LES.
This is needed to account for interpolation errors and differences in 3D ERA5 density and the 1D LES base state density.
"""
fields_era = {
    'u': era5.u[:,:,:,:],
    'v': era5.v[:,:,:,:],
    'w': era5.wls[:,:,:,:],
    'thl': era5.thl[:,:,:,:],
    'qt': era5.qt[:,:,:,:],
}

z_era = era5.z[:,:,:,:]
time_era = era5.time_sec

sigma_h = 10_000

create_era5_input(
    fields_era,
    era5.lons.data,   # Strip off array masks.
    era5.lats.data,
    z_era,
    time_era,
    vgrid.z,
    vgrid.zsize,
    bs.rho,
    bs.rhoh,
    dom0,
    sigma_h,
    name_suffix='era5',
    output_dir=work_dir,
    ntasks=16,
    dtype=TF)
