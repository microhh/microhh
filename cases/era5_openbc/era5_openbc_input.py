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
import netCDF4 as nc4
import numpy as np
import ls2d

sys.path.append('/home/bart/meteo/models/microhhpy')
from microhhpy.spatial import Domain, plot_domains, calc_vertical_grid_2nd
from microhhpy.openbc import create_era5_input

import microhhpy.io as io
import microhhpy.thermo as thermo

from microhhpy.logger import logger
logger.setLevel('DEBUG')


"""
Settings
"""
start_date = datetime(year=2020, month=2, day=5, hour=12)
end_date   = datetime(year=2020, month=2, day=5, hour=16)
work_dir = 'test/'
TF = np.float64

sw_openbc = True

n_ghost = 3     # Ghost cells used by MicroHH (advection dependent!)
n_sponge = 5    # Sponge cells used at the lateral boundaries.


# (LS)2D settings.
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


"""
Vertical grid definition.

NOTE: vertical grid definition in (LS)2D is not identical to MicroHH's grid.
      Use `microhhpy.spatial.calculate_vertical_grid_2nd()` to make sure that the grid 
      matches the one from MicroHH, otherwise the initial fields won't be divergence free.
"""
_g = ls2d.grid.Grid_linear_stretched(kmax=128, dz0=20, alpha=0.01)
gd = calc_vertical_grid_2nd(_g.z, _g.zsize)


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
Read ERA5 data, for now (for simplicity) with (LS)2D.
"""
ls2d.download_era5(settings)
era5 = ls2d.Read_era5(settings)
era5.calculate_forcings(n_av=3, method='2nd')
era5_1d = era5.get_les_input(gd['z'])


"""
Default vertical profile input.
"""
def add_variable(name, dims, nc_group, data):
    fld = nc_group.createVariable(name, TF, dims)
    fld[:] = data

nc_main = nc4.Dataset(f'{work_dir}/era5_openbc_input.nc', mode='w', datamodel='NETCDF4', clobber=True)
nc_main.createDimension('z', gd['ktot'])
nc_init = nc_main.createGroup('init')

add_variable('z',    ('z'), nc_main, gd['z'])
add_variable('thl',  ('z'), nc_init, era5_1d.thl[0,:])
add_variable('qt',   ('z'), nc_init, era5_1d.qt[0,:])
add_variable('u',    ('z'), nc_init, era5_1d.u[0,:])
add_variable('v',    ('z'), nc_init, era5_1d.v[0,:])
add_variable('ug',   ('z'), nc_init, era5_1d.ug[0,:])
add_variable('vg',   ('z'), nc_init, era5_1d.vg[0,:])
add_variable('w_ls', ('z'), nc_init, era5_1d.wls[0,:])

nc_main.close()


"""
Calculate boundary conditions thermodynamics.
"""
thl = era5_1d.thl.mean(axis=0).values
qt  = era5_1d.qt .mean(axis=0).values

# Neumann top boundary condition.
stop_thl = (thl[-1] - thl[-2]) / (gd['z'][-1] - gd['z'][-2])
stop_qt  = (qt [-1] - qt [-2]) / (gd['z'][-1] - gd['z'][-2])

sst = float(era5_1d.sst.mean())
ps  = float(era5_1d.ps.mean())

# Dirichlet lower boundary condition.
sbot_thl = sst / thermo.exner(ps)
sbot_qt  = 0.95 * thermo.qsat(ps, sst)


"""
Create .ini file.
"""
ini = io.read_ini('era5_openbc.ini.base')

ini['grid']['itot'] = dom0.itot
ini['grid']['jtot'] = dom0.jtot
ini['grid']['ktot'] = gd['ktot']

ini['grid']['xsize'] = dom0.xsize
ini['grid']['ysize'] = dom0.ysize
ini['grid']['zsize'] = gd['zsize']

ini['buffer']['zstart'] = 0.75*gd['zsize']
ini['thermo']['pbot'] = ps

ini['boundary']['stop[thl]'] = stop_thl
ini['boundary']['stop[qt]' ] = stop_qt
ini['boundary']['sbot[thl]'] = sbot_thl
ini['boundary']['sbot[qt]' ] = sbot_qt

ini['time']['endtime'] = (end_date - start_date).total_seconds()

# Open-bounary specific settings.
ini['pres']['sw_openbc'] = sw_openbc
ini['boundary_lateral']['sw_openbc'] = sw_openbc
ini['boundary_lateral']['n_sponge'] = n_sponge


# Check if all values are set.
if io.check_ini(ini):
    raise Exception('Some ini values are None!')

io.save_ini(ini, f'{work_dir}/era5_openbc.ini')


"""
Calculate and save base state density.
"""
bs = thermo.calc_moist_basestate(
    era5_1d['thl'][0,:].values,
    era5_1d['qt'][0,:].values,
    ps,
    gd['z'],
    gd['zsize'],
    dtype=TF)

thermo.save_moist_basestate(bs, 'test/thermo_basestate_era5.0000000')

# Check with basestate generated by init phase MicroHH: ~idential.
#bs2 = thermo.read_moist_basestate('test/thermo_basestate.0000000')


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

# Standard dev. of Gaussian filter applied to interpolated fields (m).
sigma_h = 10_000

create_era5_input(
    fields_era,
    era5.lons.data,   # Strip off array masks.
    era5.lats.data,
    z_era,
    time_era,
    gd['z'],
    gd['zsize'],
    bs['rho'],
    bs['rhoh'],
    dom0,
    sigma_h,
    name_suffix='era5',
    output_dir=work_dir,
    ntasks=16,
    dtype=TF)
