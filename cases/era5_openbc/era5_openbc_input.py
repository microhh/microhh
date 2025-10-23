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
from pathlib import Path
import argparse
import glob

# Third-party.
import netCDF4 as nc4
import numpy as np
import ls2d

# MicroHH Python package.
from microhhpy.spatial import Domain, plot_domains, calc_vertical_grid_2nd
from microhhpy.real import create_input_from_regular_latlon
from microhhpy.real import regrid_les

import microhhpy.io as io
import microhhpy.thermo as thermo

from microhhpy.logger import logger
logger.setLevel('DEBUG')


"""
User input
"""
parser = argparse.ArgumentParser(description='Nested ERA5 input.')
parser.add_argument('-d', '--domain', type=int, required=True, help='Domain number')
args = parser.parse_args()


"""
Settings
"""
float_type = np.float64

start_date = datetime(year=2020, month=2, day=5, hour=12)
end_date   = datetime(year=2020, month=2, day=5, hour=15)

# All domains are put in a sub-folder `work_dir/domX`.
work_dir = 'test/'

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

zstart_buffer = 0.75 * gd['zsize']


"""
Define projection used for LES coordinates (m) to real world (lat/lon) transforms.
"""
# Outer domain, nested in ERA5.
dom0 = Domain(
    xsize=64_000,
    ysize=64_000,
    itot=64,
    jtot=64,
    n_ghost=3,
    n_sponge=5,
    lbc_freq=3600,                  # Always 3600 for ERA5!
    lon=settings['central_lon'],
    lat=settings['central_lat'],
    anchor='center',
    proj_str='+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs +type=crs'
)

# Inner domains(s), nested in parent LES domain.
dom1 = Domain(
    xsize = 32_000,
    ysize = 32_000,
    itot = 64,
    jtot = 64,
    n_ghost = 3,
    n_sponge = 3,
    lbc_freq = 60,
    center_in_parent=True,
    parent = dom0
)

#dom0 = Domain(
#    xsize=8_000,
#    ysize=8_000,
#    itot=8,
#    jtot=8,
#    n_ghost=3,
#    n_sponge=2,
#    lbc_freq=3600,                  # Always 3600 for ERA5!
#    lon=settings['central_lon'],
#    lat=settings['central_lat'],
#    anchor='center',
#    proj_str='+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs +type=crs'
#    )
#
## Inner domains(s), nested in parent LES domain.
#dom1 = Domain(
#    xsize = 4_000,
#    ysize = 4_000,
#    itot = 8,
#    jtot = 8,
#    n_ghost = 3,
#    n_sponge = 2,
#    lbc_freq = 60,
#    center_in_parent=True,
#    parent = dom0
#)

domains = [dom0, dom1]
for i in range(len(domains)-1):
    domains[i].child = domains[i+1]

#plot_domains([dom0, dom1], use_projection=True)

domain = domains[args.domain]
child = domain.child
parent = domain.parent
exp_dir = f'{work_dir}/dom{args.domain}/'


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
    fld = nc_group.createVariable(name, float_type, dims)
    fld[:] = data

nc_main = nc4.Dataset(f'{exp_dir}/era5_openbc_input.nc', mode='w', datamodel='NETCDF4', clobber=True)
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

ini['grid']['itot'] = domain.itot
ini['grid']['jtot'] = domain.jtot
ini['grid']['ktot'] = gd['ktot']

ini['grid']['xsize'] = domain.xsize
ini['grid']['ysize'] = domain.ysize
ini['grid']['zsize'] = gd['zsize']

ini['buffer']['zstart'] = zstart_buffer
ini['thermo']['pbot'] = ps

ini['boundary']['stop[thl]'] = stop_thl
ini['boundary']['stop[qt]' ] = stop_qt
ini['boundary']['sbot[thl]'] = sbot_thl
ini['boundary']['sbot[qt]' ] = sbot_qt

ini['time']['endtime'] = (end_date - start_date).total_seconds()

ini['cross']['xz'] = domain.ysize/2
ini['cross']['yz'] = domain.xsize/2

# Open-bounary specific settings.
ini['boundary_lateral']['n_sponge'] = domain.n_sponge
ini['boundary_lateral']['tau_sponge'] = domain.n_sponge * domain.dx / (10 * 5)
ini['boundary_lateral']['loadfreq'] = domain.lbc_freq

# Hack; how to best define this.
if args.domain == 0:
    ini['boundary_lateral']['slist'] = ['thl', 'qt']
    ini['buffer']['loadfreq'] = 3600
else:
    ini['boundary_lateral']['slist'] = ['thl', 'qt', 'qr', 'nr']
    ini['buffer']['loadfreq'] = 600

# Output LBCs for child domain.
if child is not None:
    ini['subdomain']['sw_subdomain'] = True
    ini['subdomain']['xstart'] = child.xstart_in_parent
    ini['subdomain']['ystart'] = child.ystart_in_parent
    ini['subdomain']['xend'] = child.xstart_in_parent + child.xsize 
    ini['subdomain']['yend'] = child.ystart_in_parent + child.ysize
    ini['subdomain']['grid_ratio_ij'] = int(domain.dx / child.dx)
    ini['subdomain']['grid_ratio_k'] = 1
    ini['subdomain']['n_ghost'] = child.n_ghost
    ini['subdomain']['n_sponge'] = child.n_sponge
    ini['subdomain']['savetime_bcs'] = child.lbc_freq
    ini['subdomain']['sw_save_wtop'] = True
    ini['subdomain']['sw_save_buffer'] = True
    ini['subdomain']['savetime_buffer'] = 600
    ini['subdomain']['zstart_buffer'] = zstart_buffer

else:
    ini['subdomain']['sw_subdomain'] = False
    ini['subdomain']['sw_save_wtop'] = False
    ini['subdomain']['sw_save_buffer'] = False

# Check if all values are set.
if io.check_ini(ini):
    raise Exception('Some ini values are None!')

io.save_ini(ini, f'{exp_dir}/era5_openbc.ini')


"""
Calculate and save base state density.
"""
bs = thermo.calc_moist_basestate(
    era5_1d['thl'][0,:].values,
    era5_1d['qt'][0,:].values,
    ps,
    gd['z'],
    gd['zsize'],
    dtype=float_type)

# Only save the density part for the dynamic core.
thermo.save_basestate_density(bs['rho'], bs['rhoh'], f'{exp_dir}/rhoref_overwrite.0000000')


if args.domain == 0:
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

    p_era = era5.p[:,:,:,:]
    z_era = era5.z[:,:,:,:]
    time_era = era5.time_sec

    # Standard dev. of Gaussian filter applied to interpolated fields (m).
    sigma_h = 10_000

    create_input_from_regular_latlon(
        fields_era,
        era5.lons.data,   # Strip off array masks.
        era5.lats.data,
        z_era,
        p_era,
        time_era,
        gd['z'],
        gd['zsize'],
        0.75*gd['zsize'],
        bs['rho'],
        bs['rhoh'],
        dom0,
        sigma_h,
        perturb_size=4,
        perturb_amplitude={'thl': 0.1, 'qt': 0.1e-3},
        name_suffix='overwrite',
        output_dir=exp_dir,
        ntasks=8,
        #ntasks=1,
        float_type=float_type)

else:
    """
    Regrid initial fields, and copy boundary conditions from parent domain.
    To prevent file name issues with in- and output, the simulations save
    output with `_out` appended. These need to be renamed...
    """
    parent_exp_dir = f'{work_dir}/dom{args.domain-1}/'

    # Regrid t=0 restart files.
    fields_3d = {
            'u': 0, 
            'v': 0,
            'w': 0,
            'thl': 0,
            'qt': 0,
            'qr': 0,
            'nr': 0}

    # Regrid all 2D pressure @ TOD files.
    fields_2d = {
            'phydro_tod': '*'}

    regrid_les(
            fields_3d,
            fields_2d,
            parent.xsize,
            parent.ysize,
            gd['z'],
            gd['zh'],
            parent.itot,
            parent.jtot,
            domain.xsize,
            domain.ysize,
            gd['z'],
            gd['zh'],
            domain.itot,
            domain.jtot,
            domain.xstart_in_parent,
            domain.ystart_in_parent,
            parent_exp_dir,
            exp_dir,
            float_type,
            method='nearest',
            name_suffix='overwrite')


    # Link boundary conditions from parent to child domain.
    # Only link the `_out` files, without `_out` they are LBCs used as input for the parent domain.
    def link_files(src_pattern):

        files = glob.glob(src_pattern)

        for f in files:
            src = Path(f).resolve()
            name = src.name.replace('_out', '')
            dst = Path(exp_dir) / name
            dst.symlink_to(src)

    link_files(f'{parent_exp_dir}/lbc_*_out.*')
    link_files(f'{parent_exp_dir}/w_top_out.*')
    link_files(f'{parent_exp_dir}/*_buffer_out.*')
