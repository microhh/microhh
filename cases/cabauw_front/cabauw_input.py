import matplotlib.pyplot as pl
import netCDF4 as nc4
import xarray as xr
import numpy as np
import os, shutil
from datetime import datetime

pl.close('all')

# Available in `microhh_root/python`:
import microhh_tools as mht

class Grid_linear_stretched:
    def __init__(self, kmax, dz0, alpha):
        self.kmax = kmax
        self.dz0  = dz0
        
        self.z = np.zeros(kmax)
        self.dz = np.zeros(kmax)
        self.zsize = None

        self.dz[:] = dz0 * (1 + alpha)**np.arange(kmax)
        zh         = np.zeros(kmax+1)
        zh[1:]     = np.cumsum(self.dz)
        self.z[:]  = 0.5 * (zh[1:] + zh[:-1])
        self.zsize = zh[-1]

def add_nc_var(name, dims, nc, data):
    """
    Create NetCDF variables and set values.
    """
    if dims is None:
        var = nc.createVariable(name, np.float64)
    else:
        var = nc.createVariable(name, np.float64, dims)
    var[:] = data

def link(f1, f2):
    """
    Create symbolic link from `f1` to `f2`, if `f2` does not yet exist.
    """
    if os.path.islink(f2):
        os.remove(f2)
    if os.path.exists(f1):
        os.symlink(f1, f2)
    else:
        raise Exception('Source file {} does not exist!'.format(f1))

def copy(f1, f2):
    """
    Copy `f1` to `f2`, if `f2` does not yet exist.
    """
    if os.path.exists(f2):
        os.remove(f2)
    if os.path.exists(f1):
        shutil.copy(f1, f2)
    else:
        raise Exception('Source file {} does not exist!'.format(f1))




if __name__ == '__main__':
    """
    Case switches.
    """
    # NOTE: (LS)2D data is available from 11-08-2016 00 UTC to 12-08-2016 00 UTC.
    start = datetime(year=2016, month=8, day=11, hour=0)
    end   = datetime(year=2016, month=8, day=12, hour=0)
    ls2d_file = 'ls2d_20160811.nc'

    #start = datetime(year=2017, month=12, day=11, hour=0)
    #end   = datetime(year=2017, month=12, day=12, hour=0)
    #ls2d_file = 'ls2d_20171210.nc'


    TF = np.float64              # Switch between double (float64) and single (float32) precision.
    use_htessel = False          # False = prescribed surface H+LE fluxes from ERA5.
    use_rrtmgp = False           # False = prescribed radiation from ERA5.
    use_homogeneous_z0 = True    # False = checkerboard pattern roughness lengths.
    use_homogeneous_ls = True    # False = checkerboard pattern (some...) land-surface fields.
    single_column = True

    # Switch between the two default RRTMGP g-point sets.
    gpt_set = '128_112' # or '256_224'

    # List of scalars (for limiter)
    scalars = ['qt']
    for scalar in ['i', 'r', 's', 'g', 'h']: 
        scalars.append('q{}'.format(scalar))
        scalars.append('n{}'.format(scalar))

    # Option to link or copy the LSM + RRTMGP lookup tables.
    copy_or_link = copy

    # Link required files (if not present)
    if use_htessel:
        copy_or_link('../../misc/van_genuchten_parameters.nc', 'van_genuchten_parameters.nc')
    if use_rrtmgp:
        if gpt_set == '256_224':
            copy_or_link('../../rte-rrtmgp-cpp/rte-rrtmgp/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc', 'coefficients_lw.nc')
            copy_or_link('../../rte-rrtmgp-cpp/rte-rrtmgp/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc', 'coefficients_sw.nc')
        elif gpt_set == '128_112':
            copy_or_link('../../rte-rrtmgp-cpp/rte-rrtmgp/rrtmgp/data/rrtmgp-data-lw-g128-210809.nc', 'coefficients_lw.nc')
            copy_or_link('../../rte-rrtmgp-cpp/rte-rrtmgp/rrtmgp/data/rrtmgp-data-sw-g112-210809.nc', 'coefficients_sw.nc')
        else:
            raise Exception('\"{}\" is not a valid g-point option...'.format(gpt_set))

        copy_or_link('../../rte-rrtmgp-cpp/rte-rrtmgp/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc', 'cloud_coefficients_lw.nc')
        copy_or_link('../../rte-rrtmgp-cpp/rte-rrtmgp/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc', 'cloud_coefficients_sw.nc')

    """
    Create vertical grid for LES
    """
    vgrid = Grid_linear_stretched(160, 25, 0.014)

    """
    Read / interpolate (LS)2D initial conditions and forcings
    """
    ls2d = xr.open_dataset(ls2d_file)
    ls2d = ls2d.sel(lay=slice(0,135), lev=slice(0,136))
    ls2d = ls2d.sel(time=slice(start, end))
    ls2d['time_sec'] -= ls2d['time_sec'][0]
    ls2d_z = ls2d.interp(z=vgrid.z)

    # Reverse the soil fields. Important NOTE: in MicroHH, the vertical
    # soil index 0 is the lowest level in the soil. In (LS)2D, this
    # is reversed, and soil index 0 is the top soil level....
    # Another NOTE: the soil type in (LS)2D is the ERA5 soil type,
    # which (FORTRAN....) is 1-based, so we need to subtract 1 to
    # get the correct C-indexing.
    theta_soil = ls2d_z.theta_soil[0,::-1].values
    t_soil = ls2d_z.t_soil[0,::-1].values
    index_soil = np.ones_like(ls2d.zs)*int(ls2d.type_soil-1)
    root_frac = ls2d_z.root_frac_low_veg[::-1].values

    """
    Update .ini file
    """
    ini = mht.Read_namelist('cabauw.ini.base')

    ini['grid']['ktot'] = vgrid.kmax
    ini['grid']['zsize'] = vgrid.zsize
    ini['buffer']['zstart'] = vgrid.zsize*3/4.

    if use_htessel:
        ini['boundary']['swboundary'] = 'surface_lsm'
        ini['boundary']['sbcbot'] = 'dirichlet'
    else:
        ini['boundary']['swboundary'] = 'surface'
        ini['boundary']['sbcbot'] = 'flux'
        ini['boundary']['swtimedep'] = True
        ini['boundary']['timedeplist'] = ['thl_sbot', 'qt_sbot']

    ini['boundary']['swconstantz0'] = use_homogeneous_z0
    ini['land_surface']['swhomogeneous'] = use_homogeneous_ls

    if use_rrtmgp:
        ini['radiation']['swradiation'] = 'rrtmgp'

    elif use_htessel:
        ini['radiation']['swradiation'] = 'prescribed'
        ini['radiation']['swtimedep_prescribed'] = True
    else:
        ini['radiation']['swradiation'] = '0'

    ini['time']['endtime'] = (end-start).total_seconds()
    datetime_utc = '{0:04d}-{1:02d}-{2:02d} {3:02d}:{4:02d}:{5:02d}'.format(
            start.year, start.month, start.day, start.hour, start.minute, start.second)
    ini['time']['datetime_utc'] = datetime_utc

    ini['limiter']['limitlist'] = scalars

    if single_column:
        ini['grid']['itot'] = 1
        ini['grid']['jtot'] = 1

        ini['grid']['xsize'] = 100
        ini['grid']['ysize'] = 100

        ini['advec']['swadvec'] = 0
        ini['diff']['swdiff'] = 0
        #ini['pres']['swpres'] = 0

        ini['time']['adaptivetimestep'] = 0
        ini['time']['dtmax'] = 5


    ini.save('cabauw.ini', allow_overwrite=True)

    """
    Create MicroHH input NetCDF file.
    """
    nc = nc4.Dataset('cabauw_input.nc', mode='w', datamodel='NETCDF4')
    nc.createDimension('z', vgrid.kmax)
    add_nc_var('z', ('z'), nc, vgrid.z)

    """
    Initial profiles
    """
    nc_init = nc.createGroup('init')
    add_nc_var('thl', ('z'), nc_init, ls2d_z.thl[0,:])
    add_nc_var('qt', ('z'), nc_init, ls2d_z.qt[0,:])
    add_nc_var('u', ('z'), nc_init, ls2d_z.u[0,:])
    add_nc_var('v', ('z'), nc_init, ls2d_z.v[0,:])

    if single_column:
        add_nc_var('nudgefac', ('z'), nc_init, np.ones(vgrid.kmax)/3600)
    else:
        add_nc_var('nudgefac', ('z'), nc_init, np.ones(vgrid.kmax)/10800)

    """
    Time varying forcings
    """
    nc_tdep = nc.createGroup('timedep')
    nc_tdep.createDimension('time_surface', ls2d_z.dims['time'])
    nc_tdep.createDimension('time_ls', ls2d_z.dims['time'])

    add_nc_var('time_surface', ('time_surface'), nc_tdep, ls2d_z.time_sec)
    add_nc_var('time_ls', ('time_surface'), nc_tdep, ls2d_z.time_sec)

    add_nc_var('p_sbot', ('time_surface'), nc_tdep, ls2d_z.ps)
    add_nc_var('u_geo', ('time_ls', 'z'), nc_tdep, ls2d_z.ug)
    add_nc_var('v_geo', ('time_ls', 'z'), nc_tdep, ls2d_z.vg)

    if not use_htessel:
        add_nc_var('thl_sbot', ('time_surface'), nc_tdep, ls2d_z.wth)
        add_nc_var('qt_sbot', ('time_surface'), nc_tdep, ls2d_z.wq)

    add_nc_var('u_ls', ('time_ls', 'z'), nc_tdep, ls2d_z.dtu_advec)
    add_nc_var('v_ls', ('time_ls', 'z'), nc_tdep, ls2d_z.dtv_advec)
    add_nc_var('thl_ls', ('time_ls', 'z'), nc_tdep, ls2d_z.dtthl_advec)
    add_nc_var('qt_ls', ('time_ls', 'z'), nc_tdep, ls2d_z.dtqt_advec)
    add_nc_var('w_ls', ('time_ls', 'z'), nc_tdep, ls2d_z.wls)

    add_nc_var('u_nudge', ('time_ls', 'z'), nc_tdep, ls2d_z.u)
    add_nc_var('v_nudge', ('time_ls', 'z'), nc_tdep, ls2d_z.v)
    add_nc_var('thl_nudge', ('time_ls', 'z'), nc_tdep, ls2d_z.thl)
    add_nc_var('qt_nudge', ('time_ls', 'z'), nc_tdep, ls2d_z.qt)

    """
    Radiation variables
    """
    if use_rrtmgp:
        nc_rad = nc.createGroup('radiation')
        nc_rad.createDimension('lay', ls2d_z.dims['lay'])
        nc_rad.createDimension('lev', ls2d_z.dims['lev'])

        # Radiation variables on LES grid.
        xm_air = 28.97; xm_h2o = 18.01528
        h2o = ls2d_z.qt.mean(axis=0) * xm_air / xm_h2o
        add_nc_var('h2o', ('z'), nc_init, h2o)
        add_nc_var('o3',  ('z'), nc_init, ls2d_z.o3[0,:]*1e-6)

        # Constant concentrations:
        for group in (nc_init, nc_rad):
            add_nc_var('co2', None, group, 397e-6)
            add_nc_var('ch4', None, group, 1.8315e-6)
            add_nc_var('n2o', None, group, 3.2699e-7)
            add_nc_var('n2',  None, group, 0.781)
            add_nc_var('o2',  None, group, 0.209)

        # Radiation variables on radiation grid/levels:
        add_nc_var('z_lay', ('lay'), nc_rad, ls2d_z.z_lay.mean(axis=0))
        add_nc_var('z_lev', ('lev'), nc_rad, ls2d_z.z_lev.mean(axis=0))
        add_nc_var('p_lay', ('lay'), nc_rad, ls2d_z.p_lay.mean(axis=0))
        add_nc_var('p_lev', ('lev'), nc_rad, ls2d_z.p_lev.mean(axis=0))
        add_nc_var('t_lay', ('lay'), nc_rad, ls2d_z.t_lay.mean(axis=0))
        add_nc_var('t_lev', ('lev'), nc_rad, ls2d_z.t_lev.mean(axis=0))
        add_nc_var('o3',    ('lay'), nc_rad, ls2d_z.o3_lay.mean(axis=0)*1e-6)
        add_nc_var('h2o',   ('lay'), nc_rad, ls2d_z.h2o_lay.mean(axis=0))

    """
    Land-surface and soil
    """
    if use_htessel:
        nc_soil = nc.createGroup('soil')
        nc_soil.createDimension('z', ls2d_z.dims['zs'])
        add_nc_var('z', ('z'), nc_soil, ls2d.zs[::-1])

        add_nc_var('theta_soil', ('z'), nc_soil, theta_soil)
        add_nc_var('t_soil', ('z'), nc_soil, t_soil)
        add_nc_var('index_soil', ('z'), nc_soil, index_soil)
        add_nc_var('root_frac', ('z'), nc_soil, root_frac)

    nc.close()
