import matplotlib.pyplot as pl
import netCDF4 as nc4
import xarray as xr
import numpy as np

# Available in `microhh_root/python`:
import microhh_tools as mht

pl.close('all')

def add_nc_var(name, dims, nc, data):
    """
    Create NetCDF variables and set values.
    """
    if dims is None:
        var = nc.createVariable(name, np.float64)
    else:
        var = nc.createVariable(name, np.float64, dims)
    var[:] = data

if __name__ == '__main__':
    """
    Case switches. LSM without radiation is not supported (duhh).
    """
    use_interactive_lsm = True    # True = HTESSEL, False = prescribed surface fluxes from ERA5
    use_radiation = True          # True = RRTMGP, False = no radiation...

    """
    Create vertical grid for LES
    """
    zsize = 3200
    ktot = 128
    dz = zsize/ktot
    z = np.arange(dz/2, zsize, dz)

    """
    Read / interpolate (LS)2D initial conditions and forcings
    """
    ls2d = xr.open_dataset('ls2d_20160815_0600_1800.nc')
    ls2d = ls2d.sel(lay=slice(0,135), lev=slice(0,136))
    ls2d_z = ls2d.interp(z=z)

    """
    Update .ini file
    """
    ini = mht.Read_namelist('cabauw.ini.base')

    ini['grid']['ktot'] = ktot
    ini['grid']['zsize'] = zsize
    ini['buffer']['zstart'] = zsize*3/4.

    if use_interactive_lsm:
        ini['boundary']['swboundary'] = 'surface_lsm'
        ini['boundary']['sbcbot'] = 'dirichlet'
    else:
        ini['boundary']['swboundary'] = 'surface'
        ini['boundary']['sbcbot'] = 'flux'
        ini['boundary']['swtimedep'] = True
        ini['boundary']['timedeplist'] = ['thl_sbot', 'qt_sbot']

    if use_radiation:
        ini['radiation']['swradiation'] = 'rrtmgp'
    else:
        ini['radiation']['swradiation'] = '0'

    ini.save('cabauw.ini', allow_overwrite=True)

    """
    Create MicroHH input NetCDF file.
    """
    nc = nc4.Dataset('cabauw_input.nc', mode='w', datamodel='NETCDF4')
    nc.createDimension('z', ktot)
    add_nc_var('z', ('z'), nc, z)

    """
    Initial profiles
    """
    nc_init = nc.createGroup('init')
    add_nc_var('thl', ('z'), nc_init, ls2d_z.thl[0,:])
    add_nc_var('qt', ('z'), nc_init, ls2d_z.qt[0,:])
    add_nc_var('u', ('z'), nc_init, ls2d_z.u[0,:])
    add_nc_var('v', ('z'), nc_init, ls2d_z.v[0,:])
    add_nc_var('nudgefac', ('z'), nc_init, np.ones(ktot)/10800)

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

    if not use_interactive_lsm:
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
    if use_radiation:
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
    if use_interactive_lsm:
        nc_soil = nc.createGroup('soil')
        nc_soil.createDimension('z', ls2d_z.dims['zs'])
        add_nc_var('z', ('z'), nc_soil, ls2d.zs[::-1])
        soil_index = int(ls2d.type_soil-1) # -1 = Fortran -> C indexing.

        add_nc_var('theta_soil', ('z'), nc_soil, ls2d_z.theta_soil[0,::-1])
        add_nc_var('t_soil', ('z'), nc_soil, ls2d_z.t_soil[0,::-1])
        add_nc_var('index_soil', ('z'), nc_soil, np.ones_like(ls2d.zs)*soil_index)
        add_nc_var('root_frac', ('z'), nc_soil, ls2d_z.root_frac_low_veg[::-1])

    nc.close()
