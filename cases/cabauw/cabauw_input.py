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
    Case switches.
    """
    use_htessel = True      # False = prescribed surface H+LE fluxes from ERA5.
    use_rrtmgp = False      # False = prescribed radiation from ERA5.
    use_constant_z0 = False  # False = checkerboard pattern.
    TF = np.float64         # Switch between double (float64) and single (float32) precision.

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

    if not use_rrtmgp:
        # Read ERA5 radiation, de-accumulate, and interpolate to LS2D times.
        # TODO: add to LS2D download...
        ds_rad = xr.open_dataset('era_rad_20160815_0000_2300.nc')
        ds_rad = ds_rad/3600.
        ds_rad['time'] = ds_rad['time'] - np.timedelta64(30, 'm')
        ds_rad = ds_rad.interp(time=ls2d_z.time)

    """
    Update .ini file
    """
    ini = mht.Read_namelist('cabauw.ini.base')

    ini['grid']['ktot'] = ktot
    ini['grid']['zsize'] = zsize
    ini['buffer']['zstart'] = zsize*3/4.

    if use_htessel:
        ini['boundary']['swboundary'] = 'surface_lsm'
        ini['boundary']['sbcbot'] = 'dirichlet'
    else:
        ini['boundary']['swboundary'] = 'surface'
        ini['boundary']['sbcbot'] = 'flux'
        ini['boundary']['swtimedep'] = True
        ini['boundary']['timedeplist'] = ['thl_sbot', 'qt_sbot']

    ini['boundary']['swconstantz0'] = use_constant_z0

    if use_rrtmgp:
        ini['radiation']['swradiation'] = 'rrtmgp'
    else:
        ini['radiation']['swradiation'] = 'prescribed'
        ini['radiation']['swtimedep_prescribed'] = True

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

    if not use_htessel:
        add_nc_var('thl_sbot', ('time_surface'), nc_tdep, ls2d_z.wth)
        add_nc_var('qt_sbot', ('time_surface'), nc_tdep, ls2d_z.wq)

    if not use_rrtmgp:
        add_nc_var('sw_flux_dn', ('time_surface'), nc_tdep, ds_rad.ssrd)
        add_nc_var('sw_flux_up', ('time_surface'), nc_tdep, ds_rad.ssrd-ds_rad.ssr)
        add_nc_var('lw_flux_dn', ('time_surface'), nc_tdep, ds_rad.strd)
        add_nc_var('lw_flux_up', ('time_surface'), nc_tdep, ds_rad.strd-ds_rad.str)

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
        soil_index = int(ls2d.type_soil-1) # -1 = Fortran -> C indexing.

        add_nc_var('theta_soil', ('z'), nc_soil, ls2d_z.theta_soil[0,::-1])
        add_nc_var('t_soil', ('z'), nc_soil, ls2d_z.t_soil[0,::-1])
        add_nc_var('index_soil', ('z'), nc_soil, np.ones_like(ls2d.zs)*soil_index)
        add_nc_var('root_frac', ('z'), nc_soil, ls2d_z.root_frac_low_veg[::-1])

    nc.close()

    """
    Create 2D binary input files (if needed)
    """
    if not use_constant_z0:
        """
        Create checkerboard pattern for z0m and z0h
        """
        itot = ini['grid']['itot']
        jtot = ini['grid']['jtot']

        z0m = ini['boundary']['z0m']
        z0h = ini['boundary']['z0h']

        z0m_2d = np.zeros((jtot, itot), dtype=TF)
        z0h_2d = np.zeros((jtot, itot), dtype=TF)

        blocksize_i = 8
        blocksize_j = 8

        for j in range(jtot):
            for i in range(itot):
                patch_i = i // blocksize_i % 2 == 0
                patch_j = j // blocksize_j % 2 == 0

                if (patch_i and patch_j) or (not patch_i and not patch_j):
                    z0m_2d[j,i] = z0m
                    z0h_2d[j,i] = z0h
                else:
                    z0m_2d[j,i] = 2*z0m
                    z0h_2d[j,i] = 2*z0h

        pl.figure()
        pl.imshow(z0m_2d)

        z0m_2d.tofile('z0m.0000000')
        z0h_2d.tofile('z0h.0000000')
