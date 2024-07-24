from datetime import datetime
import netCDF4 as nc4
import xarray as xr
import numpy as np
import os, shutil

# Available in `microhh_root/python`:
import microhh_tools as mht


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


def check_time_bounds(ds, start_date, end_date):
    """
    Check if start and end dates are withing Dataset bounds.
    """
    if start_date.minute != 0:
        raise Exception('Simulation has to start/end at a full hour!')

    if np.datetime64(start_date) < ds.time[0] or np.datetime64(end_date) > ds.time[-1]:
        raise Exception(f'Start or end date is out-of-bounds. Limits are {ds.time[0].values} to {ds.time[-1].values}')


copy_or_link = copy

def create_case_input(
        start_date,
        end_date,
        use_htessel,
        use_rrtmgp,
        use_rt,
        use_aerosols,
        use_tdep_aerosols,
        use_tdep_gasses,
        use_tdep_background,
        use_homogeneous_z0,
        use_homogeneous_ls,
        gpt_set,
        itot, jtot, ktot,
        xsize, ysize, zsize,
        TF,
        npx=1, npy=1):

    # Link required files (if not present)
    if use_htessel:
        copy_or_link('../../misc/van_genuchten_parameters.nc', 'van_genuchten_parameters.nc')
    if use_rrtmgp:
        if gpt_set == '256_224':
            copy_or_link('../../rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-gas-lw-g256.nc', 'coefficients_lw.nc')
            copy_or_link('../../rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-gas-sw-g224.nc', 'coefficients_sw.nc')
        elif gpt_set == '128_112':
            copy_or_link('../../rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-gas-lw-g128.nc', 'coefficients_lw.nc')
            copy_or_link('../../rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-gas-sw-g112.nc', 'coefficients_sw.nc')
        else:
            raise Exception('\"{}\" is not a valid g-point option...'.format(gpt_set))

        copy_or_link('../../rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-clouds-lw.nc', 'cloud_coefficients_lw.nc')
        copy_or_link('../../rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-clouds-sw.nc', 'cloud_coefficients_sw.nc')

    if use_aerosols:
        copy_or_link('../../rte-rrtmgp-cpp/data/aerosol_optics.nc', 'aerosol_optics.nc')

    """
    Create vertical grid for LES
    """
    dz = zsize/ktot
    z = np.arange(dz/2, zsize, dz)

    """
    Read / interpolate (LS)2D initial conditions and forcings
    """
    ls2d = xr.open_dataset('ls2d_20160815.nc')
    check_time_bounds(ls2d, start_date, end_date)

    # Remove top level ERA5 to stay within RRTMGP radiation bounds,
    # select requested time period, and interpolate to LES levels.
    ls2d = ls2d.sel(lay=slice(0,135), lev=slice(0,136), time=slice(start_date, end_date))
    ls2d_z = ls2d.interp(z=z)

    # Subtract start time.
    ls2d_z['time_sec'] = ls2d_z['time_sec'] - ls2d_z['time_sec'][0]

    # Read CAMS for aerosols and gasses other than ozone.
    cams = xr.open_dataset('cams_20160815.nc')
    check_time_bounds(cams, start_date, end_date)

    # Remove top level CAMS to stay within RRTMGP pressure bounds.
    cams = cams.sel(lay=slice(0, 135))

    # Interpolate to LES levels and ERA5 time (CAMS is 3-hourly).
    cams = cams.interp(time=ls2d.time)
    cams_z = cams.interp(z=z)

    if not use_rrtmgp:
        # Read ERA5 radiation, de-accumulate, and interpolate to LS2D times.
        # TODO: add to LS2D download...
        ds_rad = xr.open_dataset('era_rad_20160815.nc')
        ds_rad = ds_rad/3600.
        ds_rad['time'] = ds_rad['time'] - np.timedelta64(30, 'm')
        ds_rad = ds_rad.interp(time=ls2d_z.time)

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

    ini['master']['npx'] = npx
    ini['master']['npy'] = npy

    ini['grid']['itot'] = itot
    ini['grid']['jtot'] = jtot
    ini['grid']['ktot'] = ktot

    ini['grid']['xsize'] = xsize
    ini['grid']['ysize'] = ysize
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

    ini['boundary']['swconstantz0'] = use_homogeneous_z0
    ini['land_surface']['swhomogeneous'] = use_homogeneous_ls

    if use_rrtmgp and not use_rt:
        ini['radiation']['swradiation'] = 'rrtmgp'
    elif use_rrtmgp and use_rt:
        ini['radiation']['swradiation'] = 'rrtmgp_rt'
        ini['radiation']['rays_per_pixel'] = 256
        ini['radiation']['kngrid_i'] = 64
        ini['radiation']['kngrid_j'] = 64
        ini['radiation']['kngrid_k'] = 32
    else:
        ini['radiation']['swradiation'] = 'prescribed'
        ini['radiation']['swtimedep_prescribed'] = True

    ini['radiation']['swtimedep_background'] = use_tdep_background
    if use_tdep_gasses:
        ini['radiation']['timedeplist_gas'] = ['o3', 'co2', 'ch4']

    ini['aerosol']['swaerosol'] = use_aerosols
    ini['aerosol']['swtimedep'] = use_tdep_aerosols

    ini['time']['endtime'] = (end_date - start_date).total_seconds()
    d = start_date
    ini['time']['datetime_utc'] = f'{d.year}-{d.month:02d}-{d.day:02d} {d.hour:02d}:{d.minute:02d}:{d.second:02d}'

    ini.save('cabauw.ini', allow_overwrite=True)

    """
    Create MicroHH input NetCDF file.
    """
    def add_nc_var(name, dims, nc, data):
        """
        Add NetCDF variable to `nc` file or group.
        """
        if name not in nc.variables:
            if dims is None:
                var = nc.createVariable(name, np.float64)
            else:
                var = nc.createVariable(name, np.float64, dims)
            var[:] = data

    def add_nc_dim(name, size, nc):
        """
        Add NetCDF dimension, if it does not already exist.
        """
        if name not in nc.dimensions:
            nc.createDimension(name, size)

    nc = nc4.Dataset('cabauw_input.nc', mode='w', datamodel='NETCDF4')
    add_nc_dim('z', ktot, nc)
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
    add_nc_dim('time_surface', ls2d_z.dims['time'], nc_tdep)
    add_nc_dim('time_ls', ls2d_z.dims['time'], nc_tdep)

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
        add_nc_dim('lay', ls2d_z.dims['lay'], nc_rad)
        add_nc_dim('lev', ls2d_z.dims['lev'], nc_rad)

        # Radiation variables on LES grid.
        xm_air = 28.97; xm_h2o = 18.01528; eps = xm_h2o / xm_air
        qt_mean = ls2d_z.qt.mean(axis=0)
        h2o = qt_mean / (eps - eps * qt_mean)
        add_nc_var('h2o', ('z'), nc_init, h2o)
        add_nc_var('o3',  ('z'), nc_init, ls2d_z.o3[0,:]*1e-6)
        add_nc_var('co2', ('z'), nc_init, cams_z.co2[0,:]*1e-6)
        add_nc_var('ch4', ('z'), nc_init, cams_z.ch4[0,:]*1e-6)

        # Constant concentrations:
        for group in (nc_init, nc_rad):
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
        add_nc_var('h2o',   ('lay'), nc_rad, ls2d_z.h2o_lay.mean(axis=0))
        add_nc_var('o3',    ('lay'), nc_rad, ls2d_z.o3_lay.mean(axis=0)*1e-6)
        add_nc_var('co2',   ('lay'), nc_rad, cams_z.co2_lay.mean(axis=0))
        add_nc_var('ch4',   ('lay'), nc_rad, cams_z.ch4_lay.mean(axis=0))

        if use_tdep_background or use_tdep_aerosols or use_tdep_gasses:
            # NOTE: bit cheap, but ERA and CAMS are at the same time period/interval here.
            add_nc_dim('time_rad', ls2d_z.dims['time'], nc_tdep)
            add_nc_var('time_rad', ('time_rad'), nc_tdep, ls2d_z.time_sec)

            add_nc_dim('lay', ls2d_z.dims['lay'], nc_tdep)
            add_nc_dim('lev', ls2d_z.dims['lev'], nc_tdep)

        if use_tdep_gasses:
            add_nc_var('o3',  ('time_rad', 'z'), nc_tdep, ls2d_z.o3*1e-6)
            add_nc_var('co2', ('time_rad', 'z'), nc_tdep, cams_z.co2)
            add_nc_var('ch4', ('time_rad', 'z'), nc_tdep, cams_z.ch4)

        # Time dependent background profiles T, h2o, o3, ...
        if use_tdep_background:
            add_nc_var('z_lay',  ('time_rad', 'lay'), nc_tdep, ls2d_z.z_lay)
            add_nc_var('z_lev',  ('time_rad', 'lev'), nc_tdep, ls2d_z.z_lev)
            add_nc_var('p_lay',  ('time_rad', 'lay'), nc_tdep, ls2d_z.p_lay)
            add_nc_var('p_lev',  ('time_rad', 'lev'), nc_tdep, ls2d_z.p_lev)
            add_nc_var('t_lay',  ('time_rad', 'lay'), nc_tdep, ls2d_z.t_lay)
            add_nc_var('t_lev',  ('time_rad', 'lev'), nc_tdep, ls2d_z.t_lev)
            add_nc_var('h2o_bg', ('time_rad', 'lay'), nc_tdep, ls2d_z.h2o_lay)
            add_nc_var('o3_bg',  ('time_rad', 'lay'), nc_tdep, ls2d_z.o3_lay*1e-6)
            add_nc_var('co2_bg', ('time_rad', 'lay'), nc_tdep, cams_z.co2_lay)
            add_nc_var('ch4_bg', ('time_rad', 'lay'), nc_tdep, cams_z.ch4_lay)

        # Aerosols for domain and background column
        if use_aerosols:
            add_nc_var('aermr01', ('z'), nc_init, cams_z.aermr01.mean(axis=0))
            add_nc_var('aermr02', ('z'), nc_init, cams_z.aermr02.mean(axis=0))
            add_nc_var('aermr03', ('z'), nc_init, cams_z.aermr03.mean(axis=0))
            add_nc_var('aermr04', ('z'), nc_init, cams_z.aermr04.mean(axis=0))
            add_nc_var('aermr05', ('z'), nc_init, cams_z.aermr05.mean(axis=0))
            add_nc_var('aermr06', ('z'), nc_init, cams_z.aermr06.mean(axis=0))
            add_nc_var('aermr07', ('z'), nc_init, cams_z.aermr07.mean(axis=0))
            add_nc_var('aermr08', ('z'), nc_init, cams_z.aermr08.mean(axis=0))
            add_nc_var('aermr09', ('z'), nc_init, cams_z.aermr09.mean(axis=0))
            add_nc_var('aermr10', ('z'), nc_init, cams_z.aermr10.mean(axis=0))
            add_nc_var('aermr11', ('z'), nc_init, cams_z.aermr11.mean(axis=0))

            add_nc_var('aermr01', ('lay'), nc_rad, cams_z.aermr01_lay.mean(axis=0))
            add_nc_var('aermr02', ('lay'), nc_rad, cams_z.aermr02_lay.mean(axis=0))
            add_nc_var('aermr03', ('lay'), nc_rad, cams_z.aermr03_lay.mean(axis=0))
            add_nc_var('aermr04', ('lay'), nc_rad, cams_z.aermr04_lay.mean(axis=0))
            add_nc_var('aermr05', ('lay'), nc_rad, cams_z.aermr05_lay.mean(axis=0))
            add_nc_var('aermr06', ('lay'), nc_rad, cams_z.aermr06_lay.mean(axis=0))
            add_nc_var('aermr07', ('lay'), nc_rad, cams_z.aermr07_lay.mean(axis=0))
            add_nc_var('aermr08', ('lay'), nc_rad, cams_z.aermr08_lay.mean(axis=0))
            add_nc_var('aermr09', ('lay'), nc_rad, cams_z.aermr09_lay.mean(axis=0))
            add_nc_var('aermr10', ('lay'), nc_rad, cams_z.aermr10_lay.mean(axis=0))
            add_nc_var('aermr11', ('lay'), nc_rad, cams_z.aermr11_lay.mean(axis=0))

            if use_tdep_aerosols:
                add_nc_var('aermr01_bg', ('time_rad', 'lay'), nc_tdep, cams_z.aermr01_lay)
                add_nc_var('aermr02_bg', ('time_rad', 'lay'), nc_tdep, cams_z.aermr02_lay)
                add_nc_var('aermr03_bg', ('time_rad', 'lay'), nc_tdep, cams_z.aermr03_lay)
                add_nc_var('aermr04_bg', ('time_rad', 'lay'), nc_tdep, cams_z.aermr04_lay)
                add_nc_var('aermr05_bg', ('time_rad', 'lay'), nc_tdep, cams_z.aermr05_lay)
                add_nc_var('aermr06_bg', ('time_rad', 'lay'), nc_tdep, cams_z.aermr06_lay)
                add_nc_var('aermr07_bg', ('time_rad', 'lay'), nc_tdep, cams_z.aermr07_lay)
                add_nc_var('aermr08_bg', ('time_rad', 'lay'), nc_tdep, cams_z.aermr08_lay)
                add_nc_var('aermr09_bg', ('time_rad', 'lay'), nc_tdep, cams_z.aermr09_lay)
                add_nc_var('aermr10_bg', ('time_rad', 'lay'), nc_tdep, cams_z.aermr10_lay)
                add_nc_var('aermr11_bg', ('time_rad', 'lay'), nc_tdep, cams_z.aermr11_lay)

                add_nc_var('aermr01', ('time_rad', 'z'), nc_tdep, cams_z.aermr01)
                add_nc_var('aermr02', ('time_rad', 'z'), nc_tdep, cams_z.aermr02)
                add_nc_var('aermr03', ('time_rad', 'z'), nc_tdep, cams_z.aermr03)
                add_nc_var('aermr04', ('time_rad', 'z'), nc_tdep, cams_z.aermr04)
                add_nc_var('aermr05', ('time_rad', 'z'), nc_tdep, cams_z.aermr05)
                add_nc_var('aermr06', ('time_rad', 'z'), nc_tdep, cams_z.aermr06)
                add_nc_var('aermr07', ('time_rad', 'z'), nc_tdep, cams_z.aermr07)
                add_nc_var('aermr08', ('time_rad', 'z'), nc_tdep, cams_z.aermr08)
                add_nc_var('aermr09', ('time_rad', 'z'), nc_tdep, cams_z.aermr09)
                add_nc_var('aermr10', ('time_rad', 'z'), nc_tdep, cams_z.aermr10)
                add_nc_var('aermr11', ('time_rad', 'z'), nc_tdep, cams_z.aermr11)

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

    """
    Create 2D binary input files (if needed)
    """
    def get_patches(blocksize_i, blocksize_j):
        """
        Get mask for the surface patches
        """
        mask = np.zeros((jtot, itot), dtype=bool)
        mask[:] = False

        for j in range(jtot):
            for i in range(itot):
                patch_i = i // blocksize_i % 2 == 0
                patch_j = j // blocksize_j % 2 == 0

                if (patch_i and patch_j) or (not patch_i and not patch_j):
                    mask[j,i] = True

        return mask

    if not use_homogeneous_z0:
        """
        Create checkerboard pattern for z0m and z0h
        """

        z0m = ini['boundary']['z0m']
        z0h = ini['boundary']['z0h']

        z0m_2d = np.zeros((jtot, itot), dtype=TF)
        z0h_2d = np.zeros((jtot, itot), dtype=TF)

        mask = get_patches(blocksize_i=8, blocksize_j=8)

        z0m_2d[ mask] = z0m
        z0m_2d[~mask] = z0m/2.

        z0h_2d[ mask] = z0h
        z0h_2d[~mask] = z0h/2.

        z0m_2d.tofile('z0m.0000000')
        z0h_2d.tofile('z0h.0000000')

    if not use_homogeneous_ls:
        """
        Create checkerboard pattern for land-surface fields.
        """

        # Help-class to define and write the correct input for the land-surface scheme:
        # `lsm_input.py` is available in the `microhh_root/python` directory.
        from lsm_input import LSM_input

        exclude = ['z0h', 'z0m', 'water_mask', 't_bot_water']
        lsm_data = LSM_input(itot, jtot, ktot=4, TF=TF, debug=True, exclude_fields=exclude)

        # Set surface fields:
        mask = get_patches(blocksize_i=8, blocksize_j=8)

        # Patched fields:
        lsm_data.c_veg[ mask] = ini['land_surface']['c_veg']
        lsm_data.c_veg[~mask] = ini['land_surface']['c_veg']/3.

        lsm_data.lai[ mask] = ini['land_surface']['lai']
        lsm_data.lai[~mask] = ini['land_surface']['lai']/2.

        # Non-patched / homogeneous fields:
        lsm_data.gD[:,:] = ini['land_surface']['gD']
        lsm_data.rs_veg_min[:,:] = ini['land_surface']['rs_veg_min']
        lsm_data.rs_soil_min[:,:] = ini['land_surface']['rs_soil_min']
        lsm_data.lambda_stable[:,:] = ini['land_surface']['lambda_stable']
        lsm_data.lambda_unstable[:,:] = ini['land_surface']['lambda_unstable']
        lsm_data.cs_veg[:,:] = ini['land_surface']['cs_veg']

        lsm_data.t_soil[:,:,:] = t_soil[:, np.newaxis, np.newaxis]
        lsm_data.theta_soil[:,:,:] = theta_soil[:, np.newaxis, np.newaxis]
        lsm_data.index_soil[:,:,:] = index_soil[:, np.newaxis, np.newaxis]
        lsm_data.root_frac[:,:,:] = root_frac[:, np.newaxis, np.newaxis]

        # Check if all the variables have been set:
        lsm_data.check()

        # Save binary input MicroHH, and NetCDF file for visual validation/plotting/etc.
        lsm_data.save_binaries(allow_overwrite=True)
        lsm_data.save_netcdf('lsm_input.nc', allow_overwrite=True)


if __name__ == '__main__':
    """
    Case switches.
    """
    TF = np.float64              # Switch between double (float64) and single (float32) precision.
    use_htessel = True           # False = prescribed surface H+LE fluxes from ERA5.
    use_rrtmgp = True            # False = prescribed surface radiation from ERA5.
    use_rt = False               # False = 2stream solver for shortwave down, True = raytracer.
    use_homogeneous_z0 = True    # False = checkerboard pattern roughness lengths.
    use_homogeneous_ls = True    # False = checkerboard pattern (some...) land-surface fields.
    use_aerosols = False         # False = no aerosols in RRTMGP.
    use_tdep_aerosols = False    # False = time fixed RRTMGP aerosol in domain and background.
    use_tdep_gasses = False      # False = time fixed ERA5 (o3) and CAMS (co2, ch4) gasses.
    use_tdep_background = False  # False = time fixed RRTMGP T/h2o/o3 background profiles.

    """
    NOTE: `use_tdep_aerosols` and `use_tdep_gasses` specify whether the aerosols and gasses
          used by RRTMGP are updated inside the LES domain. If `use_tdep_background` is true, the
          aerosols, gasses, and the temperature & humidity are also updated on the RRTMGP background levels.
    """

    # Switch between the two default RRTMGP g-point sets.
    gpt_set = '128_112' # or '256_224'

    # Time period.
    # NOTE: Included ERA5/CAMS data is limited to 2016-08-15 06:00 - 18:00 UTC.
    start_date = datetime(year=2016, month=8, day=15, hour=6)
    end_date   = datetime(year=2016, month=8, day=15, hour=18)

    # Simple equidistant grid.
    zsize = 4000
    ktot = 160

    itot = 512
    jtot = 512

    xsize = 25600
    ysize = 25600

    # Create input files.
    create_case_input(
            start_date,
            end_date,
            use_htessel,
            use_rrtmgp,
            use_rt,
            use_aerosols,
            use_tdep_aerosols,
            use_tdep_gasses,
            use_tdep_background,
            use_homogeneous_z0,
            use_homogeneous_ls,
            gpt_set,
            itot, jtot, ktot,
            xsize, ysize, zsize,
            TF,
            npx=2,
            npy=2)
