import netCDF4 as nc4
import xarray as xr
import numpy as np
import os, shutil

# Available in `microhh_root/python`:
import microhh_tools as mht

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

copy_or_link = copy

def create_case_input(
        use_htessel,
        use_rrtmgp,
        use_rt,
        use_homogeneous_z0,
        use_homogeneous_ls,
        gpt_set,
        itot, jtot, ktot,
        xsize, ysize, zsize,
        endtime, TF):

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
    dz = zsize/ktot
    z = np.arange(dz/2, zsize, dz)

    """
    Read / interpolate (LS)2D initial conditions and forcings
    """
    ls2d = xr.open_dataset('ls2d_20160815.nc')
    ls2d = ls2d.sel(lay=slice(0,135), lev=slice(0,136))
    ls2d_z = ls2d.interp(z=z)

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

    ini['time']['endtime'] = endtime

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
    use_rrtmgp = True            # False = prescribed radiation from ERA5.
    use_rt = False               # False = 2stream solver for shortwave down, True = raytracer.
    use_homogeneous_z0 = True    # False = checkerboard pattern roughness lengths.
    use_homogeneous_ls = True    # False = checkerboard pattern (some...) land-surface fields.

    # Switch between the two default RRTMGP g-point sets.
    gpt_set = '128_112' # or '256_224'

    # Simple equidistant grid.
    zsize = 4000
    ktot = 160

    itot = 512
    jtot = 512

    xsize = 25600
    ysize = 25600

    endtime = 43200

    # Create input files.
    create_case_input(
            use_htessel,
            use_rrtmgp,
            use_rt,
            use_homogeneous_z0,
            use_homogeneous_ls,
            gpt_set,
            itot, jtot, ktot,
            xsize, ysize, zsize,
            endtime,
            TF)
