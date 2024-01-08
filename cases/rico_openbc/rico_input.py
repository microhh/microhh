import shutil
import stat
import sys
import os

import netCDF4 as nc
import xarray as xr
import numpy as np

# Custom modules, from `microhh_root/python` directory.
import microhh_lbc_tools as mlt
import microhh_tools as mht


def grid_linear_stretched(dz0, alpha, ktot):
    """
    Simple linearly stretched grid where the grid spacing
    in each level increases with a factor `1+alpha`.
    """
    dz = dz0 * (1 + alpha)**np.arange(ktot) 
    zh = np.zeros(ktot+1)
    zh[1:]= np.cumsum(dz)
    z = 0.5 * (zh[1:] + zh[:-1])
    zsize = zh[-1]

    return z, zh, zsize


def check_grid_decomposition(itot, jtot, ktot, npx, npy):
    """
    Check whether grid / MPI decomposition is valid
    """

    err = False
    if itot%npx != 0:
        print('ERROR in grid: itot%npx != 0')
        err = True

    if itot%npy != 0:
        print('ERROR in grid: itot%npy != 0')
        err = True

    if jtot%npx != 0 and npy > 1:
        print('ERROR in grid: jtot%npx != 0')
        err = True

    if jtot%npy != 0:
        print('ERROR in grid: jtot%npy != 0')
        err = True

    if ktot%npx != 0:
        print('ERROt in grid: ktot%npx != 0')
        err = True

    if err:
        print('Grid: itot={}, jtot={}, ktot={}, npx={}, npy={}'.format(
            itot, jtot, ktot, npx, npy))
        raise Exception('Invalid grid configuration!')
    else:
        print('Grid: itot={}, jtot={}, ktot={}, npx={}, npy={}: OKAY!'.format(
            itot, jtot, ktot, npx, npy))


if __name__ == '__main__':

    if len(sys.argv) != 2:
        raise Exception('Provide domain number as argument (0...N)')
    domain = int(sys.argv[1])

    """
    Case settings.
    """
    # Work directory. Each domain is placed in its own sub-directory.
    #work_path = '.'
    work_path = '/home/stratum2/scratch/rico_2i6'

    """
    # Outer domain with doubly-periodic BCs.
    d0 = mlt.Domain(
            name = 'dom_0',
            itot = 96,
            jtot = 56,
            dx = 360,
            dy = 360,
            end_time = 6*3600,
            work_path = work_path)

    # Inner domains with open BCs.
    d1 = mlt.Domain(
            name = 'dom_1',
            itot = 96,
            jtot = 56,
            dx = 120,
            dy = 120,
            center_in_parent = True,
            start_offset = 0, #7200,
            end_offset = 0,
            parent = d0,
            work_path = work_path)

    d2 = mlt.Domain(
            name = 'dom_2',
            itot = 96,
            jtot = 56,
            dx = 40,
            dy = 40,
            center_in_parent = True,
            start_offset = 0,
            end_offset = 0,
            parent = d1,
            work_path = work_path)

    d0.child = d1
    d1.child = d2
    """

    # Outer domain with doubly-periodic BCs.
    d0 = mlt.Domain(
            name = 'dom_0',
            itot = 768,
            jtot = 432,
            dx = 360,
            dy = 360,
            end_time = 72*3600,
            work_path = work_path)

    # Inner domains with open BCs.
    d1 = mlt.Domain(
            name = 'dom_1',
            itot = 1536,
            jtot = 864,
            dx = 120,
            dy = 120,
            center_in_parent = True,
            start_offset = 14400,
            end_offset = -40*3600,
            parent = d0,
            work_path = work_path)

    d2 = mlt.Domain(
            name = 'dom_2',
            itot = 1536,
            jtot = 864,
            dx = 40,
            dy = 40,
            center_in_parent = True,
            start_offset = 0,
            end_offset = 0,
            parent = d1,
            work_path = work_path)

    d0.child = d1
    d1.child = d2

    float_type = np.float64
    #microhh_path = '/home/bart/meteo/models/microhh'
    #microhh_bin = '/home/bart/meteo/models/microhh/build_dp_cpumpi/microhh'

    microhh_path = '/home/stratum2/models/microhh'
    microhh_bin = '/home/stratum2/models/microhh/build_dp_cpumpi/microhh'
    
    #case = 'gcss'  # Original RICO
    case = 'ss08' # Moist RICO from Stevens/Seifert & Seifert/Heus
    #case = 'test' # More moist mixed-layer for testing

    sw_advec = '2i5'
    sw_sponge = True
    n_ghost = 3
    n_sponge = 5 if sw_sponge else 0
    lbc_freq = 60

    """
    Generate case input.
    """
    # Select correct `Domain` instance. 
    if domain == 0:
        domain = d0
    elif domain == 1:
        domain = d1
    elif domain == 2:
        domain = d2
    else:
        raise Exception('Domain number out of range!')

    # Create work directory.
    if not os.path.exists(domain.work_dir):
        os.makedirs(domain.work_dir)

    # Read base settings.
    ini = mht.Read_namelist('rico.ini.base')
    
    # Linearly stretched vertical grid.
    ktot = 144
    z, zh, zsize = grid_linear_stretched(dz0=20, alpha=0.007, ktot=ktot)
    
    # Define fields and vertical profiles.
    thl   = np.zeros(ktot)
    qt    = np.zeros(ktot)
    u     = np.zeros(ktot)
    ugeo  = np.zeros(ktot)
    v     = np.zeros(ktot)
    vgeo  = np.zeros(ktot)
    wls   = np.zeros(ktot)
    thlls = np.zeros(ktot)
    qtls  = np.zeros(ktot)
    
    for k in range(ktot):
        # Liquid water potential temperature: same in GCSS and SS08
        if(z[k] < 740.):
            thl[k] = 297.9
        else:
            thl[k] = 297.9 + (317.0 - 297.9)/(4000. - 740.) * (z[k] - 740.) 
    
        if(case == 'gcss'):
            if(z[k] < 740.):
                qt[k] = 16.0 + (13.8 - 16.0) / 740. * z[k]
            elif(z[k] < 3260.):
                qt[k] = 13.8 + (2.4 - 13.8) / (3260. - 740.) * (z[k] - 740.) 
            else:
                qt[k] = 2.4 + (1.8 - 2.4)/(4000. - 3260.) * (z[k] - 3260.) 
    
        elif(case == 'ss08'):
            if(z[k] < 740.):
                qt[k] = 16.0 + (13.8 - 16.0) / 740. * z[k]
            elif(z[k] < 3260.):
                qt[k] = 13.8 + (4.4 - 13.8) / (3260. - 740.) * (z[k] - 740.) 
            else:
                qt[k] = 4.4 + (3.6 - 4.4)/(4000. - 3260.) * (z[k] - 3260.) 
    
        elif(case == 'test'):
            q0 = 18.
            q1 = 15.8
            if(z[k] < 740.):
                qt[k] = q0 + (q1 - q0) / 740. * z[k]
            elif(z[k] < 3260.):
                qt[k] = q1 + (2.4 - q1) / (3260. - 740.) * (z[k] - 740.) 
            else:
                qt[k] = 2.4 + (1.8 - 2.4)/(4000. - 3260.) * (z[k] - 3260.) 
    
        # Subsidence
        if(z[k] < 2260):
            wls[k] = -0.005 * (z[k] / 2260.)
        else:
            wls[k] = -0.005
    
        # U and V component wind
        u[k] = -9.9 + 2.0e-3 * z[k]
        v[k] = -3.8
        ugeo[k] = u[k]
        vgeo[k] = v[k]
    
        # Advective and radiative tendency thl
        thlls[k] = -2.5 / 86400.
    
        # Advective tendency qt
        if(z[k] < 2980):
            qtls[k] = -1.0 / 86400. + (1.3456/ 86400.) * z[k] / 2980.
        else:
            qtls[k] = 4e-6
    
    # normalize profiles to SI
    qt  /= 1000.
    qtls/= 1000.

    """
    Create NetCDF file.
    """
    def add_var(nc_group, name, dims, values):
        nc = nc_group.createVariable(name, float_type, dims)
        nc[:] = values

    # write the data to a file
    nc_file = nc.Dataset(
            f'{domain.work_dir}/rico_input.nc', mode='w',
            datamodel='NETCDF4', clobber=True)

    nc_file.createDimension('z', ktot)
    add_var(nc_file, 'z', ('z'), z)

    nc_group_init = nc_file.createGroup('init');
    add_var(nc_group_init, 'thl',    ('z'), thl)
    add_var(nc_group_init, 'qt',     ('z'), qt)
    add_var(nc_group_init, 'u',      ('z'), u)
    add_var(nc_group_init, 'u_geo',  ('z'), ugeo)
    add_var(nc_group_init, 'v',      ('z'), v)
    add_var(nc_group_init, 'v_geo',  ('z'), vgeo)
    add_var(nc_group_init, 'w_ls',   ('z'), wls)
    add_var(nc_group_init, 'thl_ls', ('z'), thlls)
    add_var(nc_group_init, 'qt_ls',  ('z'), qtls)

    nc_file.close()
    
    # Calculate SST
    ep = 287.04 / 461.5 
    
    # Surface settings
    def esat(T):
        c0 = 0.6105851e+03; c1 = 0.4440316e+02; c2 = 0.1430341e+01; c3 = 0.2641412e-01 
        c4 = 0.2995057e-03; c5 = 0.2031998e-05; c6 = 0.6936113e-08; c7 = 0.2564861e-11 
        c8 = -.3704404e-13 
        x  = max(-80.,T-273.15)
        return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
    
    def qsat(p, T):
        return ep*esat(T)/(p-(1.-ep)*esat(T))
    
    ps  = 101540.
    SST = 299.8 
    ths = SST / (ps/1.e5)**(287.04/1005.)
    qs  = qsat(ps, SST) 

    """
    Create LBCs & initial fields from parent domain.
    """
    if domain.parent is not None:
        print('Creating LBCs...')

        interpolation_method = 'linear'
        fields = ['thl', 'qt', 'qr', 'nr', 'u' ,'v', 'w']

        # Offset of current domain in parent domain.
        xstart = domain.i0_in_parent * domain.parent.dx
        ystart = domain.j0_in_parent * domain.parent.dy

        """
        Create LBCs from parent domain.
        """
        xz = {}
        yz = {}
        for fld in fields:
            xz[fld] = xr.open_dataset(f'{domain.parent.work_dir}/{fld}.xz.nc', decode_times=False)
            yz[fld] = xr.open_dataset(f'{domain.parent.work_dir}/{fld}.yz.nc', decode_times=False)

            # Select correct time period.
            t0 = domain.start_offset
            t1 = domain.parent.end_time+domain.end_offset

            xz[fld] = xz[fld].sel(time=slice(t0, t1))
            yz[fld] = yz[fld].sel(time=slice(t0, t1))

        time = xz[list(xz.keys())[0]].time.values

        # Create Xarray Dataset with correct fields/coordinates/.. of LBCs.
        time_lbc = np.arange(0, domain.end_time+1, lbc_freq)

        lbc_ds = mlt.get_lbc_xr_dataset(
                fields,
                domain.xsize,
                domain.ysize,
                domain.itot,
                domain.jtot,
                z, zh[:-1],
                time_lbc,
                n_ghost,
                n_sponge,
                dtype = float_type)

        # Add offset of child in parent domain, for easy interpolation with Xarray.
        for v in lbc_ds.variables:
            if 'x' in v:
                lbc_ds[v] = lbc_ds[v] + xstart
            if 'y' in v:
                lbc_ds[v] = lbc_ds[v] + ystart
    
        print('Interpolating LBCs...')
        for loc in ['north', 'west', 'east', 'south']:
            for fld in fields:
                print(f' - {fld}-{loc}')

                # Short cuts.
                lbc_in = lbc_ds[f'{fld}_{loc}']
                dims = lbc_in.dims
    
                # Dimensions in LBC file.
                xloc, yloc = dims[3], dims[2]
    
                # Dimensions in cross-section.
                xloc_in = 'xh' if 'xh' in xloc else 'x'
                yloc_in = 'yh' if 'yh' in yloc else 'y'
    
                # Switch between yz and xz crosses.
                cc = yz if loc in ['west','east'] else xz
    
                # Interpolate!
                ip = cc[fld].interp({yloc_in: lbc_ds[yloc], xloc_in: lbc_ds[xloc]}, method=interpolation_method)
    
                # Check if interpolation was success.
                if np.any(np.isnan(ip[fld].values)):
                    raise Exception('Interpolated BCs contain NaNs!')
    
                lbc_in[:] = ip[fld].values

        # DEBUG:
        print('Saving as NetCDF...')
        lbc_ds.to_netcdf(f'{domain.work_dir}/rico_lbc_input.nc')

        # Write binary input files for MicroHH.
        print('Saving as binaries...')
        mlt.write_dataset_as_binaries(lbc_ds, float_type, output_dir=domain.work_dir)

        """
        Interpolate initial fields from parent domain.
        """
        print('Interpolating initial fields...')

        # Lookup table with LES dimensions:
        dims = {
                'x':  np.arange(domain.dx/2, domain.xsize, domain.dx) + xstart,
                'xh': np.arange(0, domain.xsize, domain.dx) + xstart,
                'y':  np.arange(domain.dy/2, domain.ysize, domain.dy) + ystart,
                'yh': np.arange(0, domain.ysize, domain.dy) + ystart}

        for fld in fields:
            ds = xr.open_dataset(f'{domain.parent.work_dir}/{fld}.nc', decode_times=False)

            # Select start time.
            ds = ds.sel(time=domain.start_offset)

            # Interpolate to correct coordinates
            dim_y = ds[fld].dims[1]  # x or xh
            dim_x = ds[fld].dims[2]  # y or yh

            ds.interp({
                dim_x: dims[dim_x],
                dim_y: dims[dim_y]}, method=interpolation_method)

            # Save as binary file.
            ds[fld].values.astype(float_type).tofile(f'{domain.work_dir}/{fld}_0.0000000')


    """
    Calculate location of nest for cross-sections.
    """
    if domain.child is not None:
        # Approximate location child in parent:
        x0 = (domain.xsize - domain.child.xsize) / 2.
        y0 = (domain.ysize - domain.child.ysize) / 2.

        # Exact location child in parent.
        i0 = int(x0 / domain.dx)
        j0 = int(y0 / domain.dy)

        x0 = i0 * domain.dx
        y0 = j0 * domain.dy

        x1 = x0 + domain.child.xsize
        y1 = y0 + domain.child.ysize

        xz, yz = mlt.get_cross_locations_for_lbcs(
                x0, y0,
                x1, y1,
                domain.dx, domain.dy,
                domain.child.dx, domain.child.dy,
                n_ghost, n_sponge)

    """
    Create new `.ini` file.
    """
    ini['grid']['itot'] = domain.itot
    ini['grid']['jtot'] = domain.jtot
    ini['grid']['ktot'] = ktot

    ini['grid']['xsize'] = domain.xsize
    ini['grid']['ysize'] = domain.ysize
    ini['grid']['zsize'] = zsize

    ini['buffer']['zstart'] = 0.8*zsize

    ini['boundary']['sbot[thl]'] = ths
    ini['boundary']['sbot[qt]'] = qs

    ini['advec']['swadvec'] = sw_advec

    ini['time']['endtime'] = domain.end_time

    if domain.parent is None:
        # Outer domain with doubly-periodic BCs.
        ini['pres']['sw_openbc'] = False
        ini['boundary_lateral']['sw_openbc'] = False
    else:
        # Inner domain with open BCs.
        ini['pres']['sw_openbc'] = True
        ini['boundary_lateral']['sw_openbc'] = True
        ini['boundary_lateral']['sw_sponge'] = sw_sponge
        ini['boundary_lateral']['loadfreq'] = lbc_freq

    if domain.child is not None:
        ini['cross']['xz'] = list(xz)
        ini['cross']['yz'] = list(yz)

    ini['cross']['sampletime'] = lbc_freq

    ini.save(f'{domain.work_dir}/rico.ini', allow_overwrite=True)


    """
    Copy necessary help scripts.
    """
    to_copy = [    
        microhh_bin,
        f'{microhh_path}/python/microhh_tools.py',
        f'{microhh_path}/python/cross_to_nc.py',
        f'{microhh_path}/python/3d_to_nc.py']

    for f in to_copy:
        target = '{}/{}'.format(domain.work_dir, f.split('/')[-1])
        if not os.path.exists(target):
            shutil.copy(f, domain.work_dir)

#    """
#    Yes, I'm lazy: create run script.
#    """
#    runscript = f'{domain.work_dir}/init_run.sh'
#    with open(runscript, 'w') as f:
#        f.write('mpiexec -n 8 ./microhh init rico\n\n')
#
#        f.write('mv thl_0.0000000 thl.0000000\n')
#        f.write('mv qt_0.0000000 qt.0000000\n')
#        f.write('mv qr_0.0000000 qr.0000000\n')
#        f.write('mv nr_0.0000000 nr.0000000\n')
#        f.write('mv u_0.0000000 u.0000000\n')
#        f.write('mv v_0.0000000 v.0000000\n')
#        f.write('mv w_0.0000000 w.0000000\n\n')
#
#        f.write('mpiexec -n 8 ./microhh run rico\n\n')
#
#        if domain.child is not None:
#            tstart = domain.child.start_offset
#            f.write(f'python cross_to_nc.py -n 12\n')
#            f.write(f'python 3d_to_nc.py -v thl qt qr nr u v w -t0 {tstart} -t1 {tstart} -n 6')
#
#    st = os.stat(runscript)
#    os.chmod(runscript, st.st_mode | stat.S_IEXEC)
