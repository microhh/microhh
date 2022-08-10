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
    microhh_path = '/home/stratum2/models/microhh_rt'
    scratch_path = '/scratch-shared/stratum2/cabauw_rt'
    microhh_bin = '{}/build_sp_gpu/microhh'.format(microhh_path)

    """
    Shared settings
    """
    TF = np.float64          # Switch between double (float64) and single (float32) precision.
    gpt_set = '128_112'      # Switch between the two default RRTMGP g-point sets (256_224 or 128_112)

    # LES grid
    zsize = 4000
    ktot = 200
    dz = zsize/ktot
    z = np.arange(dz/2, zsize, dz)

    itot = 768
    jtot = 768
    xsize = itot*50
    ysize = jtot*50

    # Slurm settings
    wc_time_2s = '24:00:00'
    wc_time_rt = '48:00:00'

    def generate_case(
            use_rt,
            homogenize_sw=False, homogenize_lw=False,
            homogenize_hr_sw=False, homogenize_hr_lw=False, spppsqp=256):

        # Create unique job name
        job_name = ''

        if use_rt:
            job_name += 'rt'
            wc_time = wc_time_rt
        else:
            job_name += '2s'
            wc_time = wc_time_2s

        if homogenize_sw or homogenize_lw or homogenize_hr_sw or homogenize_hr_lw:
            job_name += '_hom'
        else:
            job_name += '_het'

        if homogenize_sw and not homogenize_lw:
            job_name += '_sw'
        elif homogenize_sw and homogenize_lw:
            job_name += '_swlw'

        if homogenize_hr_lw and homogenize_hr_sw:
            job_name += '_hr'
        elif homogenize_hr_sw:
            job_name += '_hrsw'
        elif homogenize_hr_lw:
            job_name += '_hrlw'

        if spppsqp != 256:
            job_name += '_{}'.format(spppsqp)

        # Link required files (if not present)
        copy('{}/misc/van_genuchten_parameters.nc'.format(microhh_path), 'van_genuchten_parameters.nc')

        if gpt_set == '256_224':
            copy('{}/rte-rrtmgp-cpp/rte-rrtmgp/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc'.format(microhh_path), 'coefficients_lw.nc')
            copy('{}/rte-rrtmgp-cpp/rte-rrtmgp/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc'.format(microhh_path), 'coefficients_sw.nc')
        elif gpt_set == '128_112':
            copy('{}/rte-rrtmgp-cpp/rte-rrtmgp/rrtmgp/data/rrtmgp-data-lw-g128-210809.nc'.format(microhh_path), 'coefficients_lw.nc')
            copy('{}/rte-rrtmgp-cpp/rte-rrtmgp/rrtmgp/data/rrtmgp-data-sw-g112-210809.nc'.format(microhh_path), 'coefficients_sw.nc')
        else:
            raise Exception('\"{}\" is not a valid g-point option...'.format(gpt_set))

        copy('{}/rte-rrtmgp-cpp/rte-rrtmgp/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc'.format(microhh_path), 'cloud_coefficients_lw.nc')
        copy('{}/rte-rrtmgp-cpp/rte-rrtmgp/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc'.format(microhh_path), 'cloud_coefficients_sw.nc')

        """
        Read / interpolate (LS)2D initial conditions and forcings
        """
        ls2d = xr.open_dataset('ls2d_20160815.nc')
        ls2d = ls2d.sel(lay=slice(0,135), lev=slice(0,136))
        ls2d_z = ls2d.interp(z=z)

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
        ini = mht.Read_namelist('cabauw_rt.ini.base')

        ini['grid']['itot'] = itot
        ini['grid']['jtot'] = jtot
        ini['grid']['ktot'] = ktot

        ini['grid']['xsize'] = xsize
        ini['grid']['ysize'] = ysize
        ini['grid']['zsize'] = zsize

        ini['buffer']['zstart'] = zsize*3/4.

        if use_rt:
            ini['radiation']['swradiation'] = 'rrtmgp_rt'
            ini['radiation']['rays_per_pixel'] = spppsqp
            ini['radiation']['kngrid_i'] = 64
            ini['radiation']['kngrid_j'] = 64
            ini['radiation']['kngrid_k'] = 32
        else:
            ini['radiation']['swradiation'] = 'rrtmgp'

        ini['radiation']['swhomogenizesfc_sw'] = homogenize_sw
        ini['radiation']['swhomogenizesfc_lw'] = homogenize_lw

        ini['radiation']['swhomogenizehr_sw'] = homogenize_hr_sw
        ini['radiation']['swhomogenizehr_lw'] = homogenize_hr_lw

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
        nc_soil = nc.createGroup('soil')
        nc_soil.createDimension('z', ls2d_z.dims['zs'])
        add_nc_var('z', ('z'), nc_soil, ls2d.zs[::-1])

        add_nc_var('theta_soil', ('z'), nc_soil, theta_soil)
        add_nc_var('t_soil', ('z'), nc_soil, t_soil)
        add_nc_var('index_soil', ('z'), nc_soil, index_soil)
        add_nc_var('root_frac', ('z'), nc_soil, root_frac)

        nc.close()

        """
        Create runscript
        """
        with open('run_gpu.slurm', 'w') as f:

            f.write('#!/bin/bash\n')
            f.write('#SBATCH --job-name={}\n'.format(job_name.replace('_', '')))
            f.write('#SBATCH --output=mhh-%j.out\n')
            f.write('#SBATCH --error=mhh-%j.err\n')
            f.write('#SBATCH --partition=gpu\n')
            f.write('#SBATCH --gpus=1\n')
            f.write('#SBATCH --cpus-per-task=18\n')
            f.write('#SBATCH -t {}\n\n'.format(wc_time))

            f.write('module purge\n')
            f.write('module load 2021\n')
            f.write('module load CMake/3.20.1-GCCcore-10.3.0\n')
            f.write('module load foss/2021a\n')
            f.write('module load netCDF/4.8.0-gompi-2021a\n')
            f.write('module load CUDA/11.3.1\n\n')

            f.write('./microhh init cabauw\n')
            f.write('./microhh run cabauw')

        """
        Move files to work directory
        """
        to_move = [
                'cabauw.ini', 'cabauw_input.nc', 'run_gpu.slurm',
                'cloud_coefficients_lw.nc', 'cloud_coefficients_sw.nc',
                'coefficients_lw.nc', 'coefficients_sw.nc', 'van_genuchten_parameters.nc']
        to_copy = [microhh_bin]

        work_dir = '{0}/{1}_{2}'.format(scratch_path, itot, job_name)
        os.mkdir(work_dir)

        for f in to_move:
            shutil.move(f, work_dir)
        for f in to_copy:
            shutil.copy(f, work_dir)

    # Generate all cases.
    # Reference `rt` and `2s` cases:
    generate_case(use_rt=True)
    generate_case(use_rt=False)

    # Sensitivity "samples per pixel per spectral quadrature point".
    # Cases `rt-s032`, `rt-s064` and `rt-s128`:
    generate_case(use_rt=True, spppsqp=32)
    generate_case(use_rt=True, spppsqp=64)
    generate_case(use_rt=True, spppsqp=128)

    # Homogenized surface solar radiation.
    # Cases `rt-hom` and `2s-hom`:
    generate_case(use_rt=True, homogenize_sw=True)
    generate_case(use_rt=False, homogenize_sw=True)

    # Homogenized heating rates - sw only
    # Cases `rt-hom-hr` and `2s-hom-hr`:
    generate_case(use_rt=True, homogenize_hr_sw=True)
    generate_case(use_rt=False, homogenize_hr_sw=True)


    # Cases not used for the paper:
    # Surface shortwave and longwave radiation homogenisation;
    generate_case(use_rt=True, homogenize_sw=True, homogenize_lw=True)
    generate_case(use_rt=False, homogenize_sw=True, homogenize_lw=True)

    # Homogenized heating rates of both longwave and shortwave:
    generate_case(use_rt=True, homogenize_hr_sw=True, homogenize_hr_lw=True)
    generate_case(use_rt=False, homogenize_hr_sw=True, homogenize_hr_lw=True)
