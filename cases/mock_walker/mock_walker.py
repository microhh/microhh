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

import shutil

import numpy as np

# `pip install microhhpy`
from microhhpy.io import read_ini, check_ini, save_ini, save_case_input
from microhhpy.thermo import qsat, exner


"""
Constants.
"""
Rd = 287.04
Rv = 461.5
g = 9.79764
Rd = 287.04
cp = 1005.
eps = 18.01528 / 28.9647   # molar mass water / molar mass air
p0 = 1e5


"""
Help functions.
"""
def calc_profiles(z, ps, mean_sst, mean_q0):
    """
    Create vertical profiles p, q, T, thl, o3 following RCEMIP-I specification.
    """

    z_q1 = 4.0e3
    z_q2 = 7.5e3
    z_t = 15.e3
    q_t = 1.e-14

    q = mean_q0 * np.exp(-z / z_q1) * np.exp(-(z / z_q2)**2)

    # CvH hack to remove moisture jump.
    q_tb = mean_q0 * np.exp(-z_t/z_q1) * np.exp(-(z_t/z_q2)**2)
    q -= q_tb + q_t

    i_above_zt = np.where(z >= z_t)
    q[i_above_zt] = q_t

    gamma = 6.7e-3
    Tv_0 = (1. + 0.608 * mean_q0) * mean_sst
    Tv = Tv_0 - gamma*z
    Tv_t = Tv_0 - gamma*z_t
    Tv[i_above_zt] = Tv_t
    T = Tv / (1. + 0.608*q)

    p = ps * (Tv / Tv_0)**(g/(Rd*gamma))
    p_tmp = ps * (Tv_t/Tv_0)**(g/(Rd*gamma)) \
          * np.exp( -( (g*(z-z_t)) / (Rd*Tv_t) ) )
    p[i_above_zt] = p_tmp[i_above_zt]

    thl = T * (p0/p)**(Rd/cp)

    g1 = 3.6478
    g2 = 0.83209
    g3 = 11.3515
    p_hpa = p/100.
    o3 = g1 * p_hpa**g2 * np.exp(-p_hpa/g3) * 1e-6

    return p, q, T, thl, o3


def mock_walker_input(
        name,
        xsize,
        ysize,
        itot,
        jtot,
        npx,
        npy,
        z,
        zsize,
        endtime,
        sw_cos_sst,
        mean_sst,
        d_sst,
        ps,
        rotated_domain,
        coef_sw,
        coef_lw,
        wc_time,
        work_dir,
        gpt_path,
        microhh_path,
        microhh_bin,
        create_slurm_script=False,
        account=None,
        float_type=np.float32):
    """
    Create input files for Mock Walker case.
    """

    # Table 2 W17.
    if mean_sst == 295:
        mean_q0 = 12.0e-3
    elif mean_sst == 300:
        mean_q0 = 18.65e-3
    elif mean_sst == 305:
        mean_q0 = 24.00e-3
    else:
        raise Exception('Unknown mean_sst')


    """
    Read ini file and set values.
    """
    ini = read_ini('mock_walker.ini.base')

    ini['master']['npx'] = npx
    ini['master']['npy'] = npy

    ini['grid']['itot'] = itot
    ini['grid']['jtot'] = jtot
    ini['grid']['ktot'] = z.size

    ini['grid']['xsize'] = xsize
    ini['grid']['ysize'] = ysize
    ini['grid']['zsize'] = zsize

    if sw_cos_sst:
        ini['boundary']['sbot_2d_list'] = ['thl', 'qt']

    ini['buffer']['zstart'] = 0.75*zsize
    ini['time']['endtime'] = endtime

    ini['thermo']['pbot'] = ps

    # For `!sw_cos_sst`:
    exn = exner(ps)
    ini['boundary']['sbot[thl]'] = mean_sst / exn
    ini['boundary']['sbot[qt]'] = qsat(ps, mean_sst)

    # Check and write to final .ini file
    check_ini(ini)
    save_ini(ini, f'{work_dir}/mock_walker.ini')


    """
    Vetical profiles on radiation grid.
    """
    z_top = 70.e3
    dz = 500.
    z_lay  = np.arange(dz/2, z_top, dz)
    z_lev = np.arange(0, z_top-dz/2, dz)
    z_lev = np.append(z_lev, z_top)

    p_lay, q_lay, T_lay, _, o3_lay = calc_profiles(z_lay, ps, mean_sst, mean_q0)
    p_lev, _,     T_lev, _, _      = calc_profiles(z_lev, ps, mean_sst, mean_q0)

    h2o_lay = q_lay / (eps - eps * q_lay)

    background_concs = {
        'co2' :  348.e-6,
        'ch4' : 1650.e-9,
        'n2o' :  306.e-9,
        'n2' : 0.7808,
        'o2' : 0.2095,
        }


    """
    Vertical profiles on LES grid.
    """
    _, qt, _, thl, o3 = calc_profiles(z, ps, mean_sst, mean_q0)
    h2o = qt / (eps - eps * qt)


    """
    Create case_input.nc.
    """
    init = {
        'z' : z,
        'thl' : thl,
        'qt' : qt,
        'h2o' : h2o,
        'o3' : o3,
        }

    radiation  = {
        'z_lay': z_lay,
        'z_lev': z_lev,
        'p_lay': p_lay,
        'p_lev': p_lev,
        't_lay': T_lay,
        't_lev': T_lev,
        'o3':    o3_lay,
        'h2o':   h2o_lay
        }

    init.update(background_concs)
    radiation.update(background_concs)

    save_case_input('mock_walker', init_profiles=init, radiation=radiation, output_dir=work_dir)


    if (sw_cos_sst):
        """
        2D surface fields with `SST(x) = <SST> - d_SST / 2 * cos(2 pi x / xsize)`.
        if rotated_domain, the gradient is put in the y-direction.
        """

        thl_sbot = np.zeros((jtot, itot), dtype=float_type)
        qt_sbot  = np.zeros((jtot, itot), dtype=float_type)

        if rotated_domain:
            dy = ysize / jtot
            y = np.arange(dy/2, ysize, dy)
            sst_y = mean_sst - d_sst / 2 * np.cos(2 * np.pi * y / ysize)

            exn = exner(ps)
            thl_y = sst_y / exn
            qt_y  = qsat(ps, sst_y)

            thl_sbot[:,:] = thl_y[:, None]
            qt_sbot[:,:] = qt_y[:, None]

        else:
            dx = xsize / itot
            x = np.arange(dx/2, xsize, dx)
            sst_x = mean_sst - d_sst / 2 * np.cos(2 * np.pi * x / xsize)

            exn = exner(ps)
            thl_x = sst_x / exn
            qt_x  = qsat(ps, sst_x)

            thl_sbot[:,:] = thl_x[None, :]
            qt_sbot[:,:] = qt_x[None, :]


        thl_sbot.tofile(f'{work_dir}/thl_bot_in.0000000')
        qt_sbot.tofile(f'{work_dir}/qt_bot_in.0000000')


    """
    Copy radiation files, executables, et cetera.
    """
    rrtmgp_path = f'{microhh_path}/rte-rrtmgp-cpp/rrtmgp-data'

    shutil.copy2(f'{gpt_path}/{coef_sw}', f'{work_dir}/coefficients_sw.nc')
    shutil.copy2(f'{gpt_path}/{coef_lw}', f'{work_dir}/coefficients_lw.nc')

    shutil.copy2(f'{rrtmgp_path}/rrtmgp-clouds-lw.nc', f'{work_dir}/cloud_coefficients_lw.nc')
    shutil.copy2(f'{rrtmgp_path}/rrtmgp-clouds-sw.nc', f'{work_dir}/cloud_coefficients_sw.nc')

    shutil.copy2(microhh_bin, f'{work_dir}/microhh')


    """
    Create SLURM script.
    """
    if create_slurm_script:

        slurm_script = f'{work_dir}/run.slurm'
        with open(slurm_script, 'w') as f:

            f.write(f'#!/bin/bash\n')
            if account is not None:
                f.write(f'#SBATCH --account={account}\n')
            f.write(f'#SBATCH --job-name={name}\n')
            f.write(f'#SBATCH --output={work_dir}/mhh-%j.out\n')
            f.write(f'#SBATCH --error={work_dir}/mhh-%j.err\n')
            f.write(f'#SBATCH --partition=par\n')
            f.write(f'#SBATCH --ntasks={npx*npy}\n')
            f.write(f'#SBATCH --cpus-per-task=1\n')
            f.write(f'#SBATCH --ntasks-per-core=1\n')
            f.write(f'#SBATCH --mem=224G\n')
            f.write(f'#SBATCH --time={wc_time}\n\n')

            f.write(f'source ~/setup_env.sh\n\n')

            f.write(f'cd {work_dir}\n')
            f.write(f'srun ./microhh init mock_walker\n')
            f.write(f'srun ./microhh run mock_walker')

        return slurm_script

    else:
        return None
