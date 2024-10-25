from datetime import datetime
import sys
import glob
import os
import shutil
import subprocess
import filecmp
import itertools

import numpy as np

from cabauw_input import create_case_input
import microhh_tools as mht

def clean_dir(path):
    files = glob.glob(f'{path}/*')
    for f in files:
        os.remove(f)


def execute(cmd, silent=True):
    if silent:
        return subprocess.call(cmd, shell=True, executable='/bin/bash', stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    else:
        return subprocess.call(cmd, shell=True, executable='/bin/bash')


def test_restart(
        architecture,
        precision,
        use_htessel = False,
        use_rrtmgp = False,
        use_homogeneous_z0 = True,
        use_homogeneous_ls = True,
        use_rt = False,
        use_aerosols = False,
        use_tdep_aerosols = False,
        use_tdep_gasses = False,
        use_tdep_background = False,
        sw_micro = '0'):

    executable = f'../../build_{precision}_{architecture}/microhh'
    name = f'build_{precision}_{architecture}'
    TF = np.float32 if precision=='sp' else np.float64

    print(f'Testing {name:>16s}, htessel={use_htessel:1}, rrtmgp={use_rrtmgp:1}, rt={use_rt:1}, homog_z0={use_homogeneous_z0:1} homog_ls={use_homogeneous_ls:1}, aerosols={use_aerosols:1}, tdep_aer={use_tdep_aerosols:1}, tdep_gas={use_tdep_gasses:1}, tdep_bg={use_tdep_background:1}', end='')

    itot = 4
    jtot = 4
    ktot = 128

    xsize = itot*50
    ysize = jtot*50
    zsize = 4000

    # NOTE: end_date is ignored for short restart test.
    start_date = datetime(year=2016, month=8, day=15, hour=12)
    end_date   = datetime(year=2016, month=8, day=15, hour=18)

    # Restart frequency. Case will run for 2 x restart_freq.
    restart_freq = 60

    gpt_set = '128_112'

    """
    Create (long..) list of restart files.
    """
    restart_files = [
            'thl', 'qt', 'u', 'v', 'w',
            'dudz_mo', 'dvdz_mo', 'dbdz_mo',
            'thermo_basestate', 'time']

    if sw_micro == '2mom_warm':
        restart_files += ['qr', 'nr']
    elif sw_micro == 'nsw6':
        restart_files += ['qr', 'qs', 'qg']

    if use_htessel:
        restart_files += ['wl_skin', 't_soil', 'theta_soil']

    for thermo_var in ['thl', 'qt']:
        restart_files.append(f'{thermo_var}_bot')

        if use_htessel:
            for tile in ['soil', 'veg', 'wet']:
                restart_files.append(f'{thermo_var}_bot_{tile}')

    if not use_htessel:
        restart_files += ['thl_gradbot', 'qt_gradbot']

        if sw_micro == '2mom_warm':
            restart_files.append('qr_gradbot', 'nr_gradbot')
        elif sw_micro == 'nsw6':
            restart_files.append('qr_gradbot', 'qs_gradbot', 'qg_gradbot')

    if not use_homogeneous_z0:
        if use_htessel:
            for tile in ['soil', 'veg', 'wet']:
                restart_files.append(f'obuk_{tile}')
        else:
            restart_files.append('obuk')


    """
    Cleanup old files and create new input.
    """
    to_rm = glob.glob('*00*')
    to_rm += glob.glob('*.txt')
    to_rm += glob.glob('*.bin')
    to_rm += ['cabauw.ini', 'cabauw_input.nc', 'stderr.log', 'stdout.log']

    for f in to_rm:
        if os.path.exists(f):
            os.remove(f)

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
            sw_micro,
            itot, jtot, ktot,
            xsize, ysize, zsize,
            TF)


    """
    Modify case for restart test.
    """
    ini = mht.Read_namelist('cabauw.ini')

    ini['time']['savetime'] = restart_freq
    ini['time']['endtime'] = restart_freq*2
    ini['cross']['swcross'] = False

    ini.save('cabauw.ini', allow_overwrite=True)


    """
    Run cold start.
    """
    ret = 0
    ret += execute(f'{executable} init cabauw', silent)
    ret += execute(f'{executable} run cabauw', silent)

    if ret > 0:
        print('\n')
        raise Exception('ERROR: cold start failed!')

    # Backup restart files from end of run.
    if not os.path.exists('restart_files'):
        os.makedirs('restart_files')
    else:
        clean_dir('restart_files')

    time_str = f'{2*restart_freq:07d}'
    for var in restart_files:
        shutil.move(f'{var}.{time_str}', f'restart_files/{var}.{time_str}')

    # DEBUG: binary dumps
    files = glob.glob('*.bin')
    for f in files:
        shutil.move(f, 'restart_files')

    """
    Run warm start.
    """
    ini = mht.Read_namelist('cabauw.ini')

    ini['time']['starttime'] = restart_freq
    ini['cross']['swcross'] = False

    ini.save('cabauw.ini', allow_overwrite=True)

    ret = 0
    ret += execute(f'{executable} run cabauw', silent)

    if ret > 0:
        print('\n')
        raise Exception('ERROR: warm start failed!')


    """
    Compare restart files.
    """
    err = 0
    for var in restart_files:

        f1 = f'{var}.{time_str}'
        f2 = f'restart_files/{var}.{time_str}'

        the_same = filecmp.cmp(f1, f2)

        if not the_same:
            err += 1

    if err == 0:
        print(f' -> OKAY!')
    else:
        print(f' -> ERROR: {err}/{len(restart_files)} files not identical!')


def test_permutations(
        archictecture_opts,
        precision_opts,
        htessel_opts = [False],
        rrtmgp_opts = [False],
        homogeneous_z0_opts = [True],
        homogeneous_ls_opts = [True],
        rt_opts = [False],
        aerosols_opts = [False],
        tdep_aerosols_opts = [False],
        tdep_gasses_opts = [False],
        tdep_background_opts = [False]):

    # Yes, this is not dangerous at all...
    arg_dict = locals()
    opts = list(arg_dict.values())
    permutations = itertools.product(*opts)

    for permutation in permutations:
        test_restart(
            permutation[0],
            permutation[1],
            permutation[2],
            permutation[3],
            permutation[4],
            permutation[5],
            permutation[6],
            permutation[7],
            permutation[8],
            permutation[9],
            permutation[10],
            sw_micro='0')


def test_base():
    print('--- Testing base case with radiation and land-surface ---')

    test_permutations(
            archictecture_opts = ['gpu', 'gpu_kl', 'cpu', 'cpumpi'],
            precision_opts = ['sp', 'dp'],
            htessel_opts = [True, False],
            rrtmgp_opts = [True, False])


def test_heterogeneous_surface():
    print('--- Testing heterogeneous land-surface options ---')

    test_permutations(
            archictecture_opts = ['gpu', 'cpu'],
            precision_opts = ['sp'],
            htessel_opts = [True, False],
            rrtmgp_opts = [False],
            homogeneous_z0_opts = [True, False],
            homogeneous_ls_opts = [True, False])


def test_aerosols_and_bg():
    print('--- Testing aerosols and time varying background ---')

    test_permutations(
            archictecture_opts = ['gpu', 'cpu'],
            precision_opts = ['sp'],
            htessel_opts = [False],
            rrtmgp_opts = [True],
            aerosols_opts = [True, False],
            tdep_aerosols_opts = [True, False])

    test_permutations(
            archictecture_opts = ['gpu', 'cpu'],
            precision_opts = ['sp'],
            htessel_opts = [False],
            rrtmgp_opts = [True],
            tdep_gasses_opts = [True, False],
            tdep_background_opts = [True, False])


def test_rt():
    print('--- Testing ray tracer ---')

    test_permutations(
            archictecture_opts = ['gpu'],
            precision_opts = ['sp'],
            htessel_opts = [True, False],
            rrtmgp_opts = [True],
            rt_opts = [True],
            aerosols_opts = [True, False])


if __name__ == '__main__':

    silent = True   # Suppress MicroHH stdout/err.

    test_base()
    test_heterogeneous_surface()
    test_aerosols_and_bg()
    test_rt()

    #test_restart('cpu', 'sp', use_htessel=False, use_rrtmgp=False, use_homogeneous_z0=False)
    #test_restart('cpu', 'sp', use_htessel=True, use_rrtmgp=False, use_homogeneous_z0=False)

    #test_restart('gpu', 'sp', use_htessel=False, use_rrtmgp=True, use_rt=True)
    #test_restart('gpu', 'sp', use_htessel=True, use_rrtmgp=True, use_rt=False)
    #test_restart('gpu', 'sp', use_htessel=True, use_rrtmgp=True, use_rt=True)
