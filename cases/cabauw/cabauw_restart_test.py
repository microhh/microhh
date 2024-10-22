from datetime import datetime
import sys
import glob
import os
import shutil
import subprocess
import filecmp

import numpy as np

from cabauw_input import create_case_input
import microhh_tools as mht

def clean_dir(path):
    files = glob.glob(f'{path}/*')
    for f in files:
        os.remove(f)


def execute(cmd):
    return subprocess.call(cmd, shell=True, executable='/bin/bash')


def test_restart(
        architecture,
        precision,
        use_kernel_launcher = False,
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

    execute_func = mht.execute if silent else execute

    itot = 8
    jtot = 8
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
    
    kl_flag = '_kl' if use_kernel_launcher else ''
    executable = f'../../build_{precision}_{architecture}{kl_flag}/microhh'

    TF = np.float32 if precision=='sp' else np.float64
    
    print(f'Testing {executable:>30s}, use_htessel={use_htessel:1}, use_rrtmgp={use_rrtmgp:1}', end='')
    
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
        restart_files.append('obuk')
    
    
    """
    Cleanup old files and create new input.
    """
    to_rm = glob.glob('*00*')
    to_rm += glob.glob('*.txt')
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
    ret += execute_func(f'{executable} init cabauw')
    ret += execute_func(f'{executable} run cabauw')

    
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
    
    
    """
    Run warm start.
    """
    ini = mht.Read_namelist('cabauw.ini')
    
    ini['time']['starttime'] = restart_freq
    ini['cross']['swcross'] = False

    ini.save('cabauw.ini', allow_overwrite=True)
    
    ret = 0
    ret += execute_func(f'{executable} run cabauw')
    
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
            #arr1 = np.fromfile(f1, dtype=TF)
            #arr2 = np.fromfile(f2, dtype=TF)
    
            #max_diff = np.abs(arr1-arr2).max()
    
            #print(f'Variable {var} not identical, max abs diff={max_diff}')
    
            err += 1
    
    if err == 0:
        print(f' -> OKAY!')
    else:
        print(f' -> ERROR: {err}/{len(restart_files)} files not identical!')




if __name__ == '__main__':

    silent = True

    #test_restart('gpu', 'sp', use_htessel=True, use_rrtmgp=True)

    prec = 'sp'
    for use_rrtmgp in (False, True):
        for use_htessel in (False, True):
            for arch in ('gpu', 'cpu', 'cpumpi'):
                test_restart(arch, prec, use_htessel=use_htessel, use_rrtmgp=use_rrtmgp)
