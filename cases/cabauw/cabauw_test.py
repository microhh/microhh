import sys, glob, os, shutil
import subprocess
import numpy as np

from cabauw_input import create_case_input
from microhh_tools import execute

#def execute(task):
#    return subprocess.call(task, shell=True, executable='/bin/bash')


def test_cabauw_case(executable):
    """
    Test Cabauw case with all options.
    """

    itot = 16
    jtot = 16
    ktot = 128

    xsize = itot*50
    ysize = jtot*50
    zsize = 4000

    endtime = 10

    # Don't test ray-tracer for now..
    use_rt = False
    gpt_set = '128_112'

    for use_htessel in [True, False]:
        for use_rrtmgp in [True, False]:

            # Only test heterogeneous mode with LSM:
            if use_htessel:
                homogeneous_ls_opts = [True, False]
                homogeneous_z0_opts = [True, False]
            else:
                homogeneous_ls_opts = [True]
                homogeneous_z0_opts = [True]

            for use_homogeneous_z0 in homogeneous_z0_opts:
                for use_homogeneous_ls in homogeneous_ls_opts:

                    print('Testing use_htessel={}, use_rrtmgp={}, gpt_set={}, homog_z0={}, homog_ls={}'.format(
                            use_htessel, use_rrtmgp, gpt_set, use_homogeneous_z0, use_homogeneous_ls))

                    # Cleanup...
                    to_rm = glob.glob('*000*')
                    to_rm += ['cabauw.ini', 'cabauw_input.nc', 'stderr.log', 'stdout.log']
                    for f in to_rm:
                        if os.path.exists(f):
                            os.remove(f)

                    # Yikes...
                    TF = np.float32 if 'sp' in executable else np.float64

                    # Create input files and .ini file from .ini.base.
                    create_case_input(
                            use_htessel,
                            use_rrtmgp,
                            use_rt,
                            use_homogeneous_z0,
                            use_homogeneous_ls,
                            gpt_set,
                            itot, jtot, ktot,
                            xsize, ysize, zsize,
                            endtime, TF)

                    # Init and run case.
                    ret = 0
                    ret += execute('{} init cabauw'.format(executable))
                    ret += execute('{} run cabauw'.format(executable))

                    if ret != 0:
                        raise Exception('Case failed for {}'.format(executable))

    
if __name__ == '__main__':
    """
    Test Cabauw case with all its options.
    """

    modes = ['gpu', 'cpu', 'cpumpi']
    precs = ['dp', 'sp']

    executables = []
    for mode in modes:
        for prec in precs:
            executable = '../../build_{}_{}/microhh'.format(prec, mode)
            print('--- Testing {} ----'.format(executable))
            test_cabauw_case(executable)
