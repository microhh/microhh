import sys, glob, os, shutil
import numpy as np

from cabauw_input import create_case_input
from microhh_tools import execute

def test_cabauw_case(executables):
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

    use_homogeneous_z0 = True
    use_homogeneous_ls = True
    use_rt = False

    # Test all executable for all Cabauw case options.
    for executable in executables:
        print('---- Testing {} ----'.format(executable))
        for use_htessel in [True, False]:
            for use_rrtmgp in [True, False]:
                for homogeneous in [True, False]:
                    for gpt_set in ['128_112', '256_224']:

                        print('Testing use_htessel={}, use_rrtmgp={}, gpt_set={}, homogeneous={}'.format(
                                use_htessel, use_rrtmgp, gpt_set, homogeneous))

                        # Cleanup...
                        to_rm = glob.glob('*000*')
                        to_rm += ['cabauw_input.nc', 'stderr.log', 'stdout.log']
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
                                homogeneous,    # z0
                                homogeneous,    # land-surface
                                gpt_set,
                                itot, jtot, ktot,
                                xsize, ysize, zsize,
                                endtime, TF)

                        # Init and run case.
                        err = 0
                        err += execute('{} init cabauw'.format(executable))
                        err += execute('{} run cabauw'.format(executable))

    
if __name__ == '__main__':
    """
    Test Cabauw case with all its options.
    """

    modes = ['cpu', 'cpumpi', 'gpu']
    precs = ['dp', 'sp']

    executables = []
    for mode in modes:
        for prec in precs:
            executables.append('../../build_{}_{}/microhh'.format(prec, mode))

    test_cabauw_case(executables)
