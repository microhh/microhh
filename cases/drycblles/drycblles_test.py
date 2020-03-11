from copy import deepcopy
import sys
import os
import copy
import shutil

sys.path.append('../../python/')
import microhh_tools as mht

# Case configuration dicts
no_opts = {}

opt_mpi = {
        'master': {'npx': 2, 'npy': 2}}

opt_small = {
        'grid': {'itot': 8, 'jtot': 8},
        'time': {'endtime': 3600, 'savetime': 1800}}

opt_nostats = {
        'stats': {'swstats': 0}}


if __name__ == '__main__':

    case_name = 'drycblles'

    kwargs = dict([arg.split('=') for arg in sys.argv[2:]])

    if len(sys.argv) > 1:
        function_name = sys.argv[1]

        if function_name == 'run_case':
            mht.run_case(case_name, no_opts, opt_mpi, **kwargs)
        elif function_name == 'run_small':
            mht.run_case(case_name, opt_small, opt_mpi, **kwargs)
        elif function_name == 'run_restart':
            mht.run_restart(case_name, opt_small, opt_mpi, None, **kwargs)
        else:
            raise Exception('\"{}\" is an invalid option'.format(function_name))

    else:
        mht.run_case(case_name, no_opts, opt_mpi)
