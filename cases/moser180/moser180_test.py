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
        'grid': {'itot': 64, 'jtot': 48, 'ktot': 32},
        'time': {'endtime': 10}}

opt_small_restart = {
        'grid': {'itot': 32, 'jtot': 16, 'ktot': 32},
        'time': {'endtime': 200, 'savetime': 100}}

# Case configuration dicts with name label for permutations.
dict_opts = {
        'default': {},
        'no_advec': {'advec': {'swadvec': 0}},
        'advec_4m': {'advec': {'swadvec': '4m'}},
        'no_diff': {'diff': {'swdiff': 0}}}


if __name__ == '__main__':

    case_name = 'moser180'

    kwargs = dict([arg.split('=') for arg in sys.argv[2:]])

    if len(sys.argv) > 1:
        function_name = sys.argv[1]

        if function_name == 'run_case':
            mht.run_case(case_name, no_opts, opt_mpi, **kwargs)
        elif function_name == 'run_small':
            mht.run_case(case_name, opt_small, opt_mpi, **kwargs)
        elif function_name == 'run_restart':
            mht.run_restart(case_name, opt_small, opt_mpi, dict_opts, **kwargs)
        else:
            raise Exception('\"{}\" is an invalid option'.format(function_name))

    else:
        mht.run_case(case_name, no_opts, opt_mpi)
