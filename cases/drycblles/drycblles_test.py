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

# Case configuration dicts with name label for permutations.
dict_opts = {
        'default': {},
        'advec_2i5': {'advec': {'swadvec': '2i5'}},
        'deardorff': {
            'diff': {'swdiff': 'tke2'},
            'boundary': {'sbot[sgstke]': 0, 'stop[sgstke]': 0},
            'advec': {'swadvec': '2i5', 'fluxlimit_list': 'sgstke'}}}

dict_res = {
    'res_32x32x32': {'grid': {'itot': 32, 'jtot': 32, 'ktot': 32}},
    'res_64x64x64': {'grid': {'itot': 64, 'jtot': 64, 'ktot': 64}},
    'res_128x128x128': {'grid': {'itot': 128, 'jtot': 128, 'ktot': 128}},
    'res_256x256x256': {'grid': {'itot': 256, 'jtot': 256, 'ktot': 256}},
    'res_512x512x512': {'grid': {'itot': 512, 'jtot': 512, 'ktot': 512}},
    'res_512x512x256': {'grid': {'itot': 512, 'jtot': 512, 'ktot': 256}},
    'res_1024x1024x256': {'grid': {'itot': 1024, 'jtot': 1024, 'ktot': 256}},
    'res_1024x1024x128': {'grid': {'itot': 1024, 'jtot': 1024, 'ktot': 128}},
    'res_2048x2048x128': {'grid': {'itot': 2048, 'jtot': 2048, 'ktot': 128}}

}


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
            mht.run_restart(case_name, opt_small, opt_mpi, dict_opts, **kwargs)
        else:
            raise Exception('\"{}\" is an invalid option'.format(function_name))

    else:
        mode='gpu'
        case_dir='.'
        experiment='local'
        executable='microhh'
        
        mht.run_permutations(
            case_name, no_opts, opt_mpi, [dict_res],
            executable=executable, mode=mode, case_dir=case_dir, experiment=experiment)
