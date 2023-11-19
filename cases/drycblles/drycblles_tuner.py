import sys
import os

sys.path.append('../../python/')
import microhh_tools as mht

base_opts = {
    'time': {'endtime': 1},
    'advec': {'swadvec': '2i5'},
    'force': 
    {
        'swlspres': 'geo',
        'fc': 1e-4,
        'swls': 1,
        'lslist': 'th',
        'swwls': 'mean',
        'swnudge': 1,
        'nudgelist': 'th'
    }
}

mpi_opts = {}

perm_opts = {
        'default': {},
        'diff_tke2': {'diff': {'swdiff': 'tke2'}},
        'timeloop_rk4': {'time': {'rkorder': 4}},
        'wls_local': {'force': {'swwls': 'local'}},
        'wls_mom': {'force': {'swwls': 'local', 'swwls_mom': 1}},
        'flux_lim': {'advec': {'fluxlimit_list': 'th'}}
        }

res_opts = {
    'res_32x32x32': {'grid': {'itot': 32, 'jtot': 32, 'ktot': 32}},
    'res_64x64x64': {'grid': {'itot': 64, 'jtot': 64, 'ktot': 64}},

    'res_128x128x128': {'grid': {'itot': 128, 'jtot': 128, 'ktot': 128}},
    'res_256x256x128': {'grid': {'itot': 256, 'jtot': 256, 'ktot': 128}},
    'res_512x512x128': {'grid': {'itot': 512, 'jtot': 512, 'ktot': 128}},
    'res_1024x1024x128': {'grid': {'itot': 1024, 'jtot': 1024, 'ktot': 128}},
    'res_2048x2048x128': {'grid': {'itot': 2048, 'jtot': 2048, 'ktot': 128}},
    'res_256x256x256': {'grid': {'itot': 256, 'jtot': 256, 'ktot': 256}},
    'res_512x512x256': {'grid': {'itot': 512, 'jtot': 512, 'ktot': 256}},
    'res_1024x1024x256': {'grid': {'itot': 1024, 'jtot': 1024, 'ktot': 256}}
}


if __name__ == '__main__':

    case_name = 'drycblles'
    mode='gpu'
    case_dir='.'


    executables = ['microhh_sp_gpu', 'microhh_dp_gpu']
    # executables = ['microhh']

    os.environ['KERNEL_LAUNCHER_TUNE'] = '*'

    for executable in executables:
        experiment='tuner_'+executable

        mht.run_permutations(
            case_name, base_opts, mpi_opts, [perm_opts, res_opts],
            executable=executable, mode=mode, case_dir=case_dir, experiment=experiment)
