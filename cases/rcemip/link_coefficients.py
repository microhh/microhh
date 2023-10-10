import os

if os.path.exists('../../rte-rrtmgp-cpp/rrtmgp-data/extensions/rrtmgp-clouds-lw.nc'):
    # Jeej, the simple version...
    microhh_root = '../..'

else:
    # Oiii, the less simple version. 
    abspath = os.path.abspath('.')
    dirs = abspath.split('/')[1:]
    if 'microhh' in dirs:
        microhh_root = '/'+'/'.join(dirs[:dirs.index('microhh')+1])
    else:
        raise Exception('Cannot determine MicroHH root...')

os.symlink('{}/rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-gas-lw-g128.nc'.format(microhh_root), 'coefficients_lw.nc')
os.symlink('{}/rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-gas-sw-g112.nc'.format(microhh_root), 'coefficients_sw.nc')

os.symlink('{}/rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-clouds-lw.nc'.format(microhh_root), 'cloud_coefficients_lw.nc')
os.symlink('{}/rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-clouds-sw.nc'.format(microhh_root), 'cloud_coefficients_sw.nc')
