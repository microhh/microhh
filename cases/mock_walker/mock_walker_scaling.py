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

import subprocess
import os

import numpy as np

# `pip install ls2d`
import ls2d

# `mock_walker.py` help functions.
from mock_walker import mock_walker_input


"""
# Eddy
project = None
partition = None
gpt_path = '/home/bart/meteo/models/coefficients_veerman/' 
microhh_path = '/home/bart/meteo/models/microhh/'
microhh_bin = '/home/bart/meteo/models/microhh/build_sp_cpumpi/microhh'
work_dir = 'test'
"""

"""
# Snellius
gpt_path = '/gpfs/work3/0/lesmodels/team_bart/coefficients_veerman'
microhh_path = '/home/bstratum/meteo/models/microhh'
microhh_bin = '/home/bstratum/meteo/models/microhh/build_sp_cpumpi/microhh'
work_dir = '/scratch-shared/bstratum/mock_walker_test'
"""

"""
# ECMWF HPC
gpt_path = '/home/nkbs/meteo/models/coefficients_veerman'
microhh_path = '/home/nkbs/meteo/models/microhh'
microhh_bin = '/home/nkbs/meteo/models/microhh/build_sp_cpumpi/microhh'
work_dir = '/scratch/nkbs/mock_walker_scaling_v1/'
"""


# LUMI
project = 'project_465002576'
partition = 'standard'
gpt_path = '/users/stratumv/meteo/models/coefficients_veerman'
microhh_path = '/users/stratumv/meteo/models/microhh'
microhh_bin = '/users/stratumv/meteo/models/microhh/build_spdp_cpumpi/microhh'
work_dir = f'/scratch/{project}/mock_walker_scaling_v1'


"""
Global settings.
"""
float_type = np.float32

sw_cos_sst = False  # False for RCEMIP-I, True for RCEMIP-II
mean_sst = 300
d_sst = 2.5
ps = 101480
endtime = 1800

#coef_sw='rrtmgp-gas-sw-g049-cf2.nc'
#coef_lw='rrtmgp-gas-lw-g056-cf2.nc'

coef_sw = 'rrtmgp-gas-sw-g112.nc'
coef_lw = 'rrtmgp-gas-lw-g128.nc'

create_slurm_script = True
wc_time = '01:00:00'


def run_weak_scaling(settings, rotated_domain):
    """
    Run weak scaling set for all configs.
    """
    name = settings['name']
    print(f'====== Running {name} ======')

    if settings['ktot'] == 128:
        # Custom 128 layer grid.
        z = np.array([0, 2_000, 20_000, 100_000])
        f = np.array([1.05, 1.012, 1.04])
        grid = ls2d.grid.Grid_stretched_manual(128, 40, z, f)
        z = grid.z
        zsize = grid.zsize
    elif settings['ktot'] == 144:
        # Official RCEMIP LES grid.
        z = np.loadtxt('rcemip_les_grid.txt')
        z = z[:-2]   # From 146 to 144 levels for domain decomposition.
        zsize = 32250
    else:
        raise Exception('Invalid vertical grid.')

    for nn_x, nn_y in settings['configs']:

        itot = settings['itot_b'] * nn_x
        jtot = settings['jtot_b'] * nn_y

        npx = settings['npx_b'] * nn_x
        npy = settings['npy_b'] * nn_y

        xsize = itot * settings['dxy']
        ysize = jtot * settings['dxy']

        case_name = f'{name}_{nn_x}x{nn_y}'
        case_path = f'{work_dir}/{name}/{nn_x}x{nn_y}'

        if not os.path.exists(case_path):
            os.makedirs(case_path)

        copy_out_to = work_dir

        slurm_script = mock_walker_input(
                case_name,
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
                case_path,
                gpt_path,
                microhh_path,
                microhh_bin,
                create_slurm_script,
                project,
                partition,
                copy_out_to,
                float_type)

        subprocess.run(['sbatch', slurm_script], check=True)


"""
200 m, 128 levels.
"""
settings = (
    dict(name='200_128_16x1', itot_b=1920, jtot_b=2048, ktot=128, npx_b=8, npy_b=16, dxy=195.3125, configs=[(1,1), (2,1), (4,1), (8,1), (16,1)]),
    dict(name='200_128_16x2', itot_b=1920, jtot_b=1024, ktot=128, npx_b=8, npy_b=16, dxy=195.3125, configs=[(1,1), (1,2), (2,2), (4,2), (8,2), (16,2)]),
    dict(name='200_128_16x4', itot_b=1920, jtot_b= 512, ktot=128, npx_b=8, npy_b=16, dxy=195.3125, configs=[(1,1), (1,2), (1,4), (2,4), (4,4), (8,4), (16,4)]),
    dict(name='200_128_16x8', itot_b=1920, jtot_b= 256, ktot=128, npx_b=8, npy_b=16, dxy=195.3125, configs=[(1,1), (1,2), (1,4), (1,8), (2,8), (4,8), (8,8), (16,8)]),
    ) 

for setting in settings:
    run_weak_scaling(setting, False)

"""
200 m, 128 levels, rotated.
"""
settings = (
    dict(name='200_128_r_1x16',  itot_b=2048, jtot_b=1920, ktot=128, npx_b=16, npy_b=8, dxy=195.3125, configs=[(1,1), (1,2), (1,4), (1,8), (1,16)]),
    dict(name='200_128_r_1x32',  itot_b=2048, jtot_b= 960, ktot=128, npx_b=16, npy_b=8, dxy=195.3125, configs=[(1,1), (1,2), (1,4), (1,8), (1,16), (1,32)]),
    dict(name='200_128_r_1x64',  itot_b=2048, jtot_b= 480, ktot=128, npx_b=16, npy_b=8, dxy=195.3125, configs=[(1,1), (1,2), (1,4), (1,8), (1,16), (1,32), (1,64)]),
    dict(name='200_128_r_1x128', itot_b=2048, jtot_b= 240, ktot=128, npx_b=16, npy_b=8, dxy=195.3125, configs=[(1,1), (1,2), (1,4), (1,8), (1,16), (1,32), (1,64), (1,128)]),
    )

for setting in settings:
    run_weak_scaling(setting, True)

"""
400 m, 128 vertical levels.
"""
settings = (
    dict(name='400_128_16x1', itot_b=960, jtot_b=1024, ktot=128, npx_b=8, npy_b=16, dxy=390.625, configs=[(1,1), (2,1), (4,1), (8,1), (16,1)]),
    dict(name='400_128_16x2', itot_b=960, jtot_b= 512, ktot=128, npx_b=8, npy_b=16, dxy=390.625, configs=[(1,1), (1,2), (2,2), (4,2), (8,2), (16,2)]),
    dict(name='400_128_16x4', itot_b=960, jtot_b= 256, ktot=128, npx_b=8, npy_b=16, dxy=390.625, configs=[(1,1), (1,2), (1,4), (2,4), (4,4), (8,4), (16,4)])
)

for setting in settings:
    run_weak_scaling(setting, False)

"""
400 m, 128 vertical levels, rotated.
"""
settings = (
    dict(name='400_128_r_1x16', itot_b=1024, jtot_b=960, ktot=128, npx_b=16, npy_b=8, dxy=390.625, configs=[(1,1), (1,2), (1,4), (1,8), (1,16)]),
    dict(name='400_128_r_1x32', itot_b=1024, jtot_b=480, ktot=128, npx_b=16, npy_b=8, dxy=390.625, configs=[(1,1), (1,2), (1,4), (1,8), (1,16), (1,32)]),
    dict(name='400_128_r_1x64', itot_b=1024, jtot_b=240, ktot=128, npx_b=16, npy_b=8, dxy=390.625, configs=[(1,1), (1,2), (1,4), (1,8), (1,16), (1,32), (1,64)]),
    )

for setting in settings:
    run_weak_scaling(setting, True)

"""
200 m, 144 levels, rotated.
"""
settings = (
    dict(name='200_144_r_1x16',  itot_b=2048, jtot_b=1920, ktot=144, npx_b=16, npy_b=8, dxy=195.3125, configs=[(1,1), (1,2), (1,4), (1,8), (1,16)]),
    dict(name='200_144_r_1x32',  itot_b=2048, jtot_b= 960, ktot=144, npx_b=16, npy_b=8, dxy=195.3125, configs=[(1,1), (1,2), (1,4), (1,8), (1,16), (1,32)]),
    dict(name='200_144_r_1x64',  itot_b=2048, jtot_b= 480, ktot=144, npx_b=16, npy_b=8, dxy=195.3125, configs=[(1,1), (1,2), (1,4), (1,8), (1,16), (1,32), (1,64)]),
    dict(name='200_144_r_1x128', itot_b=2048, jtot_b= 240, ktot=144, npx_b=16, npy_b=8, dxy=195.3125, configs=[(1,1), (1,2), (1,4), (1,8), (1,16), (1,32), (1,64), (1,128)]),
    )

for setting in settings:
    run_weak_scaling(setting, True)


"""
400 m, 144 levels, rotated.
"""
settings = (
    dict(name='400_144_r_1x16', itot_b=1024, jtot_b=960, ktot=144, npx_b=16, npy_b=8, dxy=390.625, configs=[(1,1), (1,2), (1,4), (1,8), (1,16)]),
    dict(name='400_144_r_1x32', itot_b=1024, jtot_b=480, ktot=144, npx_b=16, npy_b=8, dxy=390.625, configs=[(1,1), (1,2), (1,4), (1,8), (1,16), (1,32)]),
    dict(name='400_144_r_1x64', itot_b=1024, jtot_b=240, ktot=144, npx_b=16, npy_b=8, dxy=390.625, configs=[(1,1), (1,2), (1,4), (1,8), (1,16), (1,32), (1,64)]),
    )

for setting in settings:
    run_weak_scaling(setting, True)
