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
Case settings.
"""
float_type = np.float32

"""
# Eddy
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

# ECMWF HPC
gpt_path = '/home/nkbs/meteo/models/coefficients_veerman'
microhh_path = '/home/nkbs/meteo/models/microhh'
microhh_bin = '/home/nkbs/meteo/models/microhh/build_sp_cpumpi/microhh'
work_dir = '/scratch/nkbs/mock_walker_scaling_v1/'


z = np.array([0, 2_000, 20_000, 100_000])
f = np.array([1.05, 1.012, 1.04])
grid = ls2d.grid.Grid_stretched_manual(128, 40, z, f)

sw_cos_sst = False  # False for RCEMIP-I, True for RCEMIP-II
mean_sst = 300
d_sst = 2.5
ps = 101480

endtime = 1800

itot_base = 1920
npx_base = 8
npy_base = 16
dxy = 200

coef_sw='rrtmgp-gas-sw-g049-cf2.nc'
coef_lw='rrtmgp-gas-lw-g056-cf2.nc'

create_slurm_script = True
wc_time = '04:00:00'

# Configurations weak-scaling test.
configurations = [
    ('16x1', 2048, [(1,1), (2,1), (4,1), (8,1), (16,1)]),
    ('16x2', 1024, [(1,1), (2,2), (4,2), (8,2), (16,2)]),
    ('16x4', 512,  [(1,1), (2,2), (4,4), (8,4), (16,4)]),
]

for name, jtot_node, configs in configurations:
    print(f'Configuration: {name} nodes')

    for nn_x, nn_y in configs:
        print(f'Running {nn_x} x {nn_y}')

        itot = itot_base * nn_x
        jtot = jtot_node * nn_y

        npx = npx_base * nn_x
        npy = npy_base * nn_y

        xsize = itot * dxy
        ysize = jtot * dxy

        ny = name.split('x')[-1]
        case_name = f'{ny}x{nn_x}x{nn_y}'
        case_path = f'{work_dir}/{case_name}'

        if not os.path.exists(case_path):
            os.makedirs(case_path)

        slurm_script = mock_walker_input(
                name,
                xsize,
                ysize,
                itot,
                jtot,
                npx,
                npy,
                grid.z,
                grid.zh,
                endtime,
                sw_cos_sst,
                mean_sst,
                d_sst,
                ps,
                coef_sw,
                coef_lw,
                wc_time,
                case_path,
                gpt_path,
                microhh_path,
                microhh_bin,
                create_slurm_script,
                float_type)

        subprocess.run(['sbatch', slurm_script], check=True)
