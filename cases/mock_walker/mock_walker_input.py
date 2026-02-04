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

import os
import numpy as np

import ls2d

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

"""
# ECMWF HPC
gpt_path = '/home/nkbs/meteo/models/coefficients_veerman'
microhh_path = '/home/nkbs/meteo/models/microhh'
microhh_bin = '/home/nkbs/meteo/models/microhh/build_sp_dpfft_cpumpi/microhh'
work_dir = '/scratch/nkbs/mock_walker_xl_400m'
"""

# LUMI
project = 'project_465002576'
partition = 'small'
gpt_path = '/users/stratumv/meteo/models/coefficients_veerman'
microhh_path = '/users/stratumv/meteo/models/microhh'
microhh_bin = '/users/stratumv/meteo/models/microhh/build_spdp_cpumpi/microhh'
work_dir = f'/scratch/{project}/mock_walker_io'

sw_cos_sst = True
rotated_domain = False
mean_sst = 300
d_sst = 2.5
ps = 101480


# Small test domain.
xsize = 512*400
ysize = 512*400

itot = 512
jtot = 512

npx = 8
npy = 16


"""
# Full domain, 400 m resolution.
# 128 x 32 :  15360 x  1024 x 128 @ 128 x  32, #/core = 491520.0
xsize = 15360*400
ysize = 1024*400

itot = 15360
jtot = 1024

npx = 128
npy = 32
"""

"""
# Full domain, 200 m resolution.
# 128 x 32 :  30720 x  2048 x 128 @ 128 x  32, #/core = 1966080.0
xsize = 30720*200
ysize = 2048*200

itot = 30720
jtot = 2048

npx = 128
npy = 32
"""

#endtime = 10*24*3600
endtime = 3600

"""
# Official RCEMIP LES grid
z = np.loadtxt('rcemip_les_grid.txt')
z = z[:-2]   # From 146 to 144 levels for domain decomposition.
zsize = 32250

zh = 0.5*(z[:-1] + z[1:])
zh = np.append(0., zh)
zh = np.append(zh, zsize)
"""

z = np.array([0, 2_000, 20_000, 100_000])
f = np.array([1.05, 1.012, 1.04])
grid = ls2d.grid.Grid_stretched_manual(128, 40, z, f)
z = grid.z
zsize = grid.zsize

# G-point sets from Veerman (2024).
# For now the cheapest for testing. TBD with Martin/Menno.
coef_sw = 'rrtmgp-gas-sw-g049-cf2.nc'
coef_lw = 'rrtmgp-gas-lw-g056-cf2.nc'

name = 'mock_walker'
create_slurm_script = True
wc_time = '01:00:00'

if not os.path.exists(work_dir):
    os.makedirs(work_dir)

mock_walker_input(
        name,
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
        work_dir,
        gpt_path,
        microhh_path,
        microhh_bin,
        create_slurm_script,
        project,
        partition,
        float_type)
