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

import numpy as np

from mock_walker import mock_walker_input


"""
Case settings.
"""
float_type = np.float32

# Eddy
gpt_path = '/home/bart/meteo/models/coefficients_veerman/' 
microhh_path = '/home/bart/meteo/models/microhh/'
microhh_bin = '/home/bart/meteo/models/microhh/build_sp_cpumpi/microhh'
work_dir = 'test'

"""
# Snellius
gpt_path = '/gpfs/work3/0/lesmodels/team_bart/coefficients_veerman'
microhh_path = '/home/bstratum/meteo/models/microhh'
microhh_bin = '/home/bstratum/meteo/models/microhh/build_sp_cpumpi/microhh'
work_dir = '/scratch-shared/bstratum/mock_walker_test'
"""

sw_cos_sst = True
mean_sst = 300
d_sst = 2.5
ps = 101480

xsize = 128000
ysize = 32000

itot = 256
jtot = 64

npx = 4
npy = 2

endtime = 24*3600

# Official RCEMIP LES grid
z = np.loadtxt('rcemip_les_grid.txt')
z = z[:-2]   # From 146 to 144 levels for domain decomposition.
zsize = 32250

zh = 0.5*(z[:-1] + z[1:])
zh = np.append(0., zh)
zh = np.append(zh, zsize)

# G-point sets from Veerman (2024).
# For now the cheapest for testing. TBD with Martin/Menno.
coef_sw = 'rrtmgp-gas-sw-g049-cf2.nc'
coef_lw = 'rrtmgp-gas-lw-g056-cf2.nc'

name = 'mock_walker'
create_slurm_script = True
wc_time = '04:00:00'

mock_walker_input(
        name,
        xsize,
        ysize,
        itot,
        jtot,
        npx,
        npy,
        z,
        zh,
        endtime,
        sw_cos_sst,
        mean_sst,
        d_sst,
        ps,
        coef_sw,
        coef_lw,
        wc_time,
        work_dir,
        gpt_path,
        microhh_path,
        microhh_bin,
        create_slurm_script,
        float_type)
