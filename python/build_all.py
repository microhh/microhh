#
#  MicroHH
#  Copyright (c) 2011-2023 Chiel van Heerwaarden
#  Copyright (c) 2011-2023 Thijs Heus
#  Copyright (c) 2014-2023 Bart van Stratum
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

def execute(task):
    subprocess.call(task, shell=True, executable='/bin/bash')

TFs = ['sp', 'dp']
archs = ['cpu', 'cpumpi', 'gpu']
build_type = 'release'

for arch in archs:
    for TF in TFs:
        build_dir = '../build_{}_{}'.format(TF, arch)
        if not os.path.exists(build_dir):
            os.mkdir(build_dir)
        os.chdir(build_dir)

        usecuda = (arch == 'gpu')
        usempi  = (arch == 'cpumpi')
        usesp   = (TF == 'sp')

        cmake_cmd = 'cmake -DUSEMPI={} -DUSECUDA={} -DUSESP={} -DCMAKE_BUILD_TYPE={}  ..'.format(
            usempi, usecuda, usesp, build_type)
        execute(cmake_cmd)
        execute('make -j 4')

        os.chdir('../python')
