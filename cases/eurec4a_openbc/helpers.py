#
# MicroHH
# Copyright (c) 2011-2023 Chiel van Heerwaarden
# Copyright (c) 2011-2023 Thijs Heus
# Copyright (c) 2014-2023 Bart van Stratum
# 
# This file is part of MicroHH
# 
# MicroHH is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# MicroHH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
#

import numpy as np

def check_grid_decomposition(itot, jtot, ktot, npx, npy):
    """
    Check whether grid / MPI decomposition is valid
    """

    err = False
    if itot%npx != 0:
        print('ERROR in grid: itot%npx != 0')
        err = True

    if itot%npy != 0:
        print('ERROR in grid: itot%npy != 0')
        err = True

    if jtot%npx != 0 and npy > 1:
        print('ERROR in grid: jtot%npx != 0')
        err = True

    if jtot%npy != 0:
        print('ERROR in grid: jtot%npy != 0')
        err = True

    if ktot%npx != 0:
        print('ERROt in grid: ktot%npx != 0')
        err = True

    if err:
        print('Grid: itot={}, jtot={}, ktot={}, npx={}, npy={}'.format(
            itot, jtot, ktot, npx, npy))
        raise Exception('Invalid grid configuration!')
    else:
        print('Grid: itot={}, jtot={}, ktot={}, npx={}, npy={}: OKAY!'.format(
            itot, jtot, ktot, npx, npy))


class Grid_stretched_manual:
    def __init__(self, kmax, dz0, heights, factors):
        self.kmax = kmax
        self.dz0  = dz0

        self.z = np.zeros(kmax)
        self.zh = np.zeros(kmax+1)
        self.dz = np.zeros(kmax)
        self.zsize = None

        self.z[0]  = dz0/2.
        self.dz[0] = dz0

        def index(z, goal):
            return np.where(z-goal>0)[0][0]-1

        for k in range(1, kmax):
            self.dz[k] = self.dz[k-1] * factors[index(heights, self.z[k-1])]
            self.z[k] = self.z[k-1] + self.dz[k]

        self.zsize = self.z[kmax-1] + 0.5*self.dz[kmax-1]

        self.zh[1:-1] = self.z[1:] - self.z[:-1]
        self.zh[-1] = self.zsize

