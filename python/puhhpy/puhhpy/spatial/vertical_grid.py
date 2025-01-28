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

# Standard library

# Third-party.
import numpy as np

# Local library


class Vertical_grid_2nd:
    def __init__(self, z_in, zsize_in, remove_ghost=False, dtype=np.float64):
        """
        Calculate vertical grid, identical to definition in MicroHH,
        including one ghost cell at the bottom and top.

        Arguments:
        ----------
        z_in : np.ndarray, shape (1,)
            Array with input full level heights, like in `case_input.nc`.
        zsize : float
            Height of domain top.
        remove_ghost : bool, optional
            Clip off the ghost cells, leaving `ktot` full and `ktot+` half levels.
        dtype : np.dtype, optional
            Output datatype (np.float32 or np.float64) of arrays.
        """

        self.zsize = zsize_in

        self.kmax = z_in.size
        self.kcells = self.kmax+2

        self.kstart = 1
        self.kend = self.kmax+1

        self.z = np.zeros(self.kcells, dtype)
        self.zh = np.zeros(self.kcells, dtype)

        self.dz = np.zeros(self.kcells, dtype)
        self.dzh = np.zeros(self.kcells, dtype)

        # Full level heights
        self.z[self.kstart:self.kend] = z_in
        self.z[self.kstart-1] = -self.z[self.kstart]
        self.z[self.kend] = 2*self.zsize - self.z[self.kend-1]

        # Half level heights
        for k in range(self.kstart+1, self.kend):
            self.zh[k] = 0.5*(self.z[k-1] + self.z[k])
        self.zh[self.kstart] = 0.
        self.zh[self.kend] = self.zsize

        for k in range(1, self.kcells):
            self.dzh[k] = self.z[k] - self.z[k-1];
        self.dzh[self.kstart-1] = self.dzh[self.kstart+1]
        self.dzhi = 1./self.dzh

        for k in range(1, self.kcells-1):
            self.dz[k] = self.zh[k+1] - self.zh[k];
        self.dz[self.kstart-1] = self.dz[self.kstart];
        self.dz[self.kend] = self.dz[self.kend-1];
        self.dzi = 1./self.dz

        if remove_ghost:
            """
            Clip off the ghost cells, leaving `ktot` full levels and `ktot+1` half levels.
            """
            self.z    = self.z   [self.kstart:self.kend]
            self.dz   = self.dz  [self.kstart:self.kend]
            self.dzi  = self.dzi [self.kstart:self.kend]

            self.zh   = self.zh  [self.kstart:self.kend+1]
            self.dzh  = self.dzh [self.kstart:self.kend+1]
            self.dzhi = self.dzhi[self.kstart:self.kend+1]
