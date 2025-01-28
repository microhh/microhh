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
import puhhpy.constants as cst
from puhhpy.spatial import Vertical_grid_2nd

from .base_thermo import exner, virtual_temperature, sat_adjust


class Basestate_moist:
    def __init__(self, thl, qt, pbot, z, zsize, remove_ghost=False, dtype=np.float64):
        """
        Calculate moist thermodynamic base state,
        using a procedure identical to MicroHH.
        """
        gd = Vertical_grid_2nd(z, zsize, dtype=dtype)

        self.gd = gd
        self.remove_ghost = remove_ghost
        self.dtype = dtype
        
        self.pref = np.zeros(gd.kcells)
        self.prefh = np.zeros(gd.kcells)

        self.rho = np.zeros(gd.kcells)
        self.rhoh = np.zeros(gd.kcells)

        self.thv = np.zeros(gd.kcells)
        self.thvh = np.zeros(gd.kcells)

        self.ex = np.zeros(gd.kcells)
        self.exh = np.zeros(gd.kcells)

        # Add ghost cells to input profiles
        self.thl = np.zeros(gd.kcells, dtype)
        self.qt  = np.zeros(gd.kcells, dtype)

        self.thl[gd.kstart:gd.kend] = thl
        self.qt [gd.kstart:gd.kend] = qt

        # Calculate surface and domain top values.
        self.thl0s = self.thl[gd.kstart] - gd.z[gd.kstart] * (self.thl[gd.kstart+1] - self.thl[gd.kstart]) * gd.dzhi[gd.kstart+1]
        self.qt0s  = self.qt[gd.kstart]  - gd.z[gd.kstart] * (self.qt [gd.kstart+1] - self.qt [gd.kstart]) * gd.dzhi[gd.kstart+1]

        self.thl0t = self.thl[gd.kend-1] + (gd.zh[gd.kend] - gd.z[gd.kend-1]) * (self.thl[gd.kend-1]-self.thl[gd.kend-2])*gd.dzhi[gd.kend-1]
        self.qt0t  = self.qt[gd.kend-1]  + (gd.zh[gd.kend] - gd.z[gd.kend-1]) * (self.qt [gd.kend-1]- self.qt[gd.kend-2])*gd.dzhi[gd.kend-1]

        # Set the ghost cells for the reference temperature and moisture
        self.thl[gd.kstart-1]  = 2.*self.thl0s - self.thl[gd.kstart];
        self.thl[gd.kend]      = 2.*self.thl0t - self.thl[gd.kend-1];

        self.qt[gd.kstart-1]   = 2.*self.qt0s  - self.qt[gd.kstart];
        self.qt[gd.kend]       = 2.*self.qt0t  - self.qt[gd.kend-1];

        # Calculate profiles. 
        self.prefh[gd.kstart] = pbot
        self.exh[gd.kstart] = exner(pbot)

        _, ql, qi, _ = sat_adjust(self.thl0s, self.qt0s, pbot)

        self.thvh[gd.kstart] = virtual_temperature(
                self.exh[gd.kstart], self.thl0s, self.qt0s, ql, qi)
        self.rhoh[gd.kstart] = pbot / (cst.Rd * self.exh[gd.kstart] * self.thvh[gd.kstart])

        # Calculate the first full level pressure
        self.pref[gd.kstart] = \
            self.prefh[gd.kstart] * np.exp(-cst.grav * gd.z[gd.kstart] / (cst.Rd * self.exh[gd.kstart] * self.thvh[gd.kstart]))

        for k in range(gd.kstart+1, gd.kend+1):
            # 1. Calculate remaining values (thv and rho) at full-level[k-1]
            self.ex[k-1] = exner(self.pref[k-1])
            _, ql, qi, _ = sat_adjust(self.thl[k-1], self.qt[k-1], self.pref[k-1], self.ex[k-1])
            self.thv[k-1] = virtual_temperature(self.ex[k-1], self.thl[k-1], self.qt[k-1], ql, qi)
            self.rho[k-1] = self.pref[k-1] / (cst.Rd * self.ex[k-1] * self.thv[k-1])
            
            # 2. Calculate pressure at half-level[k]
            self.prefh[k] = self.prefh[k-1] * np.exp(-cst.grav * gd.dz[k-1] / (cst.Rd * self.ex[k-1] * self.thv[k-1]))
            self.exh[k] = exner(self.prefh[k])
            
            # 3. Use interpolated conserved quantities to calculate half-level[k] values
            thli = 0.5*(self.thl[k-1] + self.thl[k])
            qti  = 0.5*(self.qt[k-1] + self.qt[k])
            _, qli, qii, _ = sat_adjust(thli, qti, self.prefh[k], self.exh[k])
            
            self.thvh[k] = virtual_temperature(self.exh[k], thli, qti, qli, qii)
            self.rhoh[k] = self.prefh[k] / (cst.Rd * self.exh[k] * self.thvh[k])
            
            # 4. Calculate pressure at full-level[k]
            self.pref[k] = self.pref[k-1] * np.exp(-cst.grav * gd.dzh[k] / (cst.Rd * self.exh[k] * self.thvh[k]))

        self.pref[gd.kstart-1] = 2. * self.prefh[gd.kstart] - self.pref[gd.kstart]

        if remove_ghost:
            """
            Strip off the ghost cells, to leave `ktot` full levels and `ktot+1` half levels.
            """
            self.pref = self.pref[gd.kstart:gd.kend]
            self.prefh = self.prefh[gd.kstart:gd.kend+1]

            self.rho = self.rho[gd.kstart:gd.kend]
            self.rhoh = self.rhoh[gd.kstart:gd.kend+1]

            self.thv = self.thv[gd.kstart:gd.kend]
            self.thvh = self.thvh[gd.kstart:gd.kend+1]

            self.ex = self.ex[gd.kstart:gd.kend]
            self.exh = self.exh[gd.kstart:gd.kend+1]

            self.thl = self.thl[gd.kstart:gd.kend]
            self.qt  = self.qt[gd.kstart:gd.kend]


    def to_binary(self, grid_file):
        """
        Save base state in format required by MicroHH.
        """

        if self.remove_ghost:
            rho = self.rho
            rhoh = self.rhoh
        else:
            gd = self.gd
            rho = self.rho[gd.kstart:gd.kend]
            rhoh = self.rhoh[gd.kstart:gd.kend+1]

        bs = np.concatenate((rho, rhoh)).astype(self.dtype)
        bs.tofile(grid_file)