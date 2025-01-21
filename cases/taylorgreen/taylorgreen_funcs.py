import numpy as np
import microhh_tools as mht

class Parse_TaylorGreen:
    def __init__(self, time, visc, path):
        """
        Help class to:
        1. Read MicroHH grid and fields
        2. Create reference fields.
        3. Calculate model error.
        """

        # Read namelist for grid dimensions.
        ini = mht.Read_namelist(f'{path}/taylorgreen.ini')
        itot = ini['grid']['itot']
        jtot = ini['grid']['jtot']
        ktot = ini['grid']['ktot']

        # Read `grid.0000000`
        self.grid = mht.Read_grid(itot, jtot, ktot, filename=f'{path}/grid.0000000')

        dims = self.grid.dim
        self.x  = dims['x']
        self.xh = dims['xh']
        self.z  = dims['z']
        self.zh = dims['zh']

        # Grid is equidistant in all dimensions.
        self.dx = self.x[1] - self.x[0]
        self.dz = self.z[1] - self.z[0]

        """
        Read momentum and pressure fields.
        """
        N = itot * jtot* ktot
        self.u = mht.Read_binary(self.grid, f'{path}/u.xz.100.00000.{time:07d}').read(N).reshape((ktot, jtot, itot))
        self.w = mht.Read_binary(self.grid, f'{path}/w.xz.001.00000.{time:07d}').read(N).reshape((ktot, jtot, itot))
        self.p = mht.Read_binary(self.grid, f'{path}/p.xz.000.00000.{time:07d}').read(N).reshape((ktot, jtot, itot))

        """
        Create reference data.
        """
        self.u_ref = np.zeros((ktot, jtot, itot))
        self.w_ref = np.zeros((ktot, jtot, itot))
        self.p_ref = np.zeros((ktot, jtot, itot))

        for k in range(self.z.size):
            self.u_ref[k, 0, :] = np.sin(2. * np.pi * self.xh) * np.cos(2. * np.pi * self.z[k]) * np.exp(-8. * np.pi**2. * visc * time)
            self.w_ref[k, 0, :] = -np.cos(2. * np.pi * self.x) * np.sin(2. * np.pi * self.zh[k]) * np.exp(-8. * np.pi**2. * visc * time)
            self.p_ref[k, 0, :] = (0.25 * (np.cos(4. * np.pi * self.x) + np.cos(4. * np.pi * self.z[k])) - 0.25) * (np.exp(-8. * np.pi**2. * visc * time)**2.)

        """
        Calculate error.
        """
        self.u_err = 0.
        self.w_err = 0.
        self.p_err = 0.

        for k in range(self.z.size):
            self.u_err = self.u_err + np.sum(self.dx * self.dz * np.abs(self.u[k, :] - self.u_ref[k, :]))
            self.w_err = self.w_err + np.sum(self.dx * self.dz * np.abs(self.w[k, :] - self.w_ref[k, :]))
            self.p_err = self.p_err + np.sum(self.dx * self.dz * np.abs(self.p[k, :] - self.p_ref[k, :]))
