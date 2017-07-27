import numpy as np
import matplotlib.pylab as pl

from microhh_tools import *    # available in MICROHH_DIR/python; copy to this directory

# Class to create a smooth stretched vertical grid
class Grid:
    def __init__(self, kmax, nloc1, nbuf1, dz1, dz2):
        dn         = 1./kmax
        n          = np.linspace(dn, 1.-dn, kmax)
        nloc1     *= dn
        nbuf1     *= dn
        dzdn1      = dz1/dn
        dzdn2      = dz2/dn

        dzdn       = dzdn1 + 0.5*(dzdn2-dzdn1)*(1. + np.tanh((n-nloc1)/nbuf1))
        self.dz    = dzdn*dn

        self.kmax  = kmax
        self.z     = np.zeros(self.dz.size)
        stretch    = np.zeros(self.dz.size)

        self.z[0]  = 0.5*self.dz[0]
        stretch[0] = 1.

        for k in range(1, self.kmax):
              self.z [k] = self.z[k-1] + 0.5*(self.dz[k-1]+self.dz[k])
              stretch[k] = self.dz[k]/self.dz[k-1]

        self.zsize = self.z[kmax-1] + 0.5*self.dz[kmax-1]
        print('kmax={0:}, zsize={1:.4f}'.format(kmax,self.zsize))

    def plot(self):
        pl.figure()
        pl.plot(self.z, self.dz, '-x')
        pl.xlabel('z (m)')
        pl.ylabel('dz (m)') 

if __name__ == "__main__":
    # Read the .ini file
    ini = Read_namelist()

    # Create stretched grid
    grid = Grid(96, 40, 5, 0.0004, 0.0007)
    #grid = Grid(128, 40, 5, 0.0002, 0.0005)
    #grid = Grid(256, 122, 10, 0.0001, 0.0003)
    #grid = Grid(384, 180, 20, 0.00006, 0.00021)

    grid.plot()

    # Write `zsize` and `ktot` back to .ini file
    replace_namelist_value('zsize', grid.zsize)
    replace_namelist_value('ktot',  grid.kmax)

    # Create initial profiles:
    z = grid.z
    u = ini['force']['uflux'] * np.ones(z.size)
    s = z

    # Write the data to a .prof file for MicroHH
    proffile = open('sine.prof','w')
    proffile.write('{0:^20s} {1:^20s} {2:^20s}\n'.format('z', 'u', 's'))
    for k in range(z.size):
        proffile.write('{0:1.14E} {1:1.14E} {2:1.14E}\n'.format(z[k], u[k], s[k]))
    proffile.close()
