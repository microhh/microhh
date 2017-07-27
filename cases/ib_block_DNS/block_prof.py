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
        pl.subplot(121)
        pl.plot(self.z, self.dz, '-x')
        pl.xlabel('z (m)')
        pl.ylabel('dz (m)') 

        pl.subplot(122)
        pl.plot(self.z, '-x', label='full')
        pl.plot(0.5*(self.z[1:]+self.z[:-1]), '-o', label='hafl')
        pl.xlabel('z (m)')
        pl.ylabel('dz (m)') 
        pl.legend(frameon=False)

class Grid_uni:
    def __init__(self, kmax, zsize):
        dz = zsize / kmax
        self.z = np.arange(0.5*dz, zsize, dz)

if __name__ == "__main__":

    if (True):
        # Read the .ini file
        ini = Read_namelist()

        # Create stretched grid
        #grid = Grid(96, 70, 5, 0.02, 0.06)
        #grid = Grid(96, 40, 5, 0.01, 0.04)
        grid = Grid(120, 75, 5, 0.017, 0.04)
        grid.plot()

        # Write `zsize` and `ktot` back to .ini file
        replace_namelist_value('zsize', grid.zsize)
        replace_namelist_value('ktot',  grid.kmax)

        #grid = Grid_uni(96, 3)

    # Create initial profiles:
    z = grid.z
    #u = ini['force']['uflux'] * np.ones(z.size)
    u = 0.1 * np.ones(z.size)

    # Write the data to a .prof file for MicroHH
    proffile = open('block.prof','w')
    proffile.write('{0:^20s} {1:^20s}\n'.format('z', 'u'))
    for k in range(z.size):
        proffile.write('{0:1.14E} {1:1.14E}\n'.format(z[k], u[k]))
    proffile.close()
