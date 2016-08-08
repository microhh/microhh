import numpy as np
import matplotlib.pylab as pl

from microhh_tools import *

# Create stretched grid
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
        print('kmax=%i, zsize=%f'%(kmax,self.zsize))

        if(True):
            pl.figure()
            pl.plot(self.z, self.dz)

nl = Read_namelist()

if(False):
    gr = Grid(64, 40, 5, 0.0002, 0.0005)
    
    u = nl.force.uflux * np.ones(gr.z.size)
    b = gr.z 
    z = gr.z

    replace_namelist_var('zsize', gr.zsize)

if(True):
    dz = nl.grid.zsize / nl.grid.ktot
    z  = np.linspace(0.5*dz, nl.grid.zsize-0.5*dz, nl.grid.ktot)
    u  = np.ones(nl.grid.ktot)*10
    b  = z

# Write the data to a file.
proffile = open('agnesi.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s}\n'.format('z','u','b','ug'))
for k in range(nl.grid.ktot):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E}\n'.format(z[k], u[k], b[k], u[k]))
proffile.close()
