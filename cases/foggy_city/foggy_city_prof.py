import numpy as np
import matplotlib.pylab as pl
from microhh_tools import *

pl.close('all')

# Get number of vertical levels and size from .ini file
with open('foggy_city.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

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

    def plot(self):
        pl.figure()
        pl.plot(self.z, self.dz)

#g = Grid(256, 25,  5, 10.0, 10.0)     # isotropic grid
g = Grid(64, 20,  5, 40.0, 40.0)     # isotropic grid
#g = Grid(128, 70, 10, 2.0,  40.0)   # close to Wang
#g = Grid(128, 50, 20, 4.0,  20.0)   # close to Wang

g.plot()

th = np.zeros_like(g.z)
qt = np.zeros_like(g.z)

for k in range(g.z.size):
    if(g.z[k] < 200):
        th[k] = 280
        qt[k] = 6.5e-3
    else:
        th[k] = 280    + 3e-3 * (g.z[k] - 200)
        qt[k] = 6.5e-3 - 2e-6 * (g.z[k] - 200)

# BvS test
u  = np.zeros(g.z.size)
v  = u
ug = u
vg = v

# write the data to a file
proffile = open('foggy_city.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s} {5:^20s} {6:^20s}\n'.format('z','thl','qt','u','v','ug','vg'))
for k in range(g.z.size):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E} {6:1.14E}\n'.format(g.z[k], th[k], qt[k], u[k], v[k], ug[k], vg[k]))
proffile.close()

# Write 'zsize' to namelist
replace_namelist_var('zsize', g.zsize)



