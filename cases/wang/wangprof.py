import numpy as np
import matplotlib.pylab as pl

pl.close('all')

# Get number of vertical levels and size from .ini file
with open('wang.ini') as f:
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

        zsize = self.z[kmax-1] + 0.5*self.dz[kmax-1]
        print('kmax=%i, zsize=%f'%(kmax,zsize))

#g = Grid(256, 25,  5, 10.0, 10.0)     # isotropic grid
g = Grid(64, 25,  5, 40.0, 40.0)     # isotropic grid
#g = Grid(128,70, 10, 2.0,  40.0)   # close to Wang

th = np.zeros_like(g.z)

for k in range(g.z.size):
    if(g.z[k] < 500):
        th[k] = 300 + 0.5e-3 * g.z[k]
    else:
        th[k] = 300.25 + 9.8e-3 * (g.z[k] - 500)
        
# write the data to a file
proffile = open('wang.prof','w')
proffile.write('{0:^20s} {1:^20s}\n'.format('z','th'))
for k in range(g.z.size):
    proffile.write('{0:1.14E} {1:1.14E}\n'.format(g.z[k], th[k]))
proffile.close()
