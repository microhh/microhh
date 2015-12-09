import numpy as np
import matplotlib.pylab as pl

# set the height
kmax  = 64
zsize = 2.

class Grid_fields:
    def __init__(self, kmax, alpha, zsize=2):
        self.kmax = kmax

        # define the variables
        self.z = np.zeros(kmax)
        self.u = np.zeros(kmax)
        self.s = np.zeros(kmax)
        
        # create non-equidistant grid
        #alpha = 0.967
        #alpha = 0.99
        for k in range(kmax):
            eta       = -1. + 2.*((k+1)-0.5) / kmax
            self.z[k] = zsize / (2.*alpha) * np.tanh(eta*0.5*(np.log(1.+alpha) - np.log(1.-alpha))) + 0.5*zsize
            self.s[k] = self.z[k]

        # create initial parabolic shape
        dpdxls = -1.5e-6
        visc   =  1.0e-5
        for k in range(kmax):
            self.u[k] = 1./(2.*visc)*dpdxls*(self.z[k]**2. - zsize*self.z[k])

default = Grid_fields(128, 0.967)
les     = Grid_fields(48,  0.967)
case    = les 


# write the data to a file
proffile = open('moser180.prof','w')
proffile.write('{0:^14s} {1:^14s} {2:^14s}\n'.format('z', 'u', 's'))
for k in range(case.kmax):
    proffile.write('{0:1.8E} {1:1.8E} {2:1.8E}\n'.format(case.z[k], case.u[k], case.s[k]))
proffile.close()

if(True):
    default.dz = default.z[1:]-default.z[:-1]
    les.dz = les.z[1:]-les.z[:-1]
    
    pl.figure()
    pl.plot(default.z[:-1], default.dz, label='default')
    pl.plot(les.z[:-1], les.dz, label='les')
    pl.legend(frameon=False)

    pl.figure()
    pl.plot(les.u, les.z)
