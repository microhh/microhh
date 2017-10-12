import numpy as np
import matplotlib.pylab as pl
from scipy.special import erf

from microhh_tools import *

# Create stretched grid
class Stretched_grid:
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
        print('Grid: kmax=%i, zsize=%f'%(kmax,self.zsize))
        
    def plot(self):
        pl.figure()
        pl.plot(self.dz, self.z, '-x')
        pl.xlabel('dz (m)')
        pl.ylabel('z (m)')      

# Create stretched grid
gr = Stretched_grid(64, 40, 10, 2.0, 10.0)
#gr.plot()

# Write back vertical extent domain
replace_namelist_value('zsize', gr.zsize)
replace_namelist_value('ktot', gr.z.size)

# Create initial vertical profiles
u = 0.09/0.4 * np.log(gr.z)
v = np.zeros_like(gr.z)

# Write to .prof file as input for MicroHH
variables = {'z':gr.z, 'u':u, 'v':v}
write_output('ib_poly.prof', variables, gr.z.size)


# Objects...
itot=64
jtot=64
xsize=320
ysize=320

dx = xsize/itot
dy = ysize/jtot

x = np.arange(0.5*dx, xsize, dx)
y = np.arange(0.5*dy, ysize, dy)

xh = np.arange(0,xsize+0.1,dx)
yh = np.arange(0,ysize+0.1,dy)

# Square
x1 = np.array([81,141,141,81,81])
y1 = np.array([81,81,141,141,81])

# Rectangle
x2 = np.array([182,261,221,182])
y2 = np.array([202,202,142,202])

# Approx. circle
r = 40
theta = np.linspace(0,2*np.pi,8)
x3 = 110 + r * np.sin(theta)
y3 = 225 + r * np.cos(theta)

def printout(x,y,z):
    print('{0:.1f}, '.format(z), end='')
    for i in range(x.size):
        end = '' if i==x.size-1 else ', '
        print('{0:.1f}, {1:.1f}{2:}'.format(x[i],y[i],end), end='')
    print('\n')

printout(x1,y1,20)
printout(x2,y2,30)
printout(x3,y3,40)


pl.close('all')

pl.figure()
ax=pl.subplot(111)

pl.plot(x1,y1)
pl.plot(x2,y2)
pl.plot(x3,y3)

pl.xlim(0,xsize)
pl.ylim(0,ysize)
ax.set_xticks(xh)
ax.set_yticks(yh)
#pl.grid()


