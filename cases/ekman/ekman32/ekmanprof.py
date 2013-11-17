import numpy
from scipy.special import erf
from pylab import *

# set the height
kmax  = 32
zsize = 5.
dz = zsize / kmax

z = zeros(kmax)
u = zeros(kmax)
v = zeros(kmax)
s = zeros(kmax)
ug = zeros(kmax)
vg = zeros(kmax)

ug[:] = 1.
vg[:] = 0.

z = linspace(0.5*dz, zsize-0.5*dz, kmax)

# analytical solution as the starting profile to reduce run time
visc = 0.1
fc   = 1.
gamma = (fc / (2.*visc))**.5

for k in range(kmax):
  u[k] = ug[k]*(1. - exp(-gamma*z[k]) * cos(gamma*z[k]))
  v[k] = ug[k]*(     exp(-gamma*z[k]) * sin(gamma*z[k]))

# write the data to a file
proffile = open('ekman.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s}\n'.format('z','u','v','ug','vg','s'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E}\n'.format(z[k], u[k], v[k], ug[k], vg[k], s[k]))
proffile.close()

