import numpy
from scipy.special import erf

# set the height
kmax  = 512

# define the variables
# create uniform grid
#zref = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
ztmp = numpy.zeros(kmax+1)
z    = numpy.zeros(kmax)
dz   = numpy.zeros(kmax)
s    = numpy.zeros(kmax)

zrefsize = 0.5
zref = numpy.linspace(0., zrefsize, kmax+1)

# create non-equidistant grid consisting of tanh
dzratio  = 2.
dzloc    = 0.171875
dzdelta  = 0.015625

# use integrated tanh function
for k in range(kmax+1):
  ztmp[k] = zref[k] + (dzratio-1.)*dzdelta*numpy.log( numpy.exp((zref[k]-dzloc)/dzdelta)+1.)

for k in range(kmax):
  z [k] = 0.5*(ztmp[k]+ztmp[k+1])
  dz[k] = (ztmp[k+1]-ztmp[k])

b0    = 1.
delta = 4.407731e-3

for k in range(kmax):
  s[k] = z[k] + b0*erf(-0.5*z[k]/delta) + b0

# write the data to a file
proffile = open('drycbl.prof','w')
proffile.write('{0:^14s} {1:^14s}\n'.format('z','s'))
for k in range(kmax):
  proffile.write('{0:1.20E} {1:1.20E}\n'.format(z[k], s[k]))
proffile.close()

