import numpy
#from scipy.special import erf
from pylab import *

# set the height (ktot = 512)
ktot  = 384
itot  = 512
xsize = 1.

dn   = xsize/itot

n = numpy.linspace(dn, 1.-dn, ktot)

nloc1 = 80.*dn
nbuf1 = 24.*dn

nloc2 = 512.*dn
nbuf2 = 72.*dn

dz1 = 0.5*dn
dz2 = dn
dz3 = 5.*dn

# set the height (ktot = 1024)
"""
ktot = 1024
dn   = 1./ktot

n  = numpy.linspace(dn, 1.-dn, ktot)

nloc1 = 150.*dn
nbuf1 = 32.*dn

nloc2 = 1024.*dn
nbuf2 = 192.*dn

dz1 = 0.0004 #z0 is calculated as 7.37e-4
dz2 = 0.0009765625
dz3 = 0.008
"""

dzdn1 = dz1/dn
dzdn2 = dz2/dn
dzdn3 = dz3/dn

dzdn = dzdn1 + 0.5*(dzdn2-dzdn1)*(1. + numpy.tanh((n-nloc1)/nbuf1)) + 0.5*(dzdn3-dzdn2)*(1. + numpy.tanh((n-nloc2)/nbuf2))

dz = dzdn*dn

z       = numpy.zeros(numpy.size(dz))
stretch = numpy.zeros(numpy.size(dz))

z      [0] = 0.5*dz[0]
stretch[0] = 1.

for k in range(1,ktot):
  z      [k] = z[k-1] + 0.5*(dz[k-1]+dz[k])
  stretch[k] = dz[k]/dz[k-1]

zsize = z[ktot-1] + 0.5*dz[ktot-1]
print('zsize = ', zsize)

b0    = 1.
delta = 4.407731e-3
N2    = 3.

b = numpy.zeros(numpy.size(z))

for k in range(ktot):
  #b[k] = N2*z[k] + b0*erf(-0.5*z[k]/delta) + b0
  b[k] = N2*z[k]

# write the data to a file
proffile = open('thermal.prof','w')
proffile.write('{0:^20s} {1:^20s}\n'.format('z','b'))
for k in range(ktot):
  proffile.write('{0:1.14E} {1:1.14E}\n'.format(z[k], b[k]))
proffile.close()

#plot the grid
figure()
subplot(131)
plot(n,z)
subplot(132)
plot(n,dz)
subplot(133)
plot(n,stretch)
