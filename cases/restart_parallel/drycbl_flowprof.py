import numpy
from scipy.special import erf
#from pylab import *

# set the height
kmax = 512
dn   = 1./kmax

n  = numpy.linspace(dn, 1.-dn, kmax)

nloc1 = 80.*dn
nbuf1 = 16.*dn

nloc2 = 512.*dn
nbuf2 = 128.*dn

dz1 = 0.001
dz2 = 0.002
dz3 = 0.01

dzdn1 = dz1/dn
dzdn2 = dz2/dn
dzdn3 = dz3/dn

dzdn = dzdn1 + 0.5*(dzdn2-dzdn1)*(1. + numpy.tanh((n-nloc1)/nbuf1)) + 0.5*(dzdn3-dzdn2)*(1. + numpy.tanh((n-nloc2)/nbuf2))

dz = dzdn*dn

z       = numpy.zeros(numpy.size(dz))
stretch = numpy.zeros(numpy.size(dz))

z      [0] = 0.5*dz[0]
stretch[0] = 1.

for k in range(1,kmax):
  z      [k] = z[k-1] + 0.5*(dz[k-1]+dz[k])
  stretch[k] = dz[k]/dz[k-1]

zsize = z[kmax-1] + 0.5*dz[kmax-1]
print('zsize = ', zsize)

b0    = 1.
delta = 4.407731e-3
N2    = 3.

s = numpy.zeros(numpy.size(z))

for k in range(kmax):
  s[k] = N2*z[k] + b0*erf(-0.5*z[k]/delta) + b0
  #s[k] = N2*z[k]

# write the data to a file
proffile = open('drycbl_flow.prof','w')
proffile.write('{0:^20s} {1:^20s}\n'.format('z','s'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E}\n'.format(z[k], s[k]))
proffile.close()

#plot the grid
#subplot(131)
#plot(n,z)
#subplot(132)
#plot(n,dz)
#subplot(133)
#plot(n,stretch)
