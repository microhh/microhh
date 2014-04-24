import numpy
from pylab import *

# set the height
kmax  = 128
zsize = 1.
dz = zsize / kmax

# set the height
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
u = numpy.zeros(numpy.size(z))
b = numpy.zeros(numpy.size(z))

# set du = 1 and db = 1
# delta is defined as du / (du/dz)max
du = 1.
db = 1.
delta = 0.167
visc = 5.e-5
N2 = db / delta
Ri = N2 * delta**2 / du**2
Re = delta * du / visc
fac = numpy.pi**.5 / delta

print('delta', delta)
print('N2', N2)
print('Ri', Ri)
print('Re', Re)

for k in range(kmax):
  u[k] = 0.5*math.erf(fac*(z[k] - 0.5*zsize))
  b[k] = 0.5*math.erf(fac*(z[k] - 0.5*zsize))

# write the data to a file
proffile = open('shearlayer.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s}\n'.format('z','u','b'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E} {2:1.14E}\n'.format(z[k],u[k],b[k]))
proffile.close()
