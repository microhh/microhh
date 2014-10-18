###############################
# Creates stretched grid for  #
# using geometric progression #
###############################
import numpy

# set vertical grid info
kmax = 512
zmax = 1.
dn   = zmax/kmax
n = numpy.linspace(dn, 1., kmax)

# iterate to solve for r in tha sum of geometric progression
dz1 = 0.001
r   = 1.01
for t in range(0,50):
	r = ( 1 - (zmax/dz1)*(1-r) )**(1.0/kmax)

print r

# create array of dz values
dz = numpy.zeros(kmax)

for k in range(0,kmax):
	dz[k] = dz1 * r**k

# create array of heigh values
z = numpy.zeros(numpy.size(dz))

z[0] = 0.5*dz[0]

for k in range(1,kmax):
  z[k] = z[k-1] + 0.5*(dz[k-1]+dz[k])

zsize = z[kmax-1] + 0.5*dz[kmax-1]
print('zsize = ', zsize)

s = numpy.zeros(numpy.size(z))

# write the data to a file
proffile = open('prandtlslope.prof','w')
proffile.write('{0:^20s} {1:^20s}\n'.format('z','s'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E}\n'.format(z[k], s[k]))
proffile.close()
