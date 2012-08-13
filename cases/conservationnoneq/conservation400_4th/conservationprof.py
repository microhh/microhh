import numpy

# set the height
kmax    = 32
zsize   = 1.
dz      = zsize / kmax
stretch = 1.02

# define the variables
z = numpy.zeros(kmax)
s = numpy.zeros(kmax)

# non-equidistant grid
dz = zsize*(1-stretch)/(1-stretch**kmax)
z[0] = 0.5*dz
for k in range(1,kmax):
  z[k] = z[k-1] + 0.5*dz
  dz   *= stretch
  z[k] += 0.5*dz
  s[k] =  2.*z[k]

# write the data to a file
proffile = open('conservation.prof','w')
proffile.write('{0:^22s} {1:^22s}\n'.format('z', 's'))
for k in range(kmax):
  proffile.write('{0:1.16E} {1:1.16E}\n'.format(z[k], s[k]))
proffile.close()

