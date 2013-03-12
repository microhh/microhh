import numpy

# set the height
kmax  = 64
zsize = 3200.
dz = zsize / kmax

dthetadz = 0.003

# set the height
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
s = numpy.zeros(numpy.size(z))


for k in range(kmax):
  s[k] = dthetadz*z[k]

# write the data to a file
proffile = open('drycblles.prof','w')
proffile.write('{0:^20s} {1:^20s}\n'.format('z','s'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E}\n'.format(z[k], s[k]))
proffile.close()
