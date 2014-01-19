import numpy

# set the height
kmax  = 512
zsize = 0.5
dz = zsize / kmax

# set the height
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
s = numpy.zeros(numpy.size(z))
s[0:kmax/2] = 1.

# write the data to a file
proffile = open('rayleightaylor.prof','w')
proffile.write('{0:^20s} {1:^20s}\n'.format('z','s'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E}\n'.format(z[k], s[k]))
proffile.close()
