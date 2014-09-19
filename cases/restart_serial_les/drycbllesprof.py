import numpy

# set the height
kmax  = 32
zsize = 3200.
dz = zsize / kmax

dthetadz = 0.003

# set the height
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
th = numpy.zeros(numpy.size(z))

for k in range(kmax):
  th[k] = dthetadz*z[k]

# write the data to a file
proffile = open('drycblles.prof','w')
proffile.write('{0:^20s} {1:^20s}\n'.format('z','th'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E}\n'.format(z[k], th[k]))
proffile.close()

# write the data to a file
proffile = open('drycblles_restart.prof','w')
proffile.write('{0:^20s} {1:^20s}\n'.format('z','th'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E}\n'.format(z[k], th[k]))
proffile.close()
