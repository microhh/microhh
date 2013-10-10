import numpy

# set the height
kmax  = 32
zsize = 3200.
dz = zsize / kmax

dthetadz = 0.003

# set the height
z   = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
s   = numpy.zeros(numpy.size(z))
sls = numpy.zeros(numpy.size(z))
wls = numpy.zeros(numpy.size(z))

for k in range(kmax):
  s  [k] = dthetadz*z[k]
  sls[k] = 2.*(z[k]/zsize - 0.5) / 3600.
  wls[k] = -0.01*(z[k]/zsize)

# write the data to a file
proffile = open('drycblles.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s}\n'.format('z','s','sls', 'wls'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E}\n'.format(z[k], s[k], sls[k], wls[k]))
proffile.close()
