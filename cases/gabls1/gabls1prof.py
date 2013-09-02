import numpy

# set the height
kmax  = 32
zsize = 400.
dz = zsize / kmax

dthetadz = 0.01

# set the height
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
s = numpy.zeros(numpy.size(z))
u = numpy.zeros(numpy.size(z))
v = numpy.zeros(numpy.size(z))
ug = numpy.zeros(numpy.size(z))
vg = numpy.zeros(numpy.size(z))

u [:] = 8.
ug[:] = 8.

for k in range(kmax):
  if(z[k] > 100.):
    s[k] = dthetadz*(z[k]-100.)

# write the data to a file
proffile = open('gabls1.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s} {5:^20s}\n'.format('z','s','u','v','ug','vg'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E}\n'.format(z[k], s[k], u[k], v[k], ug[k], vg[k]))
proffile.close()

