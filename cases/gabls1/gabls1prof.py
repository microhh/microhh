import numpy

# set the height
kmax  = 64
zsize = 400.
dz = zsize / kmax

dthetadz = 0.01

# set the height
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
th = numpy.zeros(numpy.size(z))
u = numpy.zeros(numpy.size(z))
ug = numpy.zeros(numpy.size(z))

u [:] = 8.
ug[:] = 8.

for k in range(kmax):
  if(z[k] <= 100.):
    th[k] = 265.
  if(z[k] > 100.):
    th[k] = 265. + dthetadz*(z[k]-100.)

# write the data to a file
proffile = open('gabls1.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s}\n'.format('z','th','u','ug'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} \n'.format(z[k], th[k], u[k], ug[k]))
proffile.close()

# write surface temperature
timefile = open('gabls1.time','w')
timefile.write('t     sbot[th] \n')
timefile.write('0     265 \n')
timefile.write('32400 262.75 \n')
timefile.close()
