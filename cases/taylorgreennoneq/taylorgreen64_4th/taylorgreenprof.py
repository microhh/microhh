import numpy

# set the values
kmax  = 32
zsize = 0.5
stretch = 1.03

# equidistant grid
dz = zsize / kmax
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)

# non-equidistant grid
dz = zsize*(1-stretch)/(1-stretch**kmax)
z[0] = 0.5*dz
for k in range(1,kmax):
  z[k] = z[k-1] + 0.5*dz
  dz   *= stretch
  z[k] += 0.5*dz


# write the data to a file
proffile = open('taylorgreen.prof','w')
proffile.write('{0:^22s}\n'.format('z'))
for k in range(kmax):
  proffile.write('{0:1.16E}\n'.format(z[k]))
proffile.close()

