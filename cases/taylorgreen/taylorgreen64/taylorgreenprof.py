import numpy

# set the height
kmax  = 32
zsize = 0.5

# define the variables
dz = zsize / kmax
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)

# write the data to a file
proffile = open('taylorgreen.prof','w')
proffile.write('{0:^14s}\n'.format('z'))
for k in range(kmax):
  proffile.write('{0:1.8E}\n'.format(z[k]))
proffile.close()

