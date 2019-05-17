import numpy

# set the height
kmax  = 256
zsize = 10.24
N2 = 0.0004

# define the variables
dz = zsize / kmax
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)

b = numpy.zeros(kmax)
for k in range(kmax):
    b[k] = N2 * z[k]

# write the data to a file
proffile = open('shapiro.prof','w')
proffile.write('{0:^14s} {1:^14s}\n'.format('z', 'b'))
for k in range(kmax):
    proffile.write('{0:1.16E} {1:1.16E}\n'.format(z[k], b[k]))
proffile.close()

