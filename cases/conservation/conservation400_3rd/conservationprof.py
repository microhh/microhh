import numpy

# set the height
kmax  = 32
zsize = 1.
dz    = zsize / kmax

# define the variables
z = numpy.zeros(kmax)
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
s = numpy.zeros(kmax)

# create non-equidistant grid
alpha = 0.967
for k in range(kmax):
  eta  = -1. + 2.*((k+1)-0.5) / kmax
#z[k] = zsize / (2.*alpha) * numpy.tanh(eta*0.5*(numpy.log(1.+alpha) - numpy.log(1.-alpha))) + 0.5*zsize
  s[k] = 2.*z[k]

# write the data to a file
proffile = open('conservation.prof','w')
proffile.write('{0:^22s} {1:^22s}\n'.format('z', 's'))
for k in range(kmax):
  proffile.write('{0:1.16E} {1:1.16E}\n'.format(z[k], s[k]))
proffile.close()

