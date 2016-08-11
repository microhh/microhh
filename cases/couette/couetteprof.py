import numpy

# set the height
kmax  = 32
zsize = 2.

# define the variables
z = numpy.zeros(kmax)
u = numpy.zeros(kmax)

# create non-equidistant grid
alpha = 0.967
for k in range(kmax):
    eta  = -1. + 2.*((k+1)-0.5) / kmax
    z[k] = zsize / (2.*alpha) * numpy.tanh(eta*0.5*(numpy.log(1.+alpha) - numpy.log(1.-alpha))) + 0.5*zsize

# write the data to a file
proffile = open('couette.prof','w')
proffile.write('{0:^14s} {1:^14s}\n'.format('z', 'u'))
for k in range(kmax):
    proffile.write('{0:1.8E} {1:1.8E}\n'.format(z[k], u[k]))
proffile.close()

