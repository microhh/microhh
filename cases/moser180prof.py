import numpy

# set the height
kmax  = 128
zsize = 2.

# define the variables
z = numpy.zeros(kmax)
u = numpy.zeros(kmax)
s = numpy.arange(1., kmax+1., 1.) * zsize/(kmax+1)

# create non-equidistant grid
alpha = 0.967
for k in range(kmax):
  eta  = -1. + 2.*((k+1)-0.5) / kmax
  z[k] = zsize / (2.*alpha) * numpy.tanh(eta*0.5*(numpy.log(1.+alpha) - numpy.log(1.-alpha))) + 0.5*zsize

# create initial parabolic shape
dpdxls = -1.5e-6
visc   =  1.0e-5
for k in range(kmax):
  u[k] = 1./(2.*visc)*dpdxls*(z[k]**2. - zsize*z[k])

# write the data to a file
proffile = open('moser180.prof','w')
proffile.write('{0:^14s} {1:^14s} {2:^14s}\n'.format('z', 'u', 's'))
for k in range(kmax):
  proffile.write('{0:1.8E} {1:1.8E} {2:1.8E}\n'.format(z[k], u[k], s[k]))
proffile.close()

