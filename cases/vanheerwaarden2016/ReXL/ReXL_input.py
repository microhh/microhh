import numpy
from pylab import *
import netCDF4 as nc
float_type = 'f8'

# set the height (ktot = 1024)
kmax = 1024
dn = 1./kmax

n = numpy.linspace(dn, 1.-dn, kmax)

nloc1 = 150.*dn
nbuf1 = 42.*dn

nloc2 = 1024.*dn
nbuf2 = 96.*dn

dz1 = 0.0004  # z0 is calculated as 7.37e-4
dz2 = 0.0009765625
dz3 = 0.009

dzdn1 = dz1/dn
dzdn2 = dz2/dn
dzdn3 = dz3/dn

dzdn = dzdn1 + 0.5*(dzdn2-dzdn1)*(1. + numpy.tanh((n-nloc1)/nbuf1)) + \
    0.5*(dzdn3-dzdn2)*(1. + numpy.tanh((n-nloc2)/nbuf2))

dz = dzdn*dn

z = numpy.zeros(numpy.size(dz))
stretch = numpy.zeros(numpy.size(dz))

z[0] = 0.5*dz[0]
stretch[0] = 1.

for k in range(1, kmax):
    z[k] = z[k-1] + 0.5*(dz[k-1]+dz[k])
    stretch[k] = dz[k]/dz[k-1]

zsize = z[kmax-1] + 0.5*dz[kmax-1]
print('zsize = ', zsize)

b0 = 1.
delta = 4.407731e-3
N2 = 3.

b = numpy.zeros(numpy.size(z))

for k in range(kmax):
    b[k] = N2*z[k]

# write the data to a file
nc_file = nc.Dataset("drycbl_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", ktot)
nc_z = nc_file.createVariable("z" , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_b = nc_group_init.createVariable("b", float_type, ("z"))

nc_z[:] = z[:]
nc_b[:] = b[:]

nc_file.close()

# plot the grid
# figure()
# subplot(131)
# plot(n, z)
# subplot(132)
# plot(n, dz)
# subplot(133)
# plot(n, stretch)
