import numpy
import netCDF4 as nc

float_type = 'f8'

# Get number of vertical levels and size from .ini file
with open('moser600.ini') as f:
    for line in f:
        if (line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if (line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

# define the variables
z = numpy.zeros(kmax)
u = numpy.zeros(kmax)
s = numpy.zeros(kmax)

# create non-equidistant grid
alpha = 0.967
for k in range(kmax):
  eta  = -1. + 2.*((k+1)-0.5) / kmax
  z[k] = zsize / (2.*alpha) * numpy.tanh(eta*0.5*(numpy.log(1.+alpha) - numpy.log(1.-alpha))) + 0.5*zsize
  s[k] = z[k]

# create initial parabolic shape
dpdxls = -3.0e-6
visc   =  1.0e-5
for k in range(kmax):
  u[k] = 1./(2.*visc)*dpdxls*(z[k]**2. - zsize*z[k])


# write the data to a file
nc_file = nc.Dataset("moser600_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", kmax)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_u = nc_group_init.createVariable("u", float_type, ("z"))
nc_s = nc_group_init.createVariable("s", float_type, ("z"))

nc_z[:] = z[:]
nc_u[:] = u[:]
nc_s[:] = s[:]

nc_file.close()
