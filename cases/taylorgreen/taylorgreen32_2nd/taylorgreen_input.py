import netCDF4 as nc
import numpy

float_type = "f8"
# float_type = "f4"

# set the height
kmax  = 16
zsize = 0.5

# define the variables
dz = zsize / kmax
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)

# write the data to a file
nc_file = nc.Dataset("taylorgreen_input.nc", mode="w", datamodel="NETCDF4", clobber=False)

nc_file.createDimension("z", kmax)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))
nc_z[:] = z[:]

nc_group_init = nc_file.createGroup("init");

nc_file.close()
