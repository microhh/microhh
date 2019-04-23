import numpy

# set the height
kmax  = 64
zsize = 0.5

# define the variables
dz = zsize / kmax
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)

# write the data to a file
nc_file = nc.Dataset("taylorgreen_input.nc", mode="w", datamodel="NETCDF4", clobber=False)

nc_file.createDimension("z", kmax)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))
nc_z[:] = z[:]

nc_file.close()
