import netCDF4 as nc4
import numpy as np

# Model datatype
float_type = "f8"

# Set the height
kmax  = 512
zsize = 0.5
dz = zsize / kmax

# Set the profiles
z = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
b = np.zeros(np.size(z))
b[0:int(kmax/2)] = 1.

# Write input NetCDF file
nc_file = nc4.Dataset('rayleightaylor_input.nc', mode='w', datamodel='NETCDF4', clobber=True)
nc_file.createDimension('z', kmax)

nc_z = nc_file.createVariable('z' , float_type, ('z'))
nc_z[:] = z[:]

nc_group_init = nc_file.createGroup('init');
nc_b = nc_group_init.createVariable('b' , float_type, ('z'))
nc_b[:] = b[:]

nc_file.close()
