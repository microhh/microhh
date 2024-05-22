import netCDF4 as nc4
import numpy as np

# Model datatype
float_type = "f8"

# Set the height
kmax  = 512
zsize = 0.5

dz = zsize / kmax

# Define the variables
z = np.zeros(kmax)
b = np.zeros(kmax)

# Create non-equidistant grid
alpha = 0.7
for k in range(kmax):
    eta  = -1. + 2.*((k+1)-0.5) / kmax
    z[k] = zsize / (2.*alpha) * np.tanh(eta*0.5*(np.log(1.+alpha) - np.log(1.-alpha))) + 0.5*zsize

# Write input NetCDF file
nc_file = nc4.Dataset('rayleighbenard_input.nc', mode='w', datamodel='NETCDF4', clobber=True)
nc_file.createDimension('z', kmax)
nc_group_init = nc_file.createGroup('init');

nc_z = nc_file.createVariable('z' , float_type, ('z'))
nc_z[:] = z[:]

nc_file.close()
