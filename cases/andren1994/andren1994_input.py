import netCDF4 as nc4
import numpy as np

float_type = "f8"

# Read table A1 paper
A1 = np.loadtxt('andren1994_tableA1')
z  = A1[:,0]
u  = A1[:,1]
v  = A1[:,2]
q2 = A1[:,3]

ug = np.ones(z.size)*10.

# Save all the input data to NetCDF
nc_file = nc4.Dataset("andren1994_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", z.size)
nc_z = nc_file.createVariable("z", float_type, ("z"))
nc_z[:] = z[:]

# Create a group called "init" for the initial profiles.
nc_group_init = nc_file.createGroup("init")

nc_u  = nc_group_init.createVariable("u",     float_type, ("z"))
nc_v  = nc_group_init.createVariable("v",     float_type, ("z"))
nc_ug = nc_group_init.createVariable("u_geo", float_type, ("z"))

nc_u [:] = u [:]
nc_v [:] = v [:]
nc_ug[:] = ug[:]

nc_file.close()
