import netCDF4 as nc
import numpy as np

float_type = 'f8'

# set the height
kmax  = 64
zsize = 5.
dz = zsize / kmax

z = np.zeros(kmax)
u = np.zeros(kmax)
v = np.zeros(kmax)
u_geo = np.zeros(kmax)
v_geo = np.zeros(kmax)

u_geo[:] = 1.
v_geo[:] = 0.

z = np.linspace(0.5*dz, zsize-0.5*dz, kmax)

visc = 0.1
fc   = 1.
gamma = (fc / (2.*visc))**.5

u[:] = u_geo[:]
v[:] = v_geo[:]

# analytical solution as the starting profile to reduce run time
#for k in range(kmax):
#  u[k] = u_geo[k]*(1. - exp(-gamma*z[k]) * cos(gamma*z[k]))
#  v[k] = u_geo[k]*(     exp(-gamma*z[k]) * sin(gamma*z[k]))

# write the data to a file
nc_file = nc.Dataset("ekman_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", kmax)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_u = nc_group_init.createVariable("u", float_type, ("z"))
nc_v = nc_group_init.createVariable("v", float_type, ("z"))
nc_u_geo = nc_group_init.createVariable("u_geo", float_type, ("z"))
nc_v_geo = nc_group_init.createVariable("v_geo", float_type, ("z"))

nc_z[:] = z[:]
nc_u[:] = u[:]
nc_v[:] = v[:]
nc_u_geo[:] = u_geo[:]
nc_v_geo[:] = v_geo[:]

nc_file.close()
