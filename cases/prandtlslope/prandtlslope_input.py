import numpy as np
import netCDF4 as nc

float_type = "f8"
# float_type = "f4"

# set the height (ktot = 512)
kmax = 512
dn   = 1./kmax
n = np.linspace(dn, 1., kmax)

dz1 = 0.001
r   = 1.01
for t in range(0,50):
	r = ( 1 - (1./dz1)*(1-r) )**(1.0/kmax)

dz = np.zeros(kmax)

for k in range(0,kmax):
	dz[k] = dz1 * r**k

z = np.zeros(dz.size)

z[0] = 0.5*dz[0]

for k in range(1,kmax):
  z[k] = z[k-1] + 0.5*(dz[k-1]+dz[k])

zsize = z[kmax-1] + 0.5*dz[kmax-1]

b = np.zeros(z.size)

nc_file = nc.Dataset("prandtlslope_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", kmax)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_b = nc_group_init.createVariable("b", float_type, ("z"))

nc_z[:] = z[:]
nc_b[:] = b[:]

nc_file.close()

"""
#plot the grid
figure()
subplot(131)
plot(n,z)
subplot(132)
plot(n,dz)
subplot(133)
plot(n,stretch)
"""
