import numpy as np
import netCDF4 as nc

float_type = "f8"

# Get number of vertical levels and size from .ini file
with open('conservation.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])
dz    = zsize / kmax

# define the variables
z = np.zeros(kmax)
z = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
s = np.zeros(kmax)

# create non-equidistant grid
alpha = 0.967
for k in range(kmax):
  eta  = -1. + 2.*((k+1)-0.5) / kmax
#z[k] = zsize / (2.*alpha) * numpy.tanh(eta*0.5*(numpy.log(1.+alpha) - numpy.log(1.-alpha))) + 0.5*zsize
  s[k] = 2.*z[k]


nc_file = nc.Dataset("conservation_input.nc", mode="w", datamodel="NETCDF4", clobber=True)
nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, ("z"))

nc_group_init = nc_file.createGroup("init")
nc_s  = nc_group_init.createVariable("s"   , float_type, ("z"))


nc_z[:] = z[:]
nc_s[:] = s[:]

nc_file.close()
