import numpy as np
import netCDF4 as nc

float_type = "f8"
# float_type = "f4"

# Get number of vertical levels and size from .ini file
with open('eady.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

dthetadz = 0.006

# set the height
z = np.linspace(0.5*dz, zsize-0.5*dz, kmax)

fc = 1.e-2
dudz = 1e-2

# Linearly stratified profile.
th = 300. + dthetadz*z
u = dudz*z
u_geo = u.copy()
#print("dthetady_ls = {0}".format(-dudz*fc))

# Write the data to a file.
nc_file = nc.Dataset("eady_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", kmax)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_u     = nc_group_init.createVariable("u"    , float_type, ("z"))
nc_u_geo = nc_group_init.createVariable("u_geo", float_type, ("z"))
nc_th    = nc_group_init.createVariable("th"   , float_type, ("z"))

nc_z    [:] = z    [:]
nc_u    [:] = u    [:]
nc_u_geo[:] = u_geo[:]
nc_th   [:] = th   [:]

nc_file.close()
