import netCDF4 as nc
import numpy as np
float_type = np.float64

# set the height
# Get number of vertical levels and size from .ini file
with open('sullivan2011.ini') as f:
    for line in f:
        if(line.split('=')[0] == 'ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0] == 'zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

z = np.linspace(0.5 * dz, zsize - 0.5 * dz, kmax)
th = np.empty(z.shape)
u = np.empty(z.shape)
ug = np.empty(z.shape)

for k in range(kmax):
    # temperature
    if(z[k] <= 974.):
        th[k] = 300.
    elif(z[k] <= 1074):
        th[k] = 300. + (z[k] - 974.) * 0.08
    else:
        th[k] = 308. + (z[k] - 1074.) * 0.003

u[:] = 1.
ug[:] = 1.

nc_file = nc.Dataset(
    'sullivan2011_input.nc',
    mode='w',
    datamodel='NETCDF4',
    clobber=True)

nc_file.createDimension('z', kmax)
nc_z = nc_file.createVariable('z', float_type, ('z'))

nc_group_init = nc_file.createGroup('init')
nc_th = nc_group_init.createVariable('th', float_type, ('z'))
nc_u = nc_group_init.createVariable('u', float_type, ('z'))
nc_ug = nc_group_init.createVariable('u_geo', float_type, ('z'))

nc_z[:] = z[:]
nc_u[:] = u[:]
nc_ug[:] = ug[:]
nc_th[:] = th[:]

nc_file.close()
