import netCDF4 as nc
import numpy as np

float_type = 'f8'

# Get number of vertical levels and size from .ini file
with open('gabls1.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

# Set the height
dz = zsize / kmax
dthetadz = 0.01

# Create vertical profiles
z  = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
th = np.zeros(np.size(z))
u  = np.zeros(np.size(z))
ug = np.zeros(np.size(z))

u [:] = 8.
ug[:] = 8.

for k in range(kmax):
    if(z[k] <= 100.):
        th[k] = 265.
    if(z[k] > 100.):
        th[k] = 265. + dthetadz*(z[k]-100.)

# Surface forcing
time_surface = np.array([0, 32400])
th_sbot = np.array([265, 262.75])

# Save all the input data to NetCDF
nc_file = nc.Dataset('gabls1_input.nc', mode='w', datamodel='NETCDF4', clobber=True)

nc_file.createDimension('z', kmax)
nc_z    = nc_file.createVariable('z', float_type, ('z'))
nc_z[:] = z[:]

# Create a group called "init" for the initial profiles.
nc_group_init = nc_file.createGroup('init')

nc_th = nc_group_init.createVariable('th',    float_type, ('z'))
nc_u  = nc_group_init.createVariable('u'  ,   float_type, ('z'))
nc_ug = nc_group_init.createVariable('u_geo', float_type, ('z'))

nc_th [:] = th[:]
nc_u  [:] = u [:]
nc_ug [:] = ug[:]

# Create a group called "timedep" for the timedep.
nc_group_timedep = nc_file.createGroup('timedep')
nc_group_timedep.createDimension('time_surface', time_surface.size)

nc_time_surface = nc_group_timedep.createVariable('time_surface', float_type, ('time_surface'))
nc_th_sbot  = nc_group_timedep.createVariable('th_sbot', float_type, ('time_surface'))

nc_time_surface[:] = time_surface[:]
nc_th_sbot     [:] = th_sbot     [:]

nc_file.close()
