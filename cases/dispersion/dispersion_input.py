import numpy as np
import netCDF4 as nc

import microhh_tools as mht

"""
Settings
"""
float_type = 'f8'
scalars = ['s1', 's2']     # All scalars.
scalars_outflow = ['s1']   # Scalars with non-periodic lateral BCs.
swtimedep_outflow = False  # Time depedent lateral outflow profiles.

xsize = 12800
ysize = 4800
zsize = 3200

itot = 256
jtot = 96
ktot = 64

endtime = 10800

"""
Source settings.
"""
source_x0 = 0.1*xsize
source_y0 = 0.5*ysize
source_z0 = 100

# Source is a Gaussian "blob":
sigma_x = 25
sigma_y = 25
sigma_z = 25

# Switch between volume and mass.
sw_vmr = False

# Emission strength.
# - with `sw_vmr=True`, units are (kmol s-1)
# - with `sw_vmr=False`, units are (kg kg-1 s-1)
strength = 1

"""
Define vertical grid and input profiles.
"""
dz = zsize / ktot
z  = np.arange(0.5*dz, zsize, dz)
u  = np.zeros_like(z)
th = np.zeros_like(z)
s  = np.zeros_like(z)

# Well mixed profile with jump and stable stratification above.
h = 1000.          # Initial boundary layer depth (m)
th0 = 300.         # Mixed-layer temperature (K)
dth = 2.           # Inversion strength (K)
dthz = 100.        # Inversion depth (m)
dthetadz = 0.003   # Lapse rate free troposphere (K m-1)

for k in range(ktot):
    if(z[k] <= h - 0.5*dthz):
        th[k] = th0
        s[k]  = 0.
    elif(z[k] <= h + 0.5*dthz):
        th[k] = th0 + dth/dthz * (z[k]-(h-0.5*dthz))
    else:
        th[k] = th0 + dth + dthetadz*(z[k]-(h+0.5*dthz))

u[:] = 5

# Outflow profiles
if swtimedep_outflow:
    time_ls = np.array([0, endtime])
    s_lbc = np.zeros((time_ls.size, z.size))
    s_lbc[0,:] = 0
    s_lbc[1,:] = 10
else:
    s_lbc = np.zeros_like(z)

"""
Set/write new namelist.
"""
ini = mht.Read_namelist('dispersion_base.ini')

ini['grid']['itot'] = itot
ini['grid']['jtot'] = jtot
ini['grid']['ktot'] = ktot

ini['grid']['xsize'] = xsize
ini['grid']['ysize'] = ysize
ini['grid']['zsize'] = zsize

ini['buffer']['zstart'] = 0.75*zsize

ini['time']['endtime'] = endtime

ini['fields']['slist'] = scalars
ini['advec']['fluxlimit_list'] = scalars
ini['limiter']['limitlist'] = scalars
ini['boundary']['scalar_outflow'] = scalars_outflow
ini['boundary']['swtimedep_outflow'] = swtimedep_outflow

# Sources
def const_list(value):
    value = int(value) if isinstance(value, bool) else value
    return len(scalars)*[value]

ini['source']['sourcelist'] = scalars
ini['source']['source_x0'] = const_list(source_x0) 
ini['source']['source_y0'] = const_list(source_y0) 
ini['source']['source_z0'] = const_list(source_z0) 
ini['source']['sigma_x'] = const_list(sigma_x)
ini['source']['sigma_y'] = const_list(sigma_y)
ini['source']['sigma_z'] = const_list(sigma_z)
ini['source']['strength'] = const_list(strength)
ini['source']['swvmr'] = const_list(sw_vmr)
ini['source']['line_x'] = const_list(0)
ini['source']['line_y'] = const_list(0)
ini['source']['line_z'] = const_list(0)

# Statistics/crosses/...
scalar_crosses = scalars + [s+'_path' for s in scalars]
ini['cross']['crosslist'] += scalar_crosses
ini['cross']['xy'] = source_z0
ini['cross']['xz'] = source_y0

ini.save('dispersion.ini', allow_overwrite=True)

"""
Create input NetCDF file.
"""
def add_var(name, dims, values, nc_group):
    nc_var = nc_group.createVariable(name, float_type, dims)
    nc_var[:] = values

nc_file = nc.Dataset('dispersion_input.nc', mode='w', datamodel='NETCDF4')
nc_file.createDimension('z', ktot)
add_var('z',  ('z'), z,  nc_file)

nc_init = nc_file.createGroup('init');
add_var('u',  ('z'), u,  nc_init)
add_var('th', ('z'), th, nc_init)

# Same profile with concentration = 0 for all scalars:
for scalar in scalars:
    add_var(scalar,  ('z'), s, nc_init)

if not swtimedep_outflow:
    for scalar in scalars_outflow:
        add_var('{}_inflow'.format(scalar), ('z'), s_lbc, nc_init)
else:
    nc_tdep = nc_file.createGroup('timedep');
    nc_tdep.createDimension('time_ls', time_ls.size)
    add_var('time_ls', ('time_ls'), time_ls, nc_tdep)

    for scalar in scalars_outflow:
        add_var('{}_inflow'.format(scalar), ('time_ls', 'z'), s_lbc, nc_tdep)

nc_file.close()
