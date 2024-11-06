import numpy as np
import netCDF4 as nc

import microhh_tools as mht


"""
Settings
"""
float_type = np.float32     # np.float32 for -USESP=true, else np.float64.

dust_bins = np.array([0, 2, 10, 20, 58, 83, 440])  # (μm)

xsize = 12800
ysize = 6400
zsize = 3200

itot = 256
jtot = 128
ktot = 64

endtime = 10800

dx = xsize / itot
dy = ysize / jtot
dz = zsize / ktot


"""
Define vertical grid and input profiles.
"""
x = np.arange(dx/2, xsize, dx)
y = np.arange(dy/2, ysize, dy)
z = np.arange(dz/2, zsize, dz)

u  = np.zeros_like(z)
th = np.zeros_like(z)

# Linearly statified profile.
th = 290 + 0.003 * z
u = np.ones(ktot) * 2

# All scalars start at concentration zero, with
# inflow of air with concentration zero at lateral boundaries.
dust_list = [f'{dust_bins[i]}-{dust_bins[i+1]}um' for i in range(dust_bins.size - 1)]

scalars = {}
for scalar in dust_list:
    scalars[scalar] = np.zeros(ktot)


"""
Define surface emissions dust scalars + gravitational settling velocities.
Following: DOI: 10.1063/1.5022089
"""
rho_p = 1500   # Density particles [kg m-3]
rho_a = 1.225  # Reference density air [kg m-3]
nu = 1e-5      # Kinematic viscosity air [m2 s-1]
g = 9.81       # Gravitational acceleration [m s-2]

dust_diameter = 0.5*(dust_bins[1:] + dust_bins[:-1]) * 1e-6

tau_p = dust_diameter**2 * rho_p / (18 * nu * rho_a)
w_s = -tau_p * g

# Create circular field with dust emissions.
x0 = 0.15 * xsize
y0 = 0.5 * ysize
r = 750

Y, X = np.meshgrid(y, x, indexing='ij')
D = np.sqrt((X - x0)**2 + (Y - y0)**2)
field_mask = D < r

field_flux = np.zeros((jtot, itot), dtype=float_type)
field_flux[field_mask] = 1.
for scalar in dust_list:
    field_flux.tofile('{}_bot_in.0000000'.format(scalar))


"""
Set/write new namelist.
"""
ini = mht.Read_namelist('dust.ini.base')

ini['grid']['itot'] = itot
ini['grid']['jtot'] = jtot
ini['grid']['ktot'] = ktot

ini['grid']['xsize'] = xsize
ini['grid']['ysize'] = ysize
ini['grid']['zsize'] = zsize

ini['buffer']['zstart'] = 0.75*zsize

ini['time']['endtime'] = endtime

ini['fields']['slist'] = dust_list
ini['advec']['fluxlimit_list'] = dust_list
ini['limiter']['limitlist'] = dust_list
ini['boundary']['scalar_outflow'] = dust_list
ini['boundary']['sbot_2d_list'] = dust_list

ini['dust']['dustlist'] = dust_list
for i in range(len(dust_list)):
    ini['dust'][f'ws[{dust_list[i]}]'] = w_s[i]

# Statistics/crosses/...
scalar_crosses = dust_list + [s+'_path' for s in dust_list]
ini['cross']['crosslist'] = scalar_crosses + ['th', 'u', 'v', 'w']
ini['cross']['xz'] = y0
ini['cross']['yz'] = x0

ini.save('dust.ini', allow_overwrite=True)


"""
Create input NetCDF file.
"""
def add_var(name, dims, values, nc_group):
    nc_var = nc_group.createVariable(name, float_type, dims)
    nc_var[:] = values

nc_file = nc.Dataset('dust_input.nc', mode='w', datamodel='NETCDF4')
nc_file.createDimension('z', ktot)
add_var('z',  ('z'), z,  nc_file)

nc_init = nc_file.createGroup('init');
add_var('th', ('z'), th, nc_init)
add_var('u',  ('z'), u,  nc_init)

# Same profile with concentration = 0 for all scalars:
for scalar, prof in scalars.items():
    add_var(scalar, ('z'), prof, nc_init)
    add_var('{}_inflow'.format(scalar), ('z'), prof, nc_init)

nc_file.close()
