import numpy as np
import netCDF4 as nc

import microhh_tools as mht
from emission_input import Emission_input

# `microhhpy` is for now available at https://github.com/bartvstratum/microhhpy
# See installation notes in the README.md file.
# In the future, this will be uploaded to PyPI.
from microhhpy.thermo import Basestate_dry


"""
Settings
"""
TF = np.float32
scalars = ['s1', 's2']       # All scalars.
scalars_outflow = ['s1']     # Scalars with non-periodic lateral BCs.
swtimedep_emission = False    # Time dependent emissions.

xsize = 6400
ysize = 4800
zsize = 3200

itot = 128
jtot = 96
ktot = 64

endtime = 7200

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

# Emission strength.
if swtimedep_emission:
    strengths = np.array([0, 1, 2])
    times = np.array([0, 3600, 7200])
else:
    strengths = np.array([1])
    times = np.array([0]) 

"""
Define horizontal and vertical grid and input profiles.
"""
dx = xsize / itot
dy = ysize / jtot
dz = zsize / ktot

x  = np.arange(0.5*dx, xsize, dx)
y  = np.arange(0.5*dy, ysize, dy)
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
s_lbc = np.zeros_like(z)


"""
3D emission input for more complex scenarios like cities.
"""
bs = Basestate_dry(th, 1e5, z, zsize, remove_ghost=True, dtype=TF)

# Create emission input instance.
emiss = Emission_input(scalars, times, x, y, z, np.full(ktot, dz), bs.rho, TF=TF)

# Add sources.
for time, strength in zip(times, strengths):
    emiss.add_gaussian('s1', strength, time, source_x0, source_y0, source_z0, sigma_x, sigma_y, sigma_z)
    emiss.add_gaussian('s2', strength, time, source_x0, source_y0, source_z0, sigma_x, sigma_y, sigma_z)

    #emiss.add_point('s1', strength, time, source_x0, source_y0, source_z0)
    #emiss.add_point('s2', strength, time, source_x0, source_y0, source_z0)

# Clip to required vertical extent.
emiss.clip()

# Save emissions as binary input files for MicroHH.
emiss.to_binary(path='.')


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
ini['boundary']['swtimedep_outflow'] = False

# Sources
def const_list(value):
    value = int(value) if isinstance(value, bool) else value
    return len(scalars)*[value]

ini['source']['swsource'] = '3d'
ini['source']['sourcelist'] = scalars
ini['source']['ktot'] = emiss.kmax
ini['source']['swtimedep'] = swtimedep_emission
ini['source']['loadtime'] = 3600

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
    nc_var = nc_group.createVariable(name, TF, dims)
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

# Same in/outflow concentration = 0 for all scalars:
for scalar in scalars_outflow:
    add_var('{}_inflow'.format(scalar), ('z'), s_lbc, nc_init)

nc_file.close()
