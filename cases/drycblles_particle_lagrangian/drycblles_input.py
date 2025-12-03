import numpy as np
import netCDF4 as nc
import h5py

# Available in microhh/python.
import microhh_tools as mht


"""
Settings
"""
float_type = np.float32     # np.float32 for -USESP=true, else np.float64.

xsize = 3200
ysize = 1600

itot = 128
jtot = 64


"""
Define vertical grid.
"""
def grid(dz0, alpha, ktot):
    dz = dz0 * alpha**np.arange(ktot)
    zh = np.zeros(ktot+1)
    zh[1:] = np.cumsum(dz)
    z = 0.5 * (zh[:-1] + zh[1:])
    return z, zh[-1]

#z, zsize = grid(20, 1.03, 64)
z, zsize = grid(10, 1.013, 128)
ktot = z.size


"""
Initial profiles.
"""
th = 290 + 0.003 * z
u = np.zeros(ktot) 


"""
Define initial location particles.
"""
n_particles = 10000

particle_id = np.arange(n_particles, dtype=np.int32)

"""
x0 = xsize / 2
y0 = ysize / 2
z0 = 20
size = 0.5

#xp  = np.random.uniform(x0-size, x0+size, n_particles).astype(float_type)
xp  = np.linspace(x0-size, x0+size, n_particles).astype(float_type)
yp  = np.random.uniform(y0-size, y0+size, n_particles).astype(float_type)
zp  = np.random.uniform(z0-size, z0+size, n_particles).astype(float_type)
"""

xp  = np.random.uniform(0, xsize, n_particles).astype(float_type)
yp  = np.random.uniform(0, ysize, n_particles).astype(float_type)
zp  = np.random.uniform(0, zsize, n_particles).astype(float_type)

# If this fails; try downgrading both NetCDF4 and H5PY.
# See this issue: https://github.com/h5py/h5py/issues/2453
# pip uninstall numpy netCDF4 h5py
# pip install numpy==1.26.4
# pip install h5py==3.10.0
# pip install netCDF4==1.6.5
with h5py.File('particles.0000000.h5', 'w') as f:
    f.create_dataset('particle_id', data=particle_id)
    f['particle_id'].make_scale('particles')
    
    f.create_dataset('x', data=xp)
    f['x'].dims[0].attach_scale(f['particle_id'])
    
    f.create_dataset('y', data=yp)
    f['y'].dims[0].attach_scale(f['particle_id'])
    
    f.create_dataset('z', data=zp)
    f['z'].dims[0].attach_scale(f['particle_id'])


"""
Set/write new namelist.
"""
ini = mht.Read_namelist('drycblles.ini.base')

ini['grid']['itot'] = itot
ini['grid']['jtot'] = jtot
ini['grid']['ktot'] = ktot

ini['grid']['xsize'] = xsize
ini['grid']['ysize'] = ysize
ini['grid']['zsize'] = zsize

ini['buffer']['zstart'] = 0.75*zsize

ini['cross']['xz'] = 0.5*ysize

ini.save('drycblles.ini', allow_overwrite=True)


"""
Create input NetCDF file.
"""
def add_var(name, dims, values, nc_group, dtype=float_type):
    nc_var = nc_group.createVariable(name, dtype, dims)
    nc_var[:] = values


nc_file = nc.Dataset('drycblles_input.nc', mode='w', datamodel='NETCDF4')
nc_file.createDimension('z', ktot)
add_var('z',  ('z'), z,  nc_file)

nc_init = nc_file.createGroup('init');
add_var('th', ('z'), th, nc_init)
add_var('u',  ('z'), u,  nc_init)

nc_file.close()
