import matplotlib.pyplot as pl
import numpy as np
import netCDF4 as nc

import microhh_tools as mht
import microhh_lbc_tools as mlt

TF = np.float64

ini = mht.Read_namelist('drycbl.ini')

dx = ini['grid']['xsize'] / ini['grid']['itot']
dy = ini['grid']['ysize'] / ini['grid']['jtot']
dz = ini['grid']['zsize'] / ini['grid']['ktot']

x = np.arange(dx/2, ini['grid']['xsize'], dx)
y = np.arange(dy/2, ini['grid']['ysize'], dy)
z = np.arange(dz/2, ini['grid']['zsize'], dz)
zh = np.arange(0, ini['grid']['zsize'], dz)

"""
Thermodyamic structure
"""
thl = 290 + 0.006 * z
qt = np.zeros_like(thl)

"""
Spatial `p_top` field.
"""
fc = 0.0001

ug_top = 0
vg_top = 5

dpdx = 0.01 # vg_top * fc # * rho?
dpdy = 0.01 #-ug_top * fc # * rho?

print(f'dpdx={dpdx*1000} Pa/km, dpdy={dpdy*1000} Pa/km.')

p_bottom_left = 70_000

xx,yy = np.meshgrid(x,y)
p_top = p_bottom_left + xx * dpdx + yy * dpdy

print(f'p_top_min={p_top.min()} Pa, p_top_max={p_top.max()} Pa')

p_top.astype(TF).tofile('phydro_tod.0000000')


"""
Case input.
"""
u = np.ones_like(thl) * ug_top
v = np.ones_like(thl) * vg_top

nc_file = nc.Dataset('drycbl_input.nc', mode='w', datamodel='NETCDF4', clobber=True)

nc_file.createDimension('z', ini['grid']['ktot'])
nc_z  = nc_file.createVariable('z' , TF, ('z'))

nc_group_init = nc_file.createGroup('init');

nc_thl = nc_group_init.createVariable('thl', TF, ('z'))
nc_qt  = nc_group_init.createVariable('qt', TF, ('z'))

nc_u = nc_group_init.createVariable('u', TF, ('z'))
nc_v = nc_group_init.createVariable('v', TF, ('z'))

nc_z[:] = z[:]
nc_thl[:] = thl[:]
nc_qt[:] = qt[:]
nc_u[:] = u[:]
nc_v[:] = v[:]

nc_file.close()


"""
Lateral boundaries.
"""
fields = ['thl', 'qt', 'u', 'v', 'w']
time = np.array([0])

lbc_ds = mlt.get_lbc_xr_dataset(
        fields,
        ini['grid']['xsize'],
        ini['grid']['ysize'],
        ini['grid']['itot'],
        ini['grid']['jtot'],
        z, zh, time,
        n_ghost=1,
        n_sponge=1,
        x_offset=0,
        y_offset=0,
        dtype=TF)

for edge in ['west', 'north', 'east', 'south']:
    lbc_ds[f'thl_{edge}'][:,:,:,:] = thl[None,:,None,None]
    lbc_ds[f'qt_{edge}' ][:,:,:,:] = qt [None,:,None,None]
    lbc_ds[f'u_{edge}'  ][:,:,:,:] = u  [None,:,None,None]
    lbc_ds[f'v_{edge}'  ][:,:,:,:] = v  [None,:,None,None]

mlt.write_lbcs_as_binaries(lbc_ds, dtype=TF)