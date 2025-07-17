import netCDF4 as nc
import numpy as np
import sys

import microhhpy.io as io


# Case settings.
domain = sys.argv[1]
dtype = np.float64


"""
Time.
"""
endtime = 3600
lbc_freq = 30


"""
Grid & nesting settings.
"""
# Grid settings outer domain.
itot = 64
jtot = 64
ktot = 64

xsize = 3200
ysize = 3200
zsize = 3200

dx = xsize / itot
dy = ysize / jtot
dz = zsize / ktot

# Nest settings.
xstart_sub = 800
ystart_sub = 800

xend_sub = 2400
yend_sub = 2400

xsize_sub = xend_sub - xstart_sub
ysize_sub = yend_sub - ystart_sub

grid_ratio = 2

itot_sub = int(xsize_sub / dx * grid_ratio + 0.5)
jtot_sub = int(ysize_sub / dy * grid_ratio + 0.5)

dx_sub = dx / grid_ratio
dy_sub = dy / grid_ratio

# Number of lateral buffer/sponge points.
n_sponge = 3
n_ghost = 3


"""
Define initial fields/profiles.
"""
dthetadz = 0.003

z  = np.arange(0.5*dz, zsize, dz)
zh = np.arange(0, zsize, dz)

u  = np.zeros(np.size(z)) + 2
v  = np.zeros(np.size(z)) + 1
th = 300. + dthetadz * z


"""
Write case_input.nc
"""
nc_file = nc.Dataset("drycblles_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", ktot)
nc_z  = nc_file.createVariable("z" , dtype, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_u  = nc_group_init.createVariable("u" , dtype, ("z"))
nc_v  = nc_group_init.createVariable("v" , dtype, ("z"))
nc_th = nc_group_init.createVariable("th", dtype, ("z"))

nc_z [:] = z [:]
nc_u [:] = u [:] + 1
nc_v [:] = v [:] 
nc_th[:] = th[:]

nc_file.close()


"""
Update .ini
"""
ini = io.read_ini('drycblles.ini.base')

ini['time']['endtime'] = endtime

if domain == 'outer':

    ini['grid']['itot'] = itot
    ini['grid']['jtot'] = jtot
    ini['grid']['ktot'] = ktot

    ini['grid']['xsize'] = xsize
    ini['grid']['ysize'] = ysize
    ini['grid']['zsize'] = zsize

    ini['cross']['sampletime'] = lbc_freq
    ini['cross']['xz'] = ysize / 2.
    ini['cross']['yz'] = xsize / 2.

    ini['pres']['sw_openbc'] = False
    ini['boundary_lateral']['sw_openbc'] = False

    ini['subdomain']['sw_subdomain'] = True
    ini['subdomain']['xstart'] = xstart_sub
    ini['subdomain']['xend'] = xend_sub
    ini['subdomain']['ystart'] = ystart_sub
    ini['subdomain']['yend'] = yend_sub
    ini['subdomain']['grid_ratio'] = grid_ratio
    ini['subdomain']['n_ghost'] = n_ghost
    ini['subdomain']['n_sponge'] = n_sponge
    ini['subdomain']['savetime'] = lbc_freq

elif domain == 'inner':

    ini['grid']['itot'] = itot_sub
    ini['grid']['jtot'] = jtot_sub
    ini['grid']['ktot'] = ktot

    ini['grid']['xsize'] = xsize_sub
    ini['grid']['ysize'] = ysize_sub
    ini['grid']['zsize'] = zsize

    ini['cross']['sampletime'] = lbc_freq
    ini['cross']['yz'] = xsize_sub / 2.
    ini['cross']['xz'] = ysize_sub / 2.

    ini['pres']['sw_openbc'] = True
    ini['boundary_lateral']['sw_openbc'] = True
    ini['boundary_lateral']['loadfreq'] = lbc_freq
    ini['boundary_lateral']['n_sponge'] = n_sponge

    ini['subdomain']['sw_subdomain'] = False

io.save_ini(ini, 'drycblles.ini')
