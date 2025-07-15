import matplotlib.pyplot as pl
import numpy as np
import netCDF4 as nc
import xarray as xr
from numba import jit, prange
import sys
import glob
from datetime import datetime
import asyncio

# Avaiable in `microhh/python`:
import microhh_lbc_tools as mlt
import microhh_tools as mht

# Custom tools.
from lbc_input import lbc_input

pl.close('all')

domain = sys.argv[1]
dtype = np.float64
swadvec = '2'   # needed for correct # gcs.

ini = mht.Read_namelist('drycblles.ini.base')


"""
Time.
"""
endtime = 300
lbc_freq = 30


"""
Grid & nesting settings.
"""
# Grid settings outer domain.
itot = 32
jtot = 32
ktot = 32

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
n_sponge = 5
n_ghost = 1


"""
Define initial fields/profiles.
"""
dthetadz = 0.003

z  = np.arange(0.5*dz, zsize, dz)
zh = np.arange(0, zsize, dz)

u  = np.zeros(np.size(z)) + 1
v  = np.zeros(np.size(z))
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
ini['time']['endtime'] = endtime


if domain == 'outer':

    xz, yz = mlt.get_cross_locations_for_lbcs(
            xstart_sub, ystart_sub,
            xend_sub, yend_sub,
            dx, dy,
            dx_sub, dy_sub,
            n_ghost, n_sponge)

    ini['grid']['itot'] = itot
    ini['grid']['jtot'] = jtot
    ini['grid']['ktot'] = ktot

    ini['grid']['xsize'] = xsize
    ini['grid']['ysize'] = ysize
    ini['grid']['zsize'] = zsize

    ini['cross']['sampletime'] = lbc_freq
    ini['cross']['xz'] = list(xz)
    ini['cross']['yz'] = list(yz)

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
    ini['cross']['yz'] = (xend_sub + xstart_sub)/2
    ini['cross']['xz'] = (yend_sub + ystart_sub)/2

    ini['pres']['sw_openbc'] = True
    ini['boundary_lateral']['sw_openbc'] = True
    ini['boundary_lateral']['loadfreq'] = lbc_freq

    ini['subdomain']['sw_subdomain'] = False

ini.save('drycblles.ini', allow_overwrite=True)


"""
Create lateral boundaries
"""
class Timer:
    def __init__(self):
        self.tstart = datetime.now()
        self.elapsed = 0

    def pause(self):
        self.elapsed += (datetime.now() - self.tstart).total_seconds()

    def restart(self):
        self.tstart = datetime.now()

    def stop(self):
        self.elapsed += (datetime.now() - self.tstart).total_seconds()
        print(f'Elapsed: {self.elapsed}')


if domain == 'inner':

    fields = ['u', 'v', 'w', 'th', 's']

    t = Timer()

    # Read cross-sections.
    xz = {}
    yz = {}
    for fld in fields:
        xz[fld] = xr.open_dataset(f'outer/{fld}.xz.nc', decode_times=False)
        yz[fld] = xr.open_dataset(f'outer/{fld}.yz.nc', decode_times=False)
    time = xz[list(xz.keys())[0]].time.values

    t.pause()

    # Define nest grid.
    x = np.arange(dx_sub/2, xsize_sub, dx_sub)
    xh = np.arange(0, xsize_sub, dx_sub)
    y = np.arange(dy_sub/2, ysize_sub, dy_sub)
    yh = np.arange(0, ysize_sub, dy_sub)

    lbc = lbc_input(fields, x, y, z, xh, yh, zh, time, n_ghost, n_sponge, x_offset=xstart_sub, y_offset=ystart_sub)

    t.restart()

    for loc in ['west', 'east', 'north', 'south']:
        for fld in fields:
            # Short cuts.
            lbc_in = lbc[f'{fld}_{loc}']
            dims = lbc_in.dims

            # Dimensions in LBC file.
            xloc, yloc = dims[3], dims[2]

            # Dimensions in cross-section.
            xloc_in = 'xh' if 'xh' in xloc else 'x'
            yloc_in = 'yh' if 'yh' in yloc else 'y'

            # Switch between yz and xz crosses.
            cc = yz if loc in ['west','east'] else xz

            # Interpolate!
            ip = cc[fld].interp({yloc_in: lbc[yloc], xloc_in: lbc[xloc]}, method='nearest')

            # Check if interpolation was success.
            if np.any(np.isnan(ip[fld].values)):
                raise Exception('Interpolated BCs contain NaNs!')

            lbc_in[:] = ip[fld].values

    t.stop()

    #lbc.to_netcdf('lbc_old.nc')

    # Save as binary input files.
    for fld in fields:
        for loc in ['west', 'east', 'north', 'south']:
            for t, time in enumerate(lbc.time):
                lbc_in = lbc[f'{fld}_{loc}'][t]
                lbc_in.values.astype(dtype).tofile(f'lbc_{fld}_{loc}.{int(time):07d}')
