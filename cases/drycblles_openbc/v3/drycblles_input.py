import matplotlib.pyplot as pl
import numpy as np
import netCDF4 as nc
import xarray as xr
import sys

from lbc_input import lbc_input
import microhh_tools as mht

pl.close('all')

domain = sys.argv[1]

#float_type = "f8"
float_type = "f4"

swadvec = '2'   # needed for correct # gcs.

ini = mht.Read_namelist('drycblles.ini.base')

"""
Grid & nesting settings.
"""

# Grid settings outer domain.
itot = 128
jtot = 128
ktot = 128

xsize = 6400
ysize = 6400
zsize = 3200

dx = xsize / itot
dy = ysize / jtot
dz = zsize / ktot

# Nest settings.
refinement_fac = 3

i0_nest = 32
j0_nest = 32

# Size domain in parent coordinates!
# The nest itself has `refinement_fac` times as many grid points.
itot_nest = 64
jtot_nest = 64

# Number of lateral buffer points.
nbuffer = 5

# Number of ghost cells.
if swadvec == '2':
    nghost = 1
elif swadvec == '2i4':
    nghost = 2
elif swadvec == '2i5':
    nghost = 3

xstart_nest = i0_nest*dx
ystart_nest = j0_nest*dy

xend_nest = xstart_nest + itot_nest * dx
yend_nest = ystart_nest + jtot_nest * dy

xsize_nest = xend_nest - xstart_nest
ysize_nest = yend_nest - ystart_nest

dx_nest = xsize_nest / (itot_nest * refinement_fac)
dy_nest = ysize_nest / (jtot_nest * refinement_fac)

"""
Define initial fields/profiles.
"""
dthetadz = 0.003

z  = np.arange(0.5*dz, zsize, dz)
zh = np.arange(0, zsize, dz)

u  = np.zeros(np.size(z))
v  = np.zeros(np.size(z))
th = np.zeros(np.size(z))

for k in range(ktot):
    th[k] = 300. + dthetadz*z[k]

"""
Write case_input.nc
"""
nc_file = nc.Dataset("drycblles_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", ktot)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_u  = nc_group_init.createVariable("u" , float_type, ("z"))
nc_v  = nc_group_init.createVariable("v" , float_type, ("z"))
nc_th = nc_group_init.createVariable("th", float_type, ("z"))

nc_z [:] = z [:]
nc_u [:] = u [:] + 2
nc_v [:] = v [:] + 1
nc_th[:] = th[:]

nc_file.close()

"""
Update .ini file.
"""
ini['advec']['swadvec'] = swadvec

if domain == 'outer':

    x0 = xstart_nest - nghost * dx_nest
    x1 = xstart_nest + nbuffer * dx_nest
    yzw = np.arange(x0, x1+1e-3, dx)

    x0 = xend_nest - nbuffer * dx_nest
    x1 = xend_nest + nghost * dx_nest
    yze = np.arange(x0, x1+1e-3, dx)

    yz = np.concatenate((yzw, yze))

    y0 = ystart_nest - nghost * dy_nest
    y1 = ystart_nest + nbuffer * dy_nest
    xzs = np.arange(y0, y1+1e-3, dy)

    y0 = yend_nest - nbuffer * dy_nest
    y1 = yend_nest + nghost * dy_nest
    xzn = np.arange(y0, y1+1e-3, dy)

    xz = np.concatenate((xzs, xzn))

#    #yz = np.array([
#    #    xstart_nest-2.5*dx, xstart_nest-1.5*dx, xstart_nest-0.5*dx, xstart_nest+0.5*dx,
#    #    xend_nest-0.5*dx, xend_nest+0.5*dx, xend_nest+1.5*dx, xend_nest+2.5*dx])
#    #xz = np.array([
#    #    ystart_nest-2.5*dy, ystart_nest-1.5*dy, ystart_nest-0.5*dy, ystart_nest+0.5*dy,
#    #    yend_nest-0.5*dy, yend_nest+0.5*dy, yend_nest+1.5*dy, yend_nest+2.5*dy])

    ini['grid']['itot'] = itot
    ini['grid']['jtot'] = jtot
    ini['grid']['ktot'] = ktot

    ini['grid']['xsize'] = xsize
    ini['grid']['ysize'] = ysize
    ini['grid']['zsize'] = zsize

    ini['cross']['yz'] = list(yz)
    ini['cross']['xz'] = list(xz)

    ini['pres']['swopenbc']=False
    ini['boundary']['sw_inoutflow']=False

elif domain == 'inner':

    ini['grid']['itot'] = itot_nest * refinement_fac
    ini['grid']['jtot'] = jtot_nest * refinement_fac
    ini['grid']['ktot'] = ktot

    ini['grid']['xsize'] = xsize_nest
    ini['grid']['ysize'] = ysize_nest
    ini['grid']['zsize'] = zsize

    ini['cross']['yz'] = (xend_nest+xstart_nest)/2
    ini['cross']['xz'] = (yend_nest+ystart_nest)/2

    ini['pres']['swopenbc']=True
    ini['boundary']['sw_inoutflow']=True

ini.save('drycblles.ini', allow_overwrite=True)


"""
Create lateral boundaries
"""
if domain == 'inner':

    fields = ['u', 'v', 'w', 'th', 's']

    # Read cross-sections.
    xz = {}
    yz = {}
    for fld in fields:
        xz[fld] = xr.open_dataset(f'outer/{fld}.xz.nc', decode_times=False)
        yz[fld] = xr.open_dataset(f'outer/{fld}.yz.nc', decode_times=False)
    time = xz[list(xz.keys())[0]].time.values

    # Define nest grid.
    x = np.arange(dx_nest/2, xsize_nest, dx_nest)
    xh = np.arange(0, xsize_nest, dx_nest)
    y = np.arange(dy_nest/2, ysize_nest, dy_nest)
    yh = np.arange(0, ysize_nest, dy_nest)

    lbc = lbc_input(fields, x, y, z, xh, yh, zh, time, nghost, nbuffer)

    # Add offsets for easier interpolation.
    for v in lbc.variables:
        if 'x' in v:
            lbc[v] = lbc[v] + xstart_nest
        if 'y' in v:
            lbc[v] = lbc[v] + ystart_nest

    # West boundaries
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
            ip = cc[fld].interp({yloc_in: lbc[yloc], xloc_in: lbc[xloc]})

            # Check if interpolation was success.
            if np.any(np.isnan(ip[fld].values)):
                raise Exception('Interpolated BCs contain NaNs!')

            lbc_in[:] = ip[fld].values

    lbc.to_netcdf('drycblles_lbc_input.nc')
