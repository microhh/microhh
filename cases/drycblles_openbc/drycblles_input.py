import numpy as np
import netCDF4 as nc

from lbc_input import LBC_input

float_type = "f8"
# float_type = "f4"

# Get number of vertical levels and size from .ini file
with open('drycblles.ini') as f:
    for line in f:
        if(line.split('=')[0]=='itot'):
            itot = int(line.split('=')[1])
        if(line.split('=')[0]=='jtot'):
            jtot = int(line.split('=')[1])
        if(line.split('=')[0]=='ktot'):
            ktot = int(line.split('=')[1])
        if(line.split('=')[0]=='xsize'):
            xsize = float(line.split('=')[1])
        if(line.split('=')[0]=='ysize'):
            ysize = float(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / ktot

dthetadz = 0.003

# set the height
z  = np.linspace(0.5*dz, zsize-0.5*dz, ktot)
u  = np.zeros(np.size(z))
v  = np.zeros(np.size(z))
th = np.zeros(np.size(z))

# linearly stratified profile
for k in range(ktot):
    th  [k] = 300. + dthetadz*z[k]

"""
# well mixed profile with jump
h    = 1000.
dth  = 10.
dthz = 100.

for k in range(ktot):
    if(z[k] <= h - 0.5*dthz):
        th[k] = 300.
    elif(z[k] <= h + 0.5*dthz):
        th[k] = 300. + dth/dthz * (z[k]-(h-0.5*dthz))
    else:
        th[k] = 300. + dth + dthetadz*(z[k]-(h+0.5*dthz))
    thls[k] = 2.*(z[k]/zsize - 0.5) / 3600.
    wls [k] = -0.01*(z[k]/zsize)
"""

nc_file = nc.Dataset("drycblles_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", ktot)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_u  = nc_group_init.createVariable("u" , float_type, ("z"))
nc_v  = nc_group_init.createVariable("v" , float_type, ("z"))
nc_th = nc_group_init.createVariable("th", float_type, ("z"))

nc_z [:] = z [:]
nc_u [:] = u [:]
nc_v [:] = v [:]
nc_th[:] = th[:]

nc_file.close()

"""
Create lateral boundaries
"""
dx = xsize / itot
dy = ysize / jtot

x = np.arange(dx/2, xsize, dx)
xh = np.arange(0, xsize, dx)
y = np.arange(dy/2, ysize, dy)
yh = np.arange(0, ysize, dy)

time = np.array([0, 10800])
fields = ['u','v','s','th']

lbc = LBC_input(fields, itot, jtot, ktot, time)

lbc.s_west[:] = np.cos(4*np.pi*y/ysize)
lbc.s_south[:] = np.cos(4*np.pi*x/xsize)
#lbc.s_north[:] = np.cos(4*np.pi*y/ysize)
#lbc.s_east[:] = np.cos(4*np.pi*x/xsize)

lbc.th_west[:] = th[None, :, None]
lbc.th_south[:] = th[None, :, None]
lbc.th_north[:] = th[None, :, None]
lbc.th_east[:] = th[None, :, None]

u_west = 2.1
u_east = 0.

v_south = 0.
v_north = 2.

rnd_amp = 0.01

lbc.u_west[:, :, :]  = u_west + rnd_amp * np.random.rand(2, ktot, jtot) * (zsize - z[None, :, None])/zsize
lbc.u_east[:, :, :]  = u_east + rnd_amp * np.random.rand(2, ktot, jtot) * (zsize - z[None, :, None])/zsize
lbc.u_south[:, :, :] = u_west + (u_east - u_west) * xh[None, None, :]/xsize + rnd_amp * np.random.rand(2, ktot, itot) * (zsize - z[None, :, None])/zsize
lbc.u_north[:, :, :] = u_west + (u_east - u_west) * xh[None, None, :]/xsize + rnd_amp * np.random.rand(2, ktot, itot) * (zsize - z[None, :, None])/zsize

lbc.v_west[:, :, :]  = v_south + (v_north - v_south) * yh[None, None, :]/ysize + rnd_amp * np.random.rand(2, ktot, itot) * (zsize - z[None, :, None])/zsize
lbc.v_east[:, :, :]  = v_south + (v_north - v_south) * yh[None, None, :]/ysize + rnd_amp * np.random.rand(2, ktot, itot) * (zsize - z[None, :, None])/zsize
lbc.v_south[:, :, :] = v_south + rnd_amp * np.random.rand(2, ktot, jtot) * (zsize - z[None, :, None])/zsize
lbc.v_north[:, :, :] = v_north + rnd_amp * np.random.rand(2, ktot, jtot) * (zsize - z[None, :, None])/zsize

lbc.to_netcdf('drycblles')
