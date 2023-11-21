import numpy as np
import netCDF4 as nc

from lbc_input import lbc_input

#float_type = "f8"
float_type = "f4"

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
z  = np.arange(0.5*dz, zsize, dz)
zh = np.arange(0, zsize, dz)

u  = np.zeros(np.size(z))
v  = np.zeros(np.size(z))
th = np.zeros(np.size(z))

# linearly stratified profile
for k in range(ktot):
    th[k] = 300. + dthetadz*z[k]

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

time = np.array([0])
fields = ['th','s', 'u', 'v']
nghost = 3

lbc = lbc_input(fields, x, y, z, xh, yh, zh, time, nghost)

ii = 1

lbc.th_west[:,:,:,:]  = th[None,:,None,None]
lbc.th_east[:,:,:,:]  = th[None,:,None,None]
lbc.th_south[:,:,:,:] = th[None,:,None,None]
lbc.th_north[:,:,:,:] = th[None,:,None,None]

lbc.s_west[:,:,:,:]  = 0
lbc.s_east[:,:,:,:]  = 0
lbc.s_south[:,:,:,:] = 0
lbc.s_north[:,:,:,:] = 0
lbc.s_west[:,:15,:,:] = 1

lbc.u_west[:,:,:,:]  = 0.1
lbc.u_east[:,:,:,:]  = 0
lbc.u_south[:,:,:,:] = 0
lbc.u_north[:,:,:,:] = 0

lbc.v_west[:,:,:,:]  = 0
lbc.v_east[:,:,:,:]  = 0
lbc.v_south[:,:,:,:] = 0
lbc.v_north[:,:,:,:] = 0.1

#for fld in fields:
#    for loc in ['west', 'east', 'south', 'north']:
#
#        lbc_ref = lbc[f'{fld}_{loc}']
#        dims = lbc_ref.shape
#
#        for t in range(dims[0]):
#            for k in range(dims[1]):
#                for j in range(dims[2]):
#                    for i in range(dims[3]):
#                        lbc_ref[t,k,j,i] = t*1000 + k*100 + j*10 + i

lbc.to_netcdf('drycblles_lbc_input.nc')

#lbc.s_west [0,:] = 0
#lbc.s_south[0,:] = 0
#lbc.s_north[0,:] = 0
#lbc.s_east [0,:] = 0
#
#lbc.s_west [1,:] = 10
#lbc.s_south[1,:] = 0
#lbc.s_north[1,:] = 0
#lbc.s_east [1,:] = 0

#u_west = 2.1
#u_east = 0.
#v_south = 0.
#v_north = 2.

#u_west = 0
#u_east = 0.
#v_south = 0.
#v_north = 0.

#rnd_amp = 0 #0.01
#
#def make_rand(n0, n1, n2):
#    rnd = np.random.rand(n0, n1, n2)
#    return rnd - rnd.mean()
#
#lbc.u_west[:, :, :]  = u_west + rnd_amp * make_rand(2, ktot, jtot) * (zsize - z[None, :, None])/zsize
#lbc.u_east[:, :, :]  = u_east + rnd_amp * make_rand(2, ktot, jtot) * (zsize - z[None, :, None])/zsize
#lbc.u_south[:, :, :] = u_west + (u_east - u_west) * xh[None, None, :]/xsize + rnd_amp * make_rand(2, ktot, itot) * (zsize - z[None, :, None])/zsize
#lbc.u_north[:, :, :] = u_west + (u_east - u_west) * xh[None, None, :]/xsize + rnd_amp * make_rand(2, ktot, itot) * (zsize - z[None, :, None])/zsize
#
#lbc.v_west[:, :, :]  = v_south + (v_north - v_south) * yh[None, None, :]/ysize + rnd_amp * make_rand(2, ktot, jtot) * (zsize - z[None, :, None])/zsize
#lbc.v_east[:, :, :]  = v_south + (v_north - v_south) * yh[None, None, :]/ysize + rnd_amp * make_rand(2, ktot, jtot) * (zsize - z[None, :, None])/zsize
#lbc.v_south[:, :, :] = v_south + rnd_amp * make_rand(2, ktot, itot) * (zsize - z[None, :, None])/zsize
#lbc.v_north[:, :, :] = v_north + rnd_amp * make_rand(2, ktot, itot) * (zsize - z[None, :, None])/zsize

