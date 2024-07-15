import numpy as np
import netCDF4 as nc

float_type = "f8"
# float_type = "f4"

# Get number of vertical levels and size from .ini file
with open('drycblles.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            itot = int(line.split('=')[1])
        if(line.split('=')[0]=='jtot'):
            jtot = int(line.split('=')[1])
        if(line.split('=')[0]=='ktot'):
            ktot = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / ktot

dthetadz = 0.003

# set the height
z  = np.arange(dz/2, zsize, dz)
zh = np.arange(0, zsize, dz)
u  = np.zeros(np.size(z))+1
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

## 3D buffer.
#zstart = 2400
#kstart  = int(np.where(z  >= zstart)[0][0])
#kstarth = int(np.where(zh >= zstart)[0][0])
#
#ksize = ktot - kstart
#ksizeh = ktot - kstarth
#
#u_buf  = np.zeros((ksize,  jtot, itot), float_type)
#v_buf  = np.zeros((ksize,  jtot, itot), float_type)
#w_buf  = np.zeros((ksizeh, jtot, itot), float_type)
#th_buf = np.zeros((ksizeh, jtot, itot), float_type)
#
#u_buf[:,:,:] = u[kstart:,None,None]
#v_buf[:,:,:] = v[kstart:,None,None]
#th_buf[:,:,:] = th[kstart:,None,None]
#
#u_buf.tofile('u_buffer.0000000')
#v_buf.tofile('v_buffer.0000000')
#w_buf.tofile('w_buffer.0000000')
#th_buf.tofile('th_buffer.0000000')
#
#th_buf[:] += 1.
#
#u_buf.tofile('u_buffer.0003600')
#v_buf.tofile('v_buffer.0003600')
#w_buf.tofile('w_buffer.0003600')
#th_buf.tofile('th_buffer.0003600')
#
#u_buf.tofile('u_buffer.0007200')
#v_buf.tofile('v_buffer.0007200')
#w_buf.tofile('w_buffer.0007200')
#th_buf.tofile('th_buffer.0007200')
#
#u_buf.tofile('u_buffer.0010800')
#v_buf.tofile('v_buffer.0010800')
#w_buf.tofile('w_buffer.0010800')
#th_buf.tofile('th_buffer.0010800')
