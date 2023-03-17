import numpy as np

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

dx, dy, dz = xsize/itot, ysize/jtot, zsize/ktot

xh = np.arange(0, xsize - dx/2, dx)
x = np.arange(dx/2, xsize , dx)
yh = np.arange(0, ysize - dy/2, dy)
y = np.arange(dy/2, ysize , dy)
zh = np.arange(0, zsize - dz/2, dz)
z = np.arange(dz/2, zsize , dz)

u = np.fromfile('u.0000000', dtype=np.float64).reshape(ktot, jtot, itot)
v = np.fromfile('v.0000000', dtype=np.float64).reshape(ktot, jtot, itot)
w = np.fromfile('w.0000000', dtype=np.float64).reshape(ktot, jtot, itot)
rhoref_raw = np.fromfile('rhoref.0000000', dtype=np.float64)
rhoref = rhoref_raw[:ktot]
rhorefh = rhoref_raw[ktot:]

u_west = 2.1
u_east = 2.
v_south = 2.
v_north = 2.

u[:, :, :] = (u_west + (u_east - u_west) * xh[None, None, :]/xsize) * z[:, None, None]/zsize
v[:, :, :] = (v_south + (v_north - v_south) * yh[None, :, None]/ysize) * z[:, None, None]/zsize

for k in range(1, ktot):
    hor_div = (u_east - u_west) / xsize + (v_north - v_south) / ysize
    w[k, :, :] = (rhorefh[k]*w[k-1, :, :] - rhoref[k]*hor_div*dz) / rhorefh[k+1]

u.tofile('u.0000000')
v.tofile('v.0000000')
w.tofile('w.0000000')

print(w[-1, :, :])

