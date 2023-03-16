import numpy as np

xsize, ysize, zsize = 3200., 3200., 3200
itot, jtot, ktot = 64, 64, 64

dx, dy, dz = xsize/itot, ysize/jtot, zsize/ktot

xh = np.arange(0, xsize - dx/2, dx)
x = np.arange(dx/2, xsize , dx)
yh = np.arange(0, ysize - dy/2, dy)
y = np.arange(dy/2, ysize , dy)
zh = np.arange(0, zsize - dz/2, dz)
z = np.arange(dz/2, zsize , dz)

u = np.fromfile('u.0000000', dtype=np.float32).reshape(ktot, jtot, itot)
v = np.fromfile('v.0000000', dtype=np.float32).reshape(ktot, jtot, itot)
w = np.fromfile('w.0000000', dtype=np.float32).reshape(ktot, jtot, itot)
rhoref_raw = np.fromfile('rhoref.0000000', dtype=np.float32)
rhoref = rhoref_raw[:ktot]
rhorefh = rhoref_raw[ktot:]

u_west = 2.1
u_east = 2

u[:, :, :] = u_west + (u_east - u_west) * xh[None, None, :]/xsize
v[:, :, :] = 2

for k in range(1, ktot):
    dudx = (u_east - u_west) / xsize
    w[k, :, :] = (rhorefh[k]*w[k-1, :, :] - rhoref[k]*dudx*dz) / rhorefh[k+1]

u.tofile('u.0000000')
v.tofile('v.0000000')
w.tofile('w.0000000')

print(w[:, 0, 0])

