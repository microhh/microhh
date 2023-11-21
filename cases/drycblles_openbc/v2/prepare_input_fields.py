import matplotlib.pyplot as pl
import numpy as np
import xarray as xr
from numba import jit

# from interpolation_tools import calc_w_from_uv
dtype = np.float32
ngc = 3

pl.close('all')

@jit(nopython=True, nogil=True, fastmath=True)
def calc_w_from_uv(
        w, u, v,
        u_lbc_west,
        u_lbc_east,
        v_lbc_south,
        v_lbc_north,
        rho, rhoh, dz,
        dxi, dyi,
        itot, jtot, ktot):

    for j in range(jtot):
        for i in range(itot):
            w[0,j,i] = 0.

    for k in range(ktot):
        for j in range(jtot):
            for i in range(itot):
                if i == 0:
                    um = u_lbc_west[k,j]
                else:
                    um = u[k,j,i]

                if i == itot-1:
                    up = u_lbc_east[k,j]
                else:
                    up = u[k,j,i+1]

                if j == 0:
                    vm = v_lbc_south[k,i]
                else:
                    vm = v[k,j,i]

                if j == jtot-1:
                    vp = v_lbc_north[k,i]
                else:
                    vp = v[k,j+1,i]

                w[k+1,j,i] = -(rho[k] * ((up - um) * dxi + (vp - vm) * dyi) * dz[k] - rhoh[k] * w[k,j,i]) / rhoh[k+1]


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

u = np.fromfile('u.0000000', dtype=dtype).reshape(ktot, jtot, itot)
v = np.fromfile('v.0000000', dtype=dtype).reshape(ktot, jtot, itot)
w = np.zeros((ktot+1, jtot, itot), dtype=dtype)

rhoref_raw = np.fromfile('rhoref.0000000', dtype=dtype)
rhoref = rhoref_raw[:ktot]
rhorefh = rhoref_raw[ktot:]

lbc = xr.open_dataset('drycblles_lbc_input.nc')

u_west = lbc.u_west[0,:,ngc:-ngc,-1].values
u_east = lbc.u_east[0,:,ngc:-ngc,0].values

v_south = lbc.v_south[0,:,-1,ngc:-ngc].values
v_north = lbc.v_north[0,:,0,ngc:-ngc].values

u[:, :, :] = (u_west[:,:,None] + (u_east[:,:,None] - u_west[:,:,None])    * xh[None,None,:] / xsize)
v[:, :, :] = (v_south[:,None,:] + (v_north[:,None,:] - v_south[:,None,:]) * yh[None,:,None] / ysize)

dz = np.ones(ktot)*dz

calc_w_from_uv(
        w, u, v,
        u_west,
        u_east,
        v_south,
        v_north,
        rhoref, rhorefh, dz,
        1/dx, 1/dy,
        itot, jtot, ktot)

u.tofile('u.0000000')
v.tofile('v.0000000')
w[:-1,:,:].tofile('w.0000000')

print('<w_top> = ', w[-1,:,:].mean())

pl.figure()
k = 10
pl.subplot(221)
pl.imshow(u[k,:,:], origin='lower')
pl.colorbar()

pl.subplot(222)
pl.imshow(v[k,:,:], origin='lower')
pl.colorbar()

pl.subplot(223)
pl.imshow(w[k,:,:], origin='lower')
pl.colorbar()

pl.subplot(224)
pl.imshow(w[-1,:,:], origin='lower')
pl.colorbar()
