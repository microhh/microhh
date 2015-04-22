import numpy as np
import struct as st
import pylab as pl

nx = 64
ny = 1
nz = 128

N2 = 0.0004

n   = nx*nz
fin = open("grid.{:07d}".format(0),"rb")
raw = fin.read(nx*8)
x   = np.array(st.unpack('{0}d'.format(nx), raw))
raw = fin.read(nx*8)
xh  = np.array(st.unpack('{0}d'.format(nx), raw))
raw = fin.read(ny*8)
y   = np.array(st.unpack('{0}d'.format(ny), raw))
raw = fin.read(ny*8)
yh  = np.array(st.unpack('{0}d'.format(ny), raw))
raw = fin.read(nz*8)
z   = np.array(st.unpack('{0}d'.format(nz), raw))
raw = fin.read(nz*8)
zh  = np.array(st.unpack('{0}d'.format(nz), raw))
fin.close()

fin = open("u.0010000", "rb")
raw = fin.read(nx*nz*8)
tmp = np.array(st.unpack('{0}d'.format(nx*nz), raw))
del(raw)
tmp = tmp.reshape((nz, nx))
u = np.roll(tmp, nx/4, 1)
del(tmp)
 
fin = open("w.0010000", "rb")
raw = fin.read(nx*nz*8)
tmp = np.array(st.unpack('{0}d'.format(nx*nz), raw))
del(raw)
tmp = tmp.reshape((nz, nx))
w = np.roll(tmp, nx/4, 1)
del(tmp)
 
fin = open("b.0010000", "rb")
raw = fin.read(nx*nz*8)
tmp = np.array(st.unpack('{0}d'.format(nx*nz), raw))
del(raw)
tmp = tmp.reshape((nz, nx))
b = np.roll(tmp, nx/4, 1)
del(tmp)
for k in range(nz):
    b[k,:] -= N2 * z[k]

# Load reference data
nx_hires = 513
nz_hires = 1025
fin = open("u.out", "rb")
raw = fin.read(nx_hires*nz_hires*8)
tmp = np.array(st.unpack('{0}d'.format(nx_hires*nz_hires), raw))
del(raw)
u_hires = tmp.reshape((nz_hires, nx_hires))
del(tmp)
 
fin = open("w.out", "rb")
raw = fin.read(nx_hires*nz_hires*8)
tmp = np.array(st.unpack('{0}d'.format(nx_hires*nz_hires), raw))
del(raw)
w_hires = tmp.reshape((nz_hires, nx_hires))
del(tmp)
 
fin = open("b.out", "rb")
raw = fin.read(nx_hires*nz_hires*8)
tmp = np.array(st.unpack('{0}d'.format(nx_hires*nz_hires), raw))
del(raw)
b_hires = tmp.reshape((nz_hires, nx_hires))
del(tmp)

# Get data for 64 x 128 case
stepx   = nx_hires / nx
startxh = 0
startx  = stepx / 2
stepz   = nz_hires / nz
startz  = stepz / 2
startzh = 0
u_low = u_hires[startz :nz_hires-1:stepz, startxh:nx_hires-1:stepx]
w_low = w_hires[startzh:nz_hires-1:stepz, startx :nx_hires-1:stepx]
b_low = b_hires[startz :nz_hires-1:stepz, startx :nx_hires-1:stepx]

# Visualisation of fields.
pl.close('all')

x_hires = np.linspace(0.,  5.12,  513)
z_hires = np.linspace(0., 10.24, 1025)

"""
pl.figure()
pl.pcolormesh(x_hires, z_hires, u_hires)
pl.xlim(0, 5.12)
pl.ylim(0, 5.12)
pl.colorbar()
"""
pl.figure(figsize=(10,6))
pl.subplot(131)
pl.pcolormesh(xh, z, u)
pl.xlim(0, 5.12)
#pl.ylim(0, 5.12)
pl.title('u')
pl.colorbar()
pl.subplot(132)
pl.pcolormesh(xh, z, u_low)
pl.xlim(0, 5.12)
#pl.ylim(0, 5.12)
pl.title('u_ref')
pl.colorbar()
pl.subplot(133)
pl.pcolormesh(xh, z, u - u_low)
pl.xlim(0, 5.12)
#pl.ylim(0, 5.12)
pl.title('u_error')
pl.colorbar()
pl.savefig('u.png')

pl.figure(figsize=(10,6))
pl.subplot(131)
pl.pcolormesh(x, zh, w)
pl.xlim(0, 5.12)
#pl.ylim(0, 5.12)
pl.title('w')
pl.colorbar()
pl.subplot(132)
pl.pcolormesh(x, zh, w_low)
pl.xlim(0, 5.12)
#pl.ylim(0, 5.12)
pl.title('w_ref')
pl.colorbar()
pl.subplot(133)
pl.pcolormesh(x, zh, w - w_low)
pl.xlim(0, 5.12)
#pl.ylim(0, 5.12)
pl.title('w_error')
pl.colorbar()
pl.savefig('w.png')

pl.figure(figsize=(10,6))
pl.subplot(131)
pl.pcolormesh(x, z, b)
pl.xlim(0, 5.12)
#pl.ylim(0, 5.12)
pl.title('b')
pl.colorbar()
pl.subplot(132)
pl.pcolormesh(x, z, b_low)
pl.xlim(0, 5.12)
#pl.ylim(0, 5.12)
pl.title('b_ref')
pl.colorbar()
pl.subplot(133)
pl.pcolormesh(x, z, b - b_low)
pl.xlim(0, 5.12)
#pl.ylim(0, 5.12)
pl.title('b_error')
pl.colorbar()
pl.savefig('b.png')
