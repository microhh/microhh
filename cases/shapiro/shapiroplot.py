import numpy as np
import struct as st
import pylab as pl

nx = 64
ny = 1
nz = 256

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

fin = open("u.0001000", "rb")
raw = fin.read(nx*nz*8)
tmp = np.array(st.unpack('{0}d'.format(nx*nz), raw))
del(raw)
tmp = tmp.reshape((nz, nx))
u = np.roll(tmp, nx/4, 1)
del(tmp)
 
fin = open("w.0001000", "rb")
raw = fin.read(nx*nz*8)
tmp = np.array(st.unpack('{0}d'.format(nx*nz), raw))
del(raw)
tmp = tmp.reshape((nz, nx))
w = np.roll(tmp, nx/4, 1)
del(tmp)
 
fin = open("b.0001000", "rb")
raw = fin.read(nx*nz*8)
tmp = np.array(st.unpack('{0}d'.format(nx*nz), raw))
del(raw)
tmp = tmp.reshape((nz, nx))
b = np.roll(tmp, nx/4, 1)
del(tmp)

for k in range(nz):
    b[k,:] -= N2 * z[k]

# Visualisation of fields.
pl.figure()
pl.pcolormesh(xh, z, u)
pl.colorbar()

pl.figure()
pl.pcolormesh(x, zh, w)
pl.colorbar()

pl.figure()
pl.pcolormesh(x, z, b)
pl.colorbar()
