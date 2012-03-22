from pylab import *
from struct import *

itot = 32
jtot = 32
ktot = 32

iter = 0

### DO NOT EDIT BELOW ###

igc  = 1
jgc  = 1
kgc  = 1

nx = itot+2*igc
ny = jtot+2*jgc
nz = ktot+2*kgc

n  = nx*ny*nz

fin = open("grid.{:06d}".format(0),"rb")
raw = fin.read(nx*8)
x   = array(unpack('{}d'.format(nx), raw))
raw = fin.read(nx*8)
xh  = array(unpack('{}d'.format(nx), raw))
raw = fin.read(ny*8)
y   = array(unpack('{}d'.format(ny), raw))
raw = fin.read(ny*8)
yh  = array(unpack('{}d'.format(ny), raw))
raw = fin.read(nz*8)
z   = array(unpack('{}d'.format(nz), raw))
raw = fin.read(nz*8)
zh  = array(unpack('{}d'.format(nz), raw))
fin.close()

fin = open("u.{:06d}".format(iter),"rb")
raw = fin.read(n*8)
tmp = array(unpack('{}d'.format(n), raw))
u   = tmp.reshape((nz, ny, nx))
fin.close()

fin = open("v.{:06d}".format(iter),"rb")
raw = fin.read(n*8)
tmp = array(unpack('{}d'.format(n), raw))
v   = tmp.reshape((nz, ny, nx))
fin.close()

fin = open("w.{:06d}".format(iter),"rb")
raw = fin.read(n*8)
tmp = array(unpack('{}d'.format(n), raw))
w   = tmp.reshape((nz, ny, nx))
fin.close()

fin = open("p.{:06d}".format(iter),"rb")
raw = fin.read(n*8)
tmp = array(unpack('{}d'.format(n), raw))
p   = tmp.reshape((nz, ny, nx))
fin.close()

fin = open("s.{:06d}".format(iter),"rb")
raw = fin.read(n*8)
tmp = array(unpack('{}d'.format(n), raw))
s   = tmp.reshape((nz, ny, nx))
fin.close()

