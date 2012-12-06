import numpy
import struct
from pylab import *

nx = 96
ny = 64
nz = 48
iter = 3000

fin = open("grid.{:07d}".format(0),"rb")
raw = fin.read(nx*8)
x = numpy.array(struct.unpack('<{}d'.format(nx), raw))
raw = fin.read(nx*8)
xh = numpy.array(struct.unpack('<{}d'.format(nx), raw))
raw = fin.read(ny*8)
y = numpy.array(struct.unpack('<{}d'.format(ny), raw))
raw = fin.read(ny*8)
yh = numpy.array(struct.unpack('<{}d'.format(ny), raw))
raw = fin.read(nz*8)
z = numpy.array(struct.unpack('<{}d'.format(nz), raw))
raw = fin.read(nz*8)
zh = numpy.array(struct.unpack('<{}d'.format(nz), raw))
fin.close()
 
fin = open("s.xzcross.{:07d}".format(iter),"rb")
raw = fin.read(nx*nz*8)
tmp = numpy.array(struct.unpack('<{}d'.format(nx*nz), raw))
del(raw)
s = tmp.reshape((nz, nx))
del(tmp)
fin.close()

fin = open("slngrad.xzcross.{:07d}".format(iter),"rb")
raw = fin.read(nx*nz*8)
tmp = numpy.array(struct.unpack('<{}d'.format(nx*nz), raw))
del(raw)
slngrad = tmp.reshape((nz, nx))
del(tmp)
fin.close()

figure()
pcolormesh(x, z, s)

figure()
pcolormesh(x, z, slngrad)

