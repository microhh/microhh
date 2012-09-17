import numpy
import struct
from pylab import *

# only read the first 2D slice to check the layout
nx = 64
ny = 64
nz = 64

fin = open("u.dump","rb")
raw = fin.read(nx*ny*nz*8)
tmp = numpy.array(struct.unpack('{}d'.format(nx*ny*nz), raw))
u   = tmp.reshape((nz, ny, nx))

pcolormesh(u[0,:,:])
colorbar()

