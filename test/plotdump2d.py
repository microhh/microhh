import numpy
import struct
from pylab import *

# only read the first 2D slice to check the layout
nx = 32
ny = 32

fin = open("u.dump","rb")
raw = fin.read(nx*ny*8)
tmp = numpy.array(struct.unpack('{}d'.format(nx*ny), raw))
u   = tmp.reshape((ny, nx))

pcolormesh(u)
colorbar()

