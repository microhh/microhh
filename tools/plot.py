from pylab import *
from struct import *
nx = 34
ny = 34
nz = 34
n  = nx*ny*nz

fin = open("u","rb")
raw = fin.read(n*8)
tmp = array(unpack('{}d'.format(n), raw))
u   = tmp.reshape((nz, ny, nx))

fin = open("v","rb")
raw = fin.read(n*8)
tmp = array(unpack('{}d'.format(n), raw))
v   = tmp.reshape((nz, ny, nx))

fin = open("w","rb")
raw = fin.read(n*8)
tmp = array(unpack('{}d'.format(n), raw))
w   = tmp.reshape((nz, ny, nx))

fin = open("w","rb")
raw = fin.read(n*8)
tmp = array(unpack('{}d'.format(n), raw))
p   = tmp.reshape((nz, ny, nx))

