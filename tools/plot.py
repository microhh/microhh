from pylab import *
from struct import *

itot = 32
jtot = 32
ktot = 32

igc  = 1
jgc  = 1
kgc  = 1

nx = itot+2*igc
ny = jtot+2*jgc
nz = ktot+2*kgc

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

fin = open("p","rb")
raw = fin.read(n*8)
tmp = array(unpack('{}d'.format(n), raw))
p   = tmp.reshape((nz, ny, nx))

