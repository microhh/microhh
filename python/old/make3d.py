import numpy
import struct
import netCDF4

nx = 256
ny = 192
nz = 128
iter = 0

n = nx*ny*nz

fin = open("grid.{:07d}".format(0),"rb")
raw = fin.read(nx*8)
x   = numpy.array(struct.unpack('{}d'.format(nx), raw))
raw = fin.read(nx*8)
xh  = numpy.array(struct.unpack('{}d'.format(nx), raw))
raw = fin.read(ny*8)
y   = numpy.array(struct.unpack('{}d'.format(ny), raw))
raw = fin.read(ny*8)
yh  = numpy.array(struct.unpack('{}d'.format(ny), raw))
raw = fin.read(nz*8)
z   = numpy.array(struct.unpack('{}d'.format(nz), raw))
raw = fin.read(nz*8)
zh  = numpy.array(struct.unpack('{}d'.format(nz), raw))
fin.close()

fin = open("u.{:07d}".format(iter),"rb")
raw = fin.read(n*8)
tmp = numpy.array(struct.unpack('{}d'.format(n), raw))
um  = numpy.min(tmp)
up  = numpy.max(tmp)
u   = tmp.reshape((nz, ny, nx))
fin.close()

file = netCDF4.Dataset("u.{:07d}.nc".format(iter),"w")

dim_xh = file.createDimension('xh', nx)
dim_y  = file.createDimension('y' , ny)
dim_z  = file.createDimension('z' , nz)

var_xh = file.createVariable('xh','f8',('xh',))
var_y  = file.createVariable('y' ,'f8',('y',))
var_z  = file.createVariable('z' ,'f8',('z',))
var_u  = file.createVariable('u' ,'i2',('z','y','xh',))

var_u.scale_factor = (up - um) / 2.**15.
var_u.add_offset   = um + (up - um) / 2.

var_xh[:] = xh[:]
var_y [:] = y [:]
var_z [:] = z [:]
var_u [:] = u [:,:,:]

file.close()

