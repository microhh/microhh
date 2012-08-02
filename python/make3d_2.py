import numpy
import struct
import netCDF4

nx = 1024
ny = 1
nz = 512
nt = 32
tstart = 0
tstep  = 100

smin = -3.
smax =  3.

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

file = netCDF4.Dataset("s.nc","w")

dim_x = file.createDimension('x', nx)
dim_y = file.createDimension('y', ny)
dim_z = file.createDimension('z', nz)
dim_t = file.createDimension('time', nt)

var_x = file.createVariable('x','f8',('x',))
var_y = file.createVariable('y','f8',('y',))
var_z = file.createVariable('z','f8',('z',))
var_t = file.createVariable('time','i4',('time',))
var_s = file.createVariable('s','i2',('time','z','y','x',))
var_p = file.createVariable('p','f8',('time','z','y','x',))

var_t.units = "{0:07d} iterations since {1:07d}".format(nt*tstep, tstart)
var_s.scale_factor = (smax - smin) / 2.**15.
var_s.add_offset   = smin + (smax - smin) / 2.

var_x[:] = x[:]
var_y[:] = y[:]
var_z[:] = z[:]

# loop through the files
for t in range(nt):
  time = tstart + t*tstep
  print("Processing {:07d}".format(time))

  fin = open("s.{:07d}".format(time),"rb")
  raw = fin.read(n*8)
  tmp = numpy.array(struct.unpack('{}d'.format(n), raw))
  s   = tmp.reshape((nz, ny, nx))
  fin.close()
  var_s[t,:,:,:] = s[:,:,:]

  fin = open("p.{:07d}".format(time),"rb")
  raw = fin.read(n*8)
  tmp = numpy.array(struct.unpack('{}d'.format(n), raw))
  p   = tmp.reshape((nz, ny, nx))
  fin.close()
  var_p[t,:,:,:] = p[:,:,:]

  var_t[t] = t

file.close()

