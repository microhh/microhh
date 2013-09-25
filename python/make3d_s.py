import numpy
import struct
import netCDF4

nx = 1024
ny = 1024
nz = 1024
nxsave = 1024
nysave = 1024
nzsave = 192
nt = 8
tstart = 0
tstep  = 1

smin = -0.5
smax =  4.0

#pmin = -2.
#pmax =  2.

n = nx*ny*nz

fin = open("grid.{:07d}".format(0),"rb")
raw = fin.read(nx*8)
x   = numpy.array(struct.unpack('>{}d'.format(nx), raw))
raw = fin.read(nx*8)
xh  = numpy.array(struct.unpack('>{}d'.format(nx), raw))
raw = fin.read(ny*8)
y   = numpy.array(struct.unpack('>{}d'.format(ny), raw))
raw = fin.read(ny*8)
yh  = numpy.array(struct.unpack('>{}d'.format(ny), raw))
raw = fin.read(nz*8)
z   = numpy.array(struct.unpack('>{}d'.format(nz), raw))
raw = fin.read(nz*8)
zh  = numpy.array(struct.unpack('>{}d'.format(nz), raw))
fin.close()

file = netCDF4.Dataset("s.nc","w")

dim_x = file.createDimension('x', nxsave)
dim_y = file.createDimension('y', nysave)
dim_z = file.createDimension('z', nzsave)
dim_t = file.createDimension('time', nt)

var_x = file.createVariable('x','f8',('x',))
var_y = file.createVariable('y','f8',('y',))
var_z = file.createVariable('z','f8',('z',))
var_t = file.createVariable('time','f8',('time',))
var_s = file.createVariable('s','f4',('time','z','y','x',))
#var_s = file.createVariable('s','i2',('time','z','y','x',))
#var_s.scale_factor = (smax - smin) / 2.**15.
#var_s.add_offset   = smin + (smax - smin) / 2.

var_t.units = "{0:07d} iterations since {1:07d}".format((nt-1)*tstep, tstart)
var_x[:] = x[0:nxsave]
var_y[:] = y[0:nysave]
var_z[:] = z[0:nzsave]

# loop through the files
for t in range(nt):
  time = tstart + t*tstep
  print("Processing {:07d}".format(time))

  fin  = open("time.{:07d}".format(time),"rb")
  rawt = fin.read(8)
  tval = numpy.array(struct.unpack('>{}Q'.format(1), rawt))
  fin.close()
  tvalr = float(tval[0]) / 1.e6
  var_t[t] = tvalr

  fin = open("s.{:07d}".format(time),"rb")
  for k in range(nzsave):
    print("k = {}".format(k))
    raw = fin.read(nx*ny*8)
    tmp = numpy.array(struct.unpack('>{}d'.format(nx*ny), raw))
    s   = tmp.reshape((ny, nx))
    var_s[t,k,:,:] = s[0:nysave,0:nxsave]
  fin.close()

file.close()

