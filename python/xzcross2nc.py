import numpy
import struct
import netCDF4

nx    = 96
ny    = 64
nz    = 48
niter = 30

iterstart = 0
iterstep  = 100

crossname = "slngrad"

# load the dimensions
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
 
crossfile = netCDF4.Dataset("{}.xzcross.nc".format(crossname), "w")
# create dimensions in netCDF file
dim_x = crossfile.createDimension('x', nx)
dim_z = crossfile.createDimension('z', nz)
dim_t = crossfile.createDimension('t', niter)

# create dimension variables
var_t = crossfile.createVariable('time','f4',('t',))
var_x = crossfile.createVariable('x'   ,'f4',('x',))
var_z = crossfile.createVariable('z'   ,'f4',('z',))
var_t.units = "{0:07d} iterations since {1:07d}".format(niter*iterstep, iterstart)

# save the data
var_x[:] = x[:]
var_z[:] = z[:]

# create the variable
var_s = crossfile.createVariable('s','f4',('t','z','x',))

for i in range(niter):
  iter = iterstart + i*iterstep
  print("Processing iteration: ", iter)
  fin = open("{0:}.xzcross.{1:07d}".format(crossname,iter),"rb")
  raw = fin.read(nx*nz*8)
  tmp = numpy.array(struct.unpack('<{}d'.format(nx*nz), raw))
  del(raw)
  s = tmp.reshape((nz, nx))
  del(tmp)
  var_s[i,:,:] = s[:,:]
  del(s)
  fin.close()

