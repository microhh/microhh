import numpy
import struct
import netCDF4

nx = 512
ny = 512
nz = 512
crossname = "slngrad"

timestart = 0
timestep  = 10
timeend   = 2100
index     = 0

precision = 0.01
nxsave = 512
nzsave = 400

# calculate the number of iterations
niter = (timeend-timestart) / timestep + 1

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
 
ycross = y[index]
print('Creating cross at y = {0}'.format(ycross))

crossfile = netCDF4.Dataset("{0}.xz.{1:05d}.nc".format(crossname, index), "w")
# create dimensions in netCDF file
dim_x = crossfile.createDimension('x'   , nxsave)
dim_z = crossfile.createDimension('z'   , nzsave)
dim_t = crossfile.createDimension('time', niter)

# create dimension variables
var_t = crossfile.createVariable('time','f4',('time',))
var_x = crossfile.createVariable('x'   ,'f4',('x',))
var_z = crossfile.createVariable('z'   ,'f4',('z',))
var_t.units = "{0:07d} time since {1:07d}".format(timeend, timestart)

# save the data
var_x[:] = x[0:nxsave]
var_z[:] = z[0:nzsave]

# create the variable
var_s = crossfile.createVariable(crossname,'f4',('time','z','x',))

for i in range(niter):
  iter = timestart + i*timestep
  print("Processing time: ", iter)

  var_t[i] = iter*precision

  fin = open("{0:}.xz.{1:05d}.{2:07d}".format(crossname, index, iter),"rb")
  raw = fin.read(nx*nz*8)
  tmp = numpy.array(struct.unpack('<{}d'.format(nx*nz), raw))
  del(raw)
  s = tmp.reshape((nz, nx))
  del(tmp)
  var_s[i,:,:] = s[0:nzsave,0:nxsave]
  del(s)

  fin.close()

