import numpy
import struct
import netCDF4

nx = 512
ny = 512
nz = 512
crossname = "w"

timestart = 0
timestep  = 10
timeend   = 2030
index     = 150

precision = 0.01
nxsave = 512
nysave = 512

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

if(crossname == "w"):
  zcross = zh[index]
else:
  zcross = z[index]
print('Creating cross at z = {0}'.format(zcross))
 
crossfile = netCDF4.Dataset("{0}.xy.{1:05d}.nc".format(crossname, index), "w")
# create dimensions in netCDF file
dim_x = crossfile.createDimension('x'   , nxsave)
dim_y = crossfile.createDimension('y'   , nysave)
dim_t = crossfile.createDimension('time', niter)

# create dimension variables
var_t = crossfile.createVariable('time','f4',('time',))
var_x = crossfile.createVariable('x'   ,'f4',('x',))
var_y = crossfile.createVariable('y'   ,'f4',('y',))
var_t.units = "{0:07d} time since {1:07d}".format(timeend, timestart)

# save the data
var_x[:] = x[0:nxsave]
var_y[:] = y[0:nysave]

# create the variable
var_s = crossfile.createVariable(crossname,'f4',('time','y','x',))

for i in range(niter):
  iter = timestart + i*timestep
  print("Processing time: ", iter)

  var_t[i] = iter*precision

  fin = open("{0:}.xy.{1:05d}.{2:07d}".format(crossname, index, iter),"rb")
  raw = fin.read(nx*ny*8)
  tmp = numpy.array(struct.unpack('<{}d'.format(nx*ny), raw))
  del(raw)
  s = tmp.reshape((ny, nx))
  del(tmp)
  var_s[i,:,:] = s[0:nysave,0:nxsave]
  del(s)

  fin.close()

