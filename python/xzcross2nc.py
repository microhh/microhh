import numpy
import struct
import netCDF4

nx    = 2048
ny    = 2048
nz    = 1024

B0 = 0.0032
N2 = 3.

iterstart = 100
iterstep  = 100
iterend   = 60000

nxsave = 2048
nzsave = 1024

niter = (iterend-iterstart) / iterstep + 1

#crossname = "s"
crossname = "slngrad"

# load the dimensions
fin = open("grid.{:07d}".format(0),"rb")
raw = fin.read(nx*8)
x = numpy.array(struct.unpack('>{}d'.format(nx), raw))
raw = fin.read(nx*8)
xh = numpy.array(struct.unpack('>{}d'.format(nx), raw))
raw = fin.read(ny*8)
y = numpy.array(struct.unpack('>{}d'.format(ny), raw))
raw = fin.read(ny*8)
yh = numpy.array(struct.unpack('>{}d'.format(ny), raw))
raw = fin.read(nz*8)
z = numpy.array(struct.unpack('>{}d'.format(nz), raw))
raw = fin.read(nz*8)
zh = numpy.array(struct.unpack('>{}d'.format(nz), raw))
fin.close()
 
crossfile = netCDF4.Dataset("{}.xzcross.nc".format(crossname), "w")
# create dimensions in netCDF file
dim_x = crossfile.createDimension('x'   , nxsave)
dim_z = crossfile.createDimension('z'   , nzsave)
dim_t = crossfile.createDimension('time', niter)

# create dimension variables
var_t = crossfile.createVariable('time','f4',('time',))
var_x = crossfile.createVariable('x'   ,'f4',('x',))
var_z = crossfile.createVariable('z'   ,'f4',('z',))
var_t.units = "{0:07d} iterations since {1:07d}".format(iterend, iterstart)

# save the data
var_x[:] = x[0:nxsave]
var_z[:] = z[0:nzsave]

# create the variable
var_s = crossfile.createVariable(crossname,'f4',('time','z','x',))

for i in range(niter):
  iter = iterstart + i*iterstep
  print("Processing iteration: ", iter)

  fin = open("time.{:07d}".format(iter),"rb")
  raw = fin.read(8)
  t = numpy.array(struct.unpack('>{}Q'.format(1), raw)) / 1.e6
  var_t[i] = t
  del(t)

  fin = open("{0:}.xzcross.{1:07d}".format(crossname,iter),"rb")
  raw = fin.read(nx*nz*8)
  tmp = numpy.array(struct.unpack('>{}d'.format(nx*nz), raw))
  del(raw)
  s = tmp.reshape((nz, nx))
  del(tmp)
  var_s[i,:,:] = s[0:nzsave,0:nxsave]
  del(s)

  fin.close()

