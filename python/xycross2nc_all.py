import numpy as np
import struct
import netCDF4

nx = 96
ny = 64
nz = 48

crossnames  = (['u','v','w','p'])
heights     = ([7,34,127])

timestart   = 0
timestep    = 60
timeend     = 300

t_precision = 1 # Time precision output MicroHH
d_precision = 'f8' # Precision to save data in NetCDF

nxsave = nx
nysave = ny

# calculate the number of iterations
niter = (timeend-timestart) / timestep + 1

# load the dimensions
fin = open("grid.{:07d}".format(0),"rb")
raw = fin.read(nx*8)
x = np.array(struct.unpack('<{}d'.format(nx), raw))
raw = fin.read(nx*8)
xh = np.array(struct.unpack('<{}d'.format(nx), raw))
raw = fin.read(ny*8)
y = np.array(struct.unpack('<{}d'.format(ny), raw))
raw = fin.read(ny*8)
yh = np.array(struct.unpack('<{}d'.format(ny), raw))
raw = fin.read(nz*8)
z = np.array(struct.unpack('<{}d'.format(nz), raw))
raw = fin.read(nz*8)
zh = np.array(struct.unpack('<{}d'.format(nz), raw))
fin.close()

crossfile = netCDF4.Dataset("xy.nc","w")
# create dimensions in netCDF file
dim_x  = crossfile.createDimension('x'   , nxsave)
dim_y  = crossfile.createDimension('y'   , nysave)
dim_z  = crossfile.createDimension('z'   , np.size(heights))
dim_zh = crossfile.createDimension('zh'  , np.size(heights))
dim_t  = crossfile.createDimension('time', niter)

# create dimension variables
var_t  = crossfile.createVariable('time' , d_precision, ('time',))
var_x  = crossfile.createVariable('x'    , d_precision, ('x',))
var_y  = crossfile.createVariable('y'    , d_precision, ('y',))
var_z  = crossfile.createVariable('z'    , d_precision, ('z',))
var_zh = crossfile.createVariable('zh'   , d_precision, ('zh',))
var_t.units = "{0:07d} seconds since {1:07d}".format(timeend, timestart)

# save the data
var_x[:] = x[0:nxsave]
var_y[:] = y[0:nysave]

# create the variable
var_u = crossfile.createVariable("u", d_precision,('time','z' ,'y','x',)) if "u" in crossnames else 0
var_v = crossfile.createVariable("v", d_precision,('time','z' ,'y','x',)) if "v" in crossnames else 0
var_w = crossfile.createVariable("w", d_precision,('time','zh','y','x',)) if "w" in crossnames else 0
var_p = crossfile.createVariable("p", d_precision,('time','z' ,'y','x',)) if "p" in crossnames else 0

for i in range(niter):
  iter = timestart + i*timestep
  print("Processing time: %i"%iter)

  var_t[i] = iter*t_precision

  for crossname in crossnames:
    print(".. Processing var: %s"%crossname)
    for k in range(np.size(heights)):
      index = heights[k]

      var_z[k]  = z[k]
      var_zh[k] = zh[k] 

      fin = open("{0:}.xy.{1:05d}.{2:07d}".format(crossname, index, iter),"rb")
      raw = fin.read(nx*ny*8)
      tmp = np.array(struct.unpack('<{}d'.format(nx*ny), raw))
      del(raw)
      s = tmp.reshape((ny, nx))
      del(tmp)

      if(crossname == "u"):
        var_u[i,k,:,:] = s[0:nysave,0:nxsave]
      if(crossname == "v"):
        var_v[i,k,:,:] = s[0:nysave,0:nxsave]
      if(crossname == "w"):
        var_w[i,k,:,:] = s[0:nysave,0:nxsave]
      if(crossname == "p"):
        var_p[i,k,:,:] = s[0:nysave,0:nxsave]

      del(s)

  fin.close()

