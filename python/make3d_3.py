import numpy
import struct
import netCDF4

var    = 'u'
nx     = 768
ny     = 384
nz     = 256
nt     = 1
tstart = 7200
tstep  = 100

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

file = netCDF4.Dataset("%s.nc"%var,"w")

if(var=='u'): loc = [1,0,0]
if(var=='v'): loc = [0,1,0]
if(var=='w'): loc = [0,0,1]

if(loc[0] == 0): dim_x  = file.createDimension('x', nx)  ; locx = 'x'
if(loc[1] == 0): dim_y  = file.createDimension('y', ny)  ; locy = 'y'
if(loc[2] == 0): dim_z  = file.createDimension('z', nz)  ; locz = 'z'
if(loc[0] == 1): dim_xh = file.createDimension('xh', nx) ; locx = 'xh'
if(loc[1] == 1): dim_yh = file.createDimension('yh', ny) ; locy = 'yh'
if(loc[2] == 1): dim_zh = file.createDimension('zh', nz) ; locz = 'zh'
dim_t  = file.createDimension('time', nt)

if(loc[0] == 0): var_x  = file.createVariable('x','f8',(locx,))
if(loc[1] == 0): var_y  = file.createVariable('y','f8',(locy,))
if(loc[2] == 0): var_z  = file.createVariable('z','f8',(locz,))
if(loc[0] == 1): var_xh = file.createVariable('x','f8',(locx,))
if(loc[1] == 1): var_yh = file.createVariable('y','f8',(locy,))
if(loc[2] == 1): var_zh = file.createVariable('z','f8',(locz,))
var_t  = file.createVariable('time','i4',('time',))
var_3d = file.createVariable(var,'f8',('time',locz,locy,locx,))

var_t.units = "{0:07d} iterations since {1:07d}".format(nt*tstep, tstart)

if(loc[0] == 0): var_x[:]  = x[:]
if(loc[1] == 0): var_y[:]  = y[:]
if(loc[2] == 0): var_z[:]  = z[:]
if(loc[0] == 1): var_xh[:] = xh[:]
if(loc[1] == 1): var_yh[:] = yh[:]
if(loc[2] == 1): var_zh[:] = zh[:]

# loop through the files
for t in range(nt):
  time = tstart + t*tstep
  print("Processing {:07d}".format(time))

  var_t[t] = t

  fin = open("%s.%07i"%(var,time),"rb")
  for k in range(nz):
    print("Processing %s, k=%i/%i"%(var,k+1,nz))
    raw = fin.read(nx*ny*8)
    tmp = numpy.array(struct.unpack('{}d'.format(nx*ny), raw))
    var_3d[t,k,:,:] = tmp.reshape((ny, nx))
  fin.close()

file.close()

