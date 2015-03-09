import numpy
import struct
import netCDF4

var    = 'w'
nx     = 32
ny     = 32
nz     = 32
nt     = 5
tstart = 0
tstep  = 1800

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

ncfile = netCDF4.Dataset("%s.nc"%var, "w")

if(var=='u'): loc = [1,0,0]
if(var=='v'): loc = [0,1,0]
if(var=='w'): loc = [0,0,1]
else:         loc = [0,0,0]

if(loc[0] == 0): dim_x  = ncfile.createDimension('x', nx)  ; locx = 'x'
if(loc[1] == 0): dim_y  = ncfile.createDimension('y', ny)  ; locy = 'y'
if(loc[2] == 0): dim_z  = ncfile.createDimension('z', nz)  ; locz = 'z'
if(loc[0] == 1): dim_xh = ncfile.createDimension('xh', nx) ; locx = 'xh'
if(loc[1] == 1): dim_yh = ncfile.createDimension('yh', ny) ; locy = 'yh'
if(loc[2] == 1): dim_zh = ncfile.createDimension('zh', nz) ; locz = 'zh'
dim_t  = ncfile.createDimension('time', nt)

if(loc[0] == 0): var_x  = ncfile.createVariable('x','f8',(locx,))
if(loc[1] == 0): var_y  = ncfile.createVariable('y','f8',(locy,))
if(loc[2] == 0): var_z  = ncfile.createVariable('z','f8',(locz,))
if(loc[0] == 1): var_xh = ncfile.createVariable('x','f8',(locx,))
if(loc[1] == 1): var_yh = ncfile.createVariable('y','f8',(locy,))
if(loc[2] == 1): var_zh = ncfile.createVariable('z','f8',(locz,))
var_t  = ncfile.createVariable('time','i4',('time',))
var_3d = ncfile.createVariable(var,'f8',('time',locz,locy,locx,))

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

ncfile.close()
