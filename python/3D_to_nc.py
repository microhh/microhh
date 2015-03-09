import numpy   as np
import struct  as st
import netCDF4 as nc4

# Settings -------
var    = 'p'
nx     = 32
ny     = 32
nz     = 32
nt     = 5
tstart = 0
tstep  = 1800
nxsave = nx
nysave = ny
nzsave = nz
# End settings ---

# Read grid properties from grid.0000000
n   = nx*ny*nz
fin = open("grid.{:07d}".format(0),"rb")
raw = fin.read(nx*8)
x   = np.array(st.unpack('{}d'.format(nx), raw))
raw = fin.read(nx*8)
xh  = np.array(st.unpack('{}d'.format(nx), raw))
raw = fin.read(ny*8)
y   = np.array(st.unpack('{}d'.format(ny), raw))
raw = fin.read(ny*8)
yh  = np.array(st.unpack('{}d'.format(ny), raw))
raw = fin.read(nz*8)
z   = np.array(st.unpack('{}d'.format(nz), raw))
raw = fin.read(nz*8)
zh  = np.array(st.unpack('{}d'.format(nz), raw))
fin.close()

# Create netCDF file
ncfile = nc4.Dataset("%s.nc"%var, "w")

if(var=='u'): loc = [1,0,0]
if(var=='v'): loc = [0,1,0]
if(var=='w'): loc = [0,0,1]
else:         loc = [0,0,0]

if(loc[0] == 0): dim_x  = ncfile.createDimension('x',  nxsave) ; locx = 'x'
if(loc[1] == 0): dim_y  = ncfile.createDimension('y',  nysave) ; locy = 'y'
if(loc[2] == 0): dim_z  = ncfile.createDimension('z',  nzsave) ; locz = 'z'
if(loc[0] == 1): dim_xh = ncfile.createDimension('xh', nxsave) ; locx = 'xh'
if(loc[1] == 1): dim_yh = ncfile.createDimension('yh', nysave) ; locy = 'yh'
if(loc[2] == 1): dim_zh = ncfile.createDimension('zh', nzsave) ; locz = 'zh'
dim_t  = ncfile.createDimension('time', nt)

if(loc[0] == 0): var_x  = ncfile.createVariable('x', 'f8', (locx,))
if(loc[1] == 0): var_y  = ncfile.createVariable('y', 'f8', (locy,))
if(loc[2] == 0): var_z  = ncfile.createVariable('z', 'f8', (locz,))
if(loc[0] == 1): var_xh = ncfile.createVariable('x', 'f8', (locx,))
if(loc[1] == 1): var_yh = ncfile.createVariable('y', 'f8', (locy,))
if(loc[2] == 1): var_zh = ncfile.createVariable('z', 'f8', (locz,))
var_t  = ncfile.createVariable('time', 'i4', ('time',))
var_3d = ncfile.createVariable(var, 'f8', ('time',locz,locy,locx,))

var_t.units = "time units since start"

# Write grid properties to netCDF
if(loc[0] == 0): var_x[:]  = x[:nxsave]
if(loc[1] == 0): var_y[:]  = y[:nysave]
if(loc[2] == 0): var_z[:]  = z[:nzsave]
if(loc[0] == 1): var_xh[:] = xh[:nxsave]
if(loc[1] == 1): var_yh[:] = yh[:nysave]
if(loc[2] == 1): var_zh[:] = zh[:nzsave]

# Loop through the files and read 3d field
for t in range(nt):
    time = tstart + t*tstep
    print("Processing t={:07d}".format(time))

    var_t[t] = time

    fin = open("%s.%07i"%(var,time),"rb")
    for k in range(nzsave):
        raw = fin.read(nx*ny*8)
        tmp = np.array(st.unpack('{}d'.format(nx*ny), raw))
        var_3d[t,k,:,:] = tmp.reshape((ny, nx))[:nysave,:nxsave]
    fin.close()

ncfile.close()
