import numpy   as np
import struct  as st
import netCDF4 as nc4

# Settings -------
var      = 'w'
nx       = 32
ny       = 32
nz       = 32
tstart   = 0
tstep    = 1800
tend     = 10800
nxsave   = nx
nysave   = ny
nzsave   = nz
endian   = 'little'
savetype = 'float'
# End settings ---

# Set the correct string for the endianness
if (endian == 'little'):
    en = '<'
elif (endian == 'big'):
    en = '>'
else:
    raise RuntimeError("Endianness has to be little or big")

# Set the correct string for the savetype
if (savetype == 'double'):
    sa = 'f8'
elif (savetype == 'float'):
    sa = 'f4'
else:
    raise RuntimeError("The savetype has to be float or double")

# Calculate the number of time steps
nt = (tend - tstart) / tstep + 1

# Read grid properties from grid.0000000
n   = nx*ny*nz
fin = open("grid.{:07d}".format(0),"rb")
raw = fin.read(nx*8)
x   = np.array(st.unpack('{0}{1}d'.format(en, nx), raw))
raw = fin.read(nx*8)
xh  = np.array(st.unpack('{0}{1}d'.format(en, nx), raw))
raw = fin.read(ny*8)
y   = np.array(st.unpack('{0}{1}d'.format(en, ny), raw))
raw = fin.read(ny*8)
yh  = np.array(st.unpack('{0}{1}d'.format(en, ny), raw))
raw = fin.read(nz*8)
z   = np.array(st.unpack('{0}{1}d'.format(en, nz), raw))
raw = fin.read(nz*8)
zh  = np.array(st.unpack('{0}{1}d'.format(en, nz), raw))
fin.close()

# Create netCDF file
ncfile = nc4.Dataset("%s.nc"%var, "w")

if  (var=='u'): loc = [1,0,0]
elif(var=='v'): loc = [0,1,0]
elif(var=='w'): loc = [0,0,1]
else:           loc = [0,0,0]

if(loc[0] == 0): dim_x  = ncfile.createDimension('x',  nxsave) ; locx = 'x'
if(loc[1] == 0): dim_y  = ncfile.createDimension('y',  nysave) ; locy = 'y'
if(loc[2] == 0): dim_z  = ncfile.createDimension('z',  nzsave) ; locz = 'z'
if(loc[0] == 1): dim_xh = ncfile.createDimension('xh', nxsave) ; locx = 'xh'
if(loc[1] == 1): dim_yh = ncfile.createDimension('yh', nysave) ; locy = 'yh'
if(loc[2] == 1): dim_zh = ncfile.createDimension('zh', nzsave) ; locz = 'zh'
dim_t  = ncfile.createDimension('time', nt)

if(loc[0] == 0): var_x  = ncfile.createVariable('x', '{}'.format(sa), (locx,))
if(loc[1] == 0): var_y  = ncfile.createVariable('y', '{}'.format(sa), (locy,))
if(loc[2] == 0): var_z  = ncfile.createVariable('z', '{}'.format(sa), (locz,))
if(loc[0] == 1): var_xh = ncfile.createVariable('xh', '{}'.format(sa), (locx,))
if(loc[1] == 1): var_yh = ncfile.createVariable('yh', '{}'.format(sa), (locy,))
if(loc[2] == 1): var_zh = ncfile.createVariable('zh', '{}'.format(sa), (locz,))
var_t  = ncfile.createVariable('time', 'i4', ('time',))
var_3d = ncfile.createVariable(var, '{}'.format(sa), ('time',locz,locy,locx,))

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
        tmp = np.array(st.unpack('{0}{1}d'.format(en, nx*ny), raw))
        var_3d[t,k,:,:] = tmp.reshape((ny, nx))[:nysave,:nxsave]
    fin.close()
    ncfile.sync()

ncfile.close()
