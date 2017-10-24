import numpy   as np
import struct  as st
import netCDF4 as nc4
import os

from microhh_tools import *

# Read namelist
n = Read_namelist()

# Settings -------
variables  = ['u','w']
nx         = n['grid']['itot']
ny         = n['grid']['jtot']
nz         = n['grid']['ktot']
starttime  = 0
sampletime = n['time']['savetime']
endtime    = n['time']['endtime']
iotimeprec = n['time']['iotimeprec']
endian     = 'little'
savetype   = 'float'
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
nt = int((endtime - starttime) / sampletime + 1)

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

for variable in variables:

    # Create netCDF file
    ncfile = nc4.Dataset("%s.spanwise.nc"%variable, "w")
    
    if  (variable=='u'): loc = [1,0,0]
    elif(variable=='v'): loc = [0,1,0]
    elif(variable=='w'): loc = [0,0,1]
    else:                loc = [0,0,0]
    
    locx = 'x' if loc[0] == 0 else 'xh'
    locy = 'y' if loc[1] == 0 else 'yh'
    locz = 'z' if loc[2] == 0 else 'zh'
    
    dim_x  = ncfile.createDimension(locx,  nx)
    dim_z  = ncfile.createDimension(locz,  nz)
    dim_t  = ncfile.createDimension('time', None)
    
    var_x  = ncfile.createVariable(locx, sa, (locx,))
    var_z  = ncfile.createVariable(locz, sa, (locz,))
    var_t  = ncfile.createVariable('time', 'i4', ('time',))
    var_2d = ncfile.createVariable(variable, sa, ('time',locz, locx,))
    
    var_t.units = "time units since start"
    
    # Write grid properties to netCDF
    var_x[:] = x[:] if locx=='x' else xh[:]
    var_z[:] = z[:] if locz=='z' else zh[:]
    
    # Loop through the files and read 3d field
    for t in range(nt):
        time = int((starttime + t*sampletime) / 10**iotimeprec)
        print("Processing t={:07d}".format(time))
    
        name = "%s.%07i"%(variable, time)
        
        if (os.path.exists(name)):
            var_t[t] = time * 10**iotimeprec

            fin = open("%s.%07i"%(variable, time),"rb")
            for k in range(nz):
                raw = fin.read(nx*ny*8)
                tmp = np.array(st.unpack('{0}{1}d'.format(en, nx*ny), raw)).reshape((ny, nx))
                var_2d[t,k,:] = tmp.mean(axis=0)
            fin.close()
        else:
            print('Can\'t find file {}, breaking..'.format(name))
            ncfile.sync()
            break
    
    ncfile.close()
