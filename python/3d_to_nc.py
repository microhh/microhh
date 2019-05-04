import numpy   as np
import struct  as st
import netCDF4 as nc4
import os
import sys
import microhh_tools as mht     # available in microhh/python directory

if (len(sys.argv)>1):
  os.chdir(sys.argv[1])



# Read the namelist settings
nl = mht.Read_namelist()

# Settings -------
variables  = nl['dump']['dumplist']
nx         = nl['grid']['itot']
ny         = nl['grid']['jtot']
nz         = nl['grid']['ktot']
starttime  = nl['time']['starttime']
endtime    = nl['time']['endtime']
sampletime = nl['dump']['sampletime']
iotimeprec = 0
nxsave     = nx
nysave     = ny
nzsave     = nz
endian     = 'little'
loadtype   = 'float'
savetype   = 'float'
# End settings ---

# Set the correct string for the endianness
if (endian == 'little'):
    en = '<'
elif (endian == 'big'):
    en = '>'
else:
    raise RuntimeError("Endianness has to be little or big")

# Set the correct string for the loadtype
if (loadtype == 'double'):
    TF = 8
    ra = 'd'
elif (loadtype == 'float'):
    TF = 4
    ra = 'f'
else:
    raise RuntimeError("The savetype has to be float or double")

fstring = '{0}{1}'+ra

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
raw = fin.read(nx*TF)
x   = np.array(st.unpack(fstring.format(en, nx), raw))
raw = fin.read(nx*TF)
xh  = np.array(st.unpack(fstring.format(en, nx), raw))
raw = fin.read(ny*TF)
y   = np.array(st.unpack(fstring.format(en, ny), raw))
raw = fin.read(ny*TF)
yh  = np.array(st.unpack(fstring.format(en, ny), raw))
raw = fin.read(nz*TF)
z   = np.array(st.unpack(fstring.format(en, nz), raw))
raw = fin.read(nz*TF)
zh  = np.array(st.unpack(fstring.format(en, nz), raw))
fin.close()


# Loop over the different variables
for variable in variables:
  # Create netCDF file
  try:
    ncfile = nc4.Dataset("%s.nc"%variable, "w")

    if  (variable=='u'): loc = [1,0,0]
    elif(variable=='v'): loc = [0,1,0]
    elif(variable=='w'): loc = [0,0,1]
    else:                loc = [0,0,0]

    locx = 'x' if loc[0] == 0 else 'xh'
    locy = 'y' if loc[1] == 0 else 'yh'
    locz = 'z' if loc[2] == 0 else 'zh'

    dim_x  = ncfile.createDimension(locx,  nxsave)
    dim_y  = ncfile.createDimension(locy,  nysave)
    dim_z  = ncfile.createDimension(locz,  nzsave)
    dim_t  = ncfile.createDimension('time', nt)

    var_x  = ncfile.createVariable(locx, sa, (locx,))
    var_y  = ncfile.createVariable(locy, sa, (locy,))
    var_z  = ncfile.createVariable(locz, sa, (locz,))
    var_t  = ncfile.createVariable('time', 'i4', ('time',))
    var_3d = ncfile.createVariable(variable, sa, ('time',locz, locy, locx,),zlib=True)

    var_t.units = "time units since start"

    # Write grid properties to netCDF
    var_x[:] = x[:nxsave] if locx=='x' else xh[:nxsave]
    var_y[:] = y[:nysave] if locy=='y' else yh[:nysave]
    var_z[:] = z[:nzsave] if locz=='z' else zh[:nzsave]

    # Loop through the files and read 3d field
    for t in range(nt):
        time = int((starttime + t*sampletime) / 10**iotimeprec)
        print("Processing t={:07d}".format(time))

        var_t[t] = time * 10**iotimeprec

        fin = open("%s.%07i"%(variable, time),"rb")
        for k in range(nzsave):
            raw = fin.read(nx*ny*TF)
            tmp = np.array(st.unpack(fstring.format(en, nx*ny), raw))
            var_3d[t,k,:,:] = tmp.reshape((ny, nx))[:nysave,:nxsave]
        fin.close()
        ncfile.sync()

    ncfile.close()
  except:
    print("Failed to create %s.nc"%variable)
    ncfile.close()
