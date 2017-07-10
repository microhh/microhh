import numpy   as np
import struct  as st
import netCDF4 as nc4

import microhh_tools as mht     # available in microhh/python directory

# Read the namelist settings
nl = mht.Read_namelist()

# Settings -------
variables  = nl['cross']['crosslist']
indexes    = -1   # With -1, script finds the correct indexes itself
nx         = nl['grid']['itot']
ny         = nl['grid']['jtot']
nz         = nl['grid']['ktot']
starttime  = 0
endtime    = nl['time']['endtime']
sampletime = nl['cross']['sampletime']
iotimeprec = 0
nxsave     = nx
nzsave     = nz
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

# calculate the number of iterations
niter = int((endtime-starttime) / sampletime + 1)

# load the dimensions
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

# Loop over the different variables
for crossname in variables:
    if indexes == -1:
        indexes_local = mht.get_cross_indices(crossname, 'xz')
    else:
        indexes_local = indexes

    crossfile = nc4.Dataset("{0}.xz.nc".format(crossname), "w")

    if(crossname == 'u'): loc = [1,0,0]
    elif(crossname=='v'): loc = [0,1,0]
    elif(crossname=='w'): loc = [0,0,1]
    else:                 loc = [0,0,0]

    locx = 'x' if loc[0] == 0 else 'xh'
    locy = 'y' if loc[1] == 0 else 'yh'
    locz = 'z' if loc[2] == 0 else 'zh'

    # create dimensions in netCDF file
    dim_x  = crossfile.createDimension(locx,   nxsave)
    dim_y  = crossfile.createDimension(locy,   np.size(indexes_local))
    dim_z  = crossfile.createDimension(locz,   nzsave)
    dim_t  = crossfile.createDimension('time', None)
    
    # create dimension variables
    var_t  = crossfile.createVariable('time', sa, ('time',))
    var_x  = crossfile.createVariable(locx,   sa, (locx,  ))
    var_y  = crossfile.createVariable(locy,   sa, (locy,  ))
    var_z  = crossfile.createVariable(locz,   sa, (locz,  ))
    var_t.units = "Seconds since start of experiment"
    
    # save the data
    var_x[:]  = x[:nxsave] if locx=='x' else xh[:nxsave]
    var_z[:]  = z[:nzsave] if locz=='z' else zh[:nzsave]
    
    var_s = crossfile.createVariable(crossname, sa, ('time', locz, locx, locy,))
    
    stop = False 
    for t in range(niter):
        if (stop):
            break
        for i in range(np.size(indexes_local)):
            index = indexes_local[i]
            otime = int((starttime + t*sampletime) / 10**iotimeprec)
            f_in  = "{0:}.xz.{1:05d}.{2:07d}".format(crossname, index, otime)

            try:
                fin = open(f_in, "rb")
            except:
                print('Stopping: cannot find file {}'.format(f_in))
                crossfile.sync()
                stop = True
                break
    
            print("Processing %8s, time=%7i, index=%4i"%(crossname, otime, index))

            var_t[t] = otime * 10**iotimeprec
            var_y[i] = y[index] if locy=='y' else yh[index] 
    
            raw = fin.read(nx*nz*8)
            tmp = np.array(st.unpack('{0}{1}d'.format(en, nx*nz), raw))
            del(raw)
            s = tmp.reshape((nz, nx))
            del(tmp)
            var_s[t,:,:,i] = s[:nzsave,:nxsave]
            del(s)
            fin.close()

    crossfile.close() 
