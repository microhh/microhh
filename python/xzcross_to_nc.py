import numpy   as np
import struct  as st
import netCDF4 as nc4

# Settings -------
variables  = ['u','v','w','th','p']
indexes    = [0,4]
nx         = 32
ny         = 32
nz         = 32
starttime  = 0
endtime    = 3600
sampletime = 1800
iotimeprec = 0.01
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
    dim_y  = crossfile.createDimension(locy,   np.size(indexes))
    dim_z  = crossfile.createDimension(locz,   nzsave)
    dim_t  = crossfile.createDimension('time', niter)
    
    # create dimension variables
    var_t  = crossfile.createVariable('time', sa, ('time',))
    var_x  = crossfile.createVariable(locx,   sa, (locx,  ))
    var_y  = crossfile.createVariable(locy,   sa, (locy,  ))
    var_z  = crossfile.createVariable(locz,   sa, (locz,  ))
    var_t.units = "Seconds since start of experiment"
    
    # save the data
    var_x[:]  = x[:nxsave] if locx=='x' else xh[:nxsave]
    var_z[:]  = z[:nxsave] if locz=='z' else zh[:nxsave]
    
    var_s = crossfile.createVariable(crossname, sa, ('time', locz, locx, locy,))
    
    for t in range(niter):
        for i in range(np.size(indexes)):
            index = indexes[i]
            otime = int((starttime + t*sampletime) / iotimeprec)
            print("Processing %5s, time=%7i, index=%4i"%(crossname, otime, index))
    
            var_t[t] = otime * iotimeprec
            var_y[i] = y[i] if locy=='y' else yh[i] 
    
            fin = open("{0:}.xz.{1:05d}.{2:07d}".format(crossname, index, otime), "rb")
            raw = fin.read(nx*nz*8)
            tmp = np.array(st.unpack('{0}{1}d'.format(en, nx*nz), raw))
            del(raw)
            s = tmp.reshape((nz, nx))
            del(tmp)
            var_s[t,:,:,i] = s[:nzsave,:nxsave]
            del(s)
    
            fin.close()
    crossfile.close() 
