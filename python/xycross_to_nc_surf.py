import numpy   as np
import struct  as st
import netCDF4 as nc4

# Settings -------
variables  = ['bfluxbot']
nx         = 2048
ny         = 2048
nz         = 1024
starttime  = 0.
endtime    = 10.
sampletime = 0.1
iotimeprec = -1
nxsave     = nx
nysave     = ny
endian     = 'big'
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
    crossfile = nc4.Dataset("{0}.xy.nc".format(crossname), "w")

    if(crossname == 'u'): loc = [1,0,0]
    elif(crossname=='v'): loc = [0,1,0]
    elif(crossname=='w'): loc = [0,0,1]
    else:                 loc = [0,0,0]

    locx = 'x' if loc[0] == 0 else 'xh'
    locy = 'y' if loc[1] == 0 else 'yh'
    locz = 'z' if loc[2] == 0 else 'zh'

    # create dimensions in netCDF file
    dim_x = crossfile.createDimension(locx,   nxsave)
    dim_y = crossfile.createDimension(locy,   nysave)
    dim_t = crossfile.createDimension('time', None)
    
    # create dimension variables
    var_t = crossfile.createVariable('time', sa, ('time',))
    var_x = crossfile.createVariable(locx,   sa, (locx,  ))
    var_y = crossfile.createVariable(locy,   sa, (locy,  ))
    var_t.units = "Seconds since start of experiment"
    
    # save the data
    var_x[:] = x[:nxsave] if locx=='x' else xh[:nxsave]
    var_y[:] = y[:nysave] if locy=='y' else yh[:nysave]
    
    var_s = crossfile.createVariable(crossname, sa, ('time', locx, locy,))
    
    stop = False 
    for t in range(niter):
        otime = int((starttime + t*sampletime) / 10**iotimeprec)
    
        try:
            fin = open("{0:}.xy.{1:07d}".format(crossname, otime), "rb")
        except:
            crossfile.sync()
            break
    
        print("Processing %8s, time=%7i"%(crossname, otime))

        var_t[t] = otime * 10**iotimeprec
    
        raw = fin.read(nx*ny*8)
        tmp = np.array(st.unpack('{0}{1}d'.format(en, nx*ny), raw))
        del(raw)
        s = tmp.reshape((ny, nx))
        del(tmp)
        var_s[t,:,:] = s[:nysave,:nxsave]
        del(s)
    
        fin.close()

    crossfile.close() 
