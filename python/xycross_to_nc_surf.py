import numpy   as np
import struct  as st
import netCDF4 as nc4

# Settings -------
variables  = ['rr_bot']
nx         = 32
ny         = 32
nz         = 100
starttime  = 0.
endtime    = 10800.
sampletime = 60
iotimeprec = 0
nxsave     = nx
nysave     = ny
endian     = 'little'
dtype      = 'float'    # Input data type
savetype   = 'float'    # Output data type
# End settings ---

# Set the correct string for the endianness
if (endian == 'little'):
    en = '<'
elif (endian == 'big'):
    en = '>'
else:
    raise RuntimeError("Endianness has to be little or big")

# Set the correct size for the input data type
if (dtype == 'double'):
    var_size = 8
    var_type = 'd'
elif (savetype == 'float'):
    var_size = 4
    var_type = 'f'
else:
    raise RuntimeError("The dtype has to be float or double")

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
fin = open("grid.{:07d}".format(0), "rb")
raw = fin.read(nx*var_size)
x   = np.array(st.unpack('{0}{1}{2}'.format(en, nx, var_type), raw))
raw = fin.read(nx*var_size)
xh  = np.array(st.unpack('{0}{1}{2}'.format(en, nx, var_type), raw))
raw = fin.read(ny*var_size)
y   = np.array(st.unpack('{0}{1}{2}'.format(en, ny, var_type), raw))
raw = fin.read(ny*var_size)
yh  = np.array(st.unpack('{0}{1}{2}'.format(en, ny, var_type), raw))
raw = fin.read(nz*var_size)
z   = np.array(st.unpack('{0}{1}{2}'.format(en, nz, var_type), raw))
raw = fin.read(nz*var_size)
zh  = np.array(st.unpack('{0}{1}{2}'.format(en, nz, var_type), raw))
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

    var_s = crossfile.createVariable(crossname, sa, ('time', locy, locx,))

    stop = False
    for t in range(niter):
        otime = int((starttime + t*sampletime) / 10**iotimeprec)

        file_name = "{0:}.xy.{1:07d}".format(crossname, otime)
        try:
            fin = open(file_name, "rb")
        except:
            print('Cant find file {}'.format(file_name))
            crossfile.sync()
            break

        print("Processing %8s, time=%7i"%(crossname, otime))

        var_t[t] = otime * 10**iotimeprec

        raw = fin.read(nx*ny*var_size)
        tmp = np.array(st.unpack('{0}{1}{2}'.format(en, nx*ny, var_type), raw))
        del(raw)
        s = tmp.reshape((ny, nx))
        del(tmp)
        var_s[t,:,:] = s[:nysave,:nxsave]
        del(s)

        fin.close()

    crossfile.close()
