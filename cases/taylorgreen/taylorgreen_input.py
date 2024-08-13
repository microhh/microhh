import netCDF4 as nc
import numpy

float_type = "f8"
# float_type = "f4"

# set the height
# Get number of vertical levels and size from .ini file
with open('taylorgreen.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

# define the variables
dz = zsize / kmax
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)

# write data to a file
nc_file = nc.Dataset("taylorgreen_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", kmax)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))
nc_z[:] = z[:]

nc_group_init = nc_file.createGroup("init");

nc_file.close()
