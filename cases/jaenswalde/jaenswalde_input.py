import matplotlib.pyplot as pl
import numpy as np
import netCDF4 as nc

float_type = "f8"

Mair = 28.9664
Mco2 = 44.01

# Get number of vertical levels and size from .ini file
with open('jaenswalde.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])
        if(line.split('=')[0]=='ysize'):
            ysize = float(line.split('=')[1])

dz = zsize / kmax

z = np.arange(0.5*dz, zsize, dz)
u = np.zeros(z.size)
v = np.zeros(z.size)
th = np.zeros(z.size)
co2 = np.zeros(z.size)

# well mixed profile with jump
th0 = 280.          # Bulk potential temperature
dth = 10.            # Potential temperature inversion
dthz = 300          # Inversion thickness
dthetadz = 0.006    # Free troposphere lapse rate
zi = 1200.          # Boundary layer depth
u0 = 5.             # Zonal wind speed
v0 = 0.             # Meridional wind speed
co20 = 400.         # CO2 concentration (ppm)

for k in range(kmax):
    if(z[k] <= zi - 0.5*dthz):
        th[k] = th0
    elif(z[k] <= zi + 0.5*dthz):
        th[k] = th0 + dth/dthz * (z[k]-(zi-0.5*dthz))
    else:
        th[k] = th0 + dth + dthetadz*(z[k]-(zi+0.5*dthz))

u[:] = u0
v[:] = v0
co2[:] = co20 * 1e-6 * Mco2/Mair    # ppmv -> kg kg-1

# Write input NetCDF file
nc_file = nc.Dataset("jaenswalde_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z" , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_u = nc_group_init.createVariable("u" , float_type, ("z"))
nc_v = nc_group_init.createVariable("v" , float_type, ("z"))
nc_th = nc_group_init.createVariable("th", float_type, ("z"))
nc_co2 = nc_group_init.createVariable("co2", float_type, ("z"))

nc_z [:] = z[:]
nc_u [:] = u[:]
nc_v [:] = v[:]
nc_th[:] = th[:]
nc_co2[:] = co2[:]

nc_file.close()

# Print .ini settings emissions:
# Coordinates of central cooling tower (m):
x0 = 1000
y0 = ysize/2.
z0 = 150

# x,y spacing towers:
dx = 290
dy = 120
ddx = 40

# Std-dev of plume widths:
sigma_x = 50
sigma_y = 50
sigma_z = 50

# Strength of plumes:
# CO2 emissions Jaenswalde, from:
# https://wwfeu.awsassets.panda.org/downloads/european_dirty_thirty_may_2007.pdf,
# which is probably not the best source.....
strength = 23.7e9 / 365.25 / 24. / 3600. / 9.
strength = np.round(strength, decimals=3)

# Generate lists with model input:
source_x0 = []
source_y0 = []
source_z0 = []

for j in range(-1,2):
    for i in range(-1,2):
        source_x0.append( x0 + i*dx + j*ddx )
        source_y0.append( y0 + j*dy )
        source_z0.append( z0 )

def to_cs_list(lst):
    """ From Python list to comma separated string """
    lst = [str(x) for x in lst]
    return ','.join(lst)

def constant_list(value, n):
    """ Create comma separated list with constant values """
    lst = n*[value]
    return to_cs_list(lst)

print('sourcelist={}'.format(constant_list('co2', 9)))

print('source_x0={}'.format(to_cs_list(source_x0)))
print('source_y0={}'.format(to_cs_list(source_y0)))
print('source_z0={}'.format(to_cs_list(source_z0)))

print('sigma_x={}'.format(constant_list(sigma_x, 9)))
print('sigma_y={}'.format(constant_list(sigma_y, 9)))
print('sigma_z={}'.format(constant_list(sigma_z, 9)))

print('strength={}'.format(constant_list(strength, 9)))

print('line_x={}'.format(constant_list(0, 9)))
print('line_y={}'.format(constant_list(0, 9)))
print('line_z={}'.format(constant_list(0, 9)))
