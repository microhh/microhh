import numpy as np
import netCDF4 as nc

# float_type = "f8"
float_type = "f4"

# Get number of vertical levels and size from .ini file
with open('drycblles_fire.ini') as f:
    for line in f:
        if(line.split('=')[0]=='itot'):
            itot = int(line.split('=')[1])
        if(line.split('=')[0]=='jtot'):
            jtot = int(line.split('=')[1])
        if(line.split('=')[0]=='ktot'):
            ktot = int(line.split('=')[1])
        if(line.split('=')[0]=='xsize'):
            xsize = float(line.split('=')[1])
        if(line.split('=')[0]=='ysize'):
            ysize = float(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / ktot

dthetadz = 0.003

# set the height
z  = np.linspace(0.5*dz, zsize-0.5*dz, ktot)
u  = np.zeros(np.size(z))
v  = np.zeros(np.size(z))
th = np.zeros(np.size(z))

# linearly stratified profile
for k in range(ktot):
    th  [k] = 300. + dthetadz*z[k]

"""
# well mixed profile with jump
h    = 1000.
dth  = 10.
dthz = 100.

for k in range(ktot):
    if(z[k] <= h - 0.5*dthz):
        th[k] = 300.
    elif(z[k] <= h + 0.5*dthz):
        th[k] = 300. + dth/dthz * (z[k]-(h-0.5*dthz))
    else:
        th[k] = 300. + dth + dthetadz*(z[k]-(h+0.5*dthz))
    thls[k] = 2.*(z[k]/zsize - 0.5) / 3600.
    wls [k] = -0.01*(z[k]/zsize)
"""

nc_file = nc.Dataset("drycblles_fire_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", ktot)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_u  = nc_group_init.createVariable("u" , float_type, ("z"))
nc_v  = nc_group_init.createVariable("v" , float_type, ("z"))
nc_th = nc_group_init.createVariable("th", float_type, ("z"))

nc_z [:] = z [:]
nc_u [:] = u [:]
nc_v [:] = v [:]
nc_th[:] = th[:]

nc_file.close()


# Create 2D files for surface fluxes of temperature. Mimics a fire plume to test 2D sbot.
dx = xsize / itot
dy = ysize / jtot
x = np.arange(dx/2, xsize, dx)
y = np.arange(dy/2, ysize, dy)
L = 200.

for n in range(0, 10801, 1800):
    th_bot = 0.1*np.ones((jtot, itot), dtype=float_type)
    th_bot[:, :] += 10.*max(0., (n - 3600.)/ (10800. - 3600.)) * np.exp(-(x[None, :] - xsize/2)**2 / L**2) * np.exp(-(y[:, None] - ysize/2)**2 / L**2)
    print(n, th_bot.max())
    th_bot.tofile('th_bot_in.{:07d}'.format(n))

