import matplotlib.pyplot as pl
import numpy as np
import netCDF4 as nc

pl.close('all'); pl.ion()

float_type = 'f8'

Mair = 28.9664
Mco2 = 44.01
T0 = 273.15
Rd = 287.04
Rv = 461.5
ep = Rd/Rv

def calc_esat(T):
    a_ab = 611.21; b_ab = 18.678; c_ab = 234.5; d_ab = 257.14
    return a_ab * np.exp((b_ab - ((T-T0) / c_ab)) * ((T-T0) / (d_ab + (T-T0))))

def calc_qsat(T, p):
    esat = calc_esat(T)
    return ep*esat / (p-(1.-ep)*esat)

# Get number of vertical levels and size from .ini file
with open('jaenschwalde.ini') as f:
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
qt = np.zeros(z.size)
co2 = np.zeros(z.size)

# well mixed profile with jump
th0 = 290.           # Bulk potential temperature
dth = 2.             # Potential temperature inversion
dthdz = 0.006        # Free troposphere lapse rate

qt0 = 5e-3          # Bulk specific humidity
dqt = -2e-3          # Specific humidity temperature inversion
dqtdz = -0.003e-3    # Free troposphere lapse rate

dzi = 200            # Inversion thickness
zi = 1000.           # Boundary layer depth
u0 = 3               # Zonal wind speed
v0 = 0.              # Meridional wind speed
co20 = 400e-6        # CO2 concentration (vmr, -)

sw_plume_rise = True

for k in range(kmax):
    if(z[k] <= zi - 0.5*dzi):
        th[k] = th0
        qt[k] = qt0
    elif(z[k] <= zi + 0.5*dzi):
        th[k] = th0 + dth/dzi * (z[k]-(zi-0.5*dzi))
        qt[k] = qt0 + dqt/dzi * (z[k]-(zi-0.5*dzi))
    else:
        th[k] = th0 + dth + dthdz*(z[k]-(zi+0.5*dzi))
        qt[k] = qt0 + dqt + dqtdz*(z[k]-(zi+0.5*dzi))

"""
pl.figure()
pl.subplot(121)
pl.plot(th, z)
pl.subplot(122)
pl.plot(qt*1e3, z)
"""

u[:] = u0
v[:] = v0
co2[:] = co20

# Write input NetCDF file
nc_file = nc.Dataset('jaenschwalde_input.nc', mode='w', datamodel='NETCDF4', clobber=True)

nc_file.createDimension('z', kmax)
nc_z = nc_file.createVariable('z' , float_type, ('z'))

nc_group_init = nc_file.createGroup('init');
nc_u = nc_group_init.createVariable('u' , float_type, ('z'))
nc_v = nc_group_init.createVariable('v' , float_type, ('z'))
nc_th = nc_group_init.createVariable('thl', float_type, ('z'))
nc_qt = nc_group_init.createVariable('qt', float_type, ('z'))
nc_co2 = nc_group_init.createVariable('co2', float_type, ('z'))
nc_co2_inflow = nc_group_init.createVariable('co2_inflow', float_type, ('z'))

nc_z [:] = z[:]
nc_u [:] = u[:]
nc_v [:] = v[:]
nc_th[:] = th[:]
nc_qt[:] = qt[:]
nc_co2[:] = co2[:]
nc_co2_inflow[:] = co2[:]

nc_file.close()

# Print .ini settings emissions:
# Coordinates of central cooling tower (m):
x0 = 1000
y0 = ysize/2.

# Std-dev of plume widths:
sigma_x = 25
sigma_y = 25

if sw_plume_rise:
    z0 = 120
    sigma_z = 25
else:
    z0 = 300
    sigma_z = 100

# x,y spacing towers:
dx = 290
dy = 120
ddx = 40

# Strength of plumes, from the CoCO2 simulation protocol:
strength_co2 = np.round(732.5 / 9. / Mco2, decimals=2)   # kmol(CO2) s-1

# Emission of heat and moisture. Numbers are from:
# Effective pollutant emission heights for atmospheric transport modelling based on real-world information
# Pregger and Friedrich, 2009, 10.1016/j.envpol.2008.09.027
Tp = 50+T0    # Flue gass temperature (K)
Mp = 790      # Volume-flux (m3 s-1)

# This is not very accurate...:
pp = 1e5
rhop = pp/(Rd*Tp)
rhp = 1.
qp = rhp * calc_qsat(Tp, pp)

strength_q = np.round(Mp * rhop * qp, decimals=2)
strength_T = np.round(Mp * rhop * Tp, decimals=2)

# Generate lists with model input:
source_x0 = []
source_y0 = []
source_z0 = []

for j in range(-1,2):
    for i in range(-1,2):
        source_x0.append( x0 + i*dx + j*ddx )
        source_y0.append( y0 + j*dy )
        source_z0.append( z0 )

def to_cs_list(lst, n=1):
    """ From Python list to comma separated string """
    out = ''
    for i in range(n):
        out += ','.join([str(x) for x in lst])
        if i < n-1:
            out += ','
    return out

def constant_list(value, n):
    """ Create comma separated list with constant values """
    lst = n*[value]
    return to_cs_list(lst)

if sw_plume_rise:
    source_list = constant_list('co2', 9) + ',' + \
                  constant_list('thl', 9) + ',' + \
                  constant_list('qt',  9)
    strength = constant_list(strength_co2, 9) + ',' + \
               constant_list(strength_T, 9) + ',' + \
               constant_list(strength_q, 9)
    sw_vmr = constant_list(True, 9) + ',' + \
             constant_list(False, 9) + ',' + \
             constant_list(False,  9)
    repeat = 3
else:
    source_list = constant_list('co2', 9)
    strength = constant_list(strength_co2, 9)
    sw_vmr = constant_list(True, 9)
    repeat = 1

print('sourcelist={}'.format(source_list))

print('source_x0={}'.format(to_cs_list(source_x0, repeat)))
print('source_y0={}'.format(to_cs_list(source_y0, repeat)))
print('source_z0={}'.format(to_cs_list(source_z0, repeat)))

print('sigma_x={}'.format(constant_list(sigma_x, repeat*9)))
print('sigma_y={}'.format(constant_list(sigma_y, repeat*9)))
print('sigma_z={}'.format(constant_list(sigma_z, repeat*9)))

print('strength={}'.format(strength))
print('swvmr={}'.format(sw_vmr))

print('line_x={}'.format(constant_list(0, repeat*9)))
print('line_y={}'.format(constant_list(0, repeat*9)))
print('line_z={}'.format(constant_list(0, repeat*9)))
