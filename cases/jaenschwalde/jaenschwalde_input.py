"""
Simplified Jaenschwalde setup.
Profiles are slightly idealised, from 04:00 UTC ERA5 data.
Wind is rotated to be perfectly westerly.
"""

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

# Enable resolved plume rise:
sw_plume_rise = True

# Vertical grid LES
dz = zsize / kmax
z = np.arange(0.5*dz, zsize, dz)

v_thl = np.array([285.7, 291.9, 293, 297.4, 307])
z_thl = np.array([0, 400, 2000, 2500, 5000])
thl = np.interp(z, z_thl, v_thl)

z_qt = np.array([0, 400, 2000, 2500, 5000])
v_qt = np.array([6.2, 4.93, 3.61, 1, 0.3])/1000
qt = np.interp(z, z_qt, v_qt)

z_u = np.array([0, 270, 3000, 5000])
v_u = np.array([2.3, 8.5, 0.6, 5.7])
u = np.interp(z, z_u, v_u)

v = np.zeros(kmax)
co2 = np.zeros(kmax)

# Surface fluxes, again idealised from ERA5.
t0 = 4*3600
t1 = 16*3600
td1 = 12*3600
td2 = 14*3600

time = np.linspace(t0, t1, 32)
wthl = 0.17   * np.sin(np.pi * (time-t0) / td1)
wqt  = 8.3e-5 * np.sin(np.pi * (time-t0) / td2)

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
nc_th[:] = thl[:]
nc_qt[:] = qt[:]
nc_co2[:] = co2[:]
nc_co2_inflow[:] = co2[:]

nc_group_tdep = nc_file.createGroup('timedep');
nc_group_tdep.createDimension("time_surface", time.size)
nc_time_surface = nc_group_tdep.createVariable("time_surface", float_type, ("time_surface"))
nc_thl_sbot = nc_group_tdep.createVariable("thl_sbot", float_type, ("time_surface"))
nc_qt_sbot = nc_group_tdep.createVariable("qt_sbot" , float_type, ("time_surface"))

nc_time_surface[:] = time
nc_thl_sbot[:] = wthl
nc_qt_sbot[:] = wqt

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
    # The heights and sigma are from Dominik's .csv profiles, curve fitted with Python.
    z0 = 299.68    # of 599.69 for high
    sigma_z = 122.37

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
    sw_vmr = constant_list('true', 9) + ',' + \
             constant_list('false', 9) + ',' + \
             constant_list('false',  9)
    repeat = 3
else:
    source_list = constant_list('co2', 9)
    strength = constant_list(strength_co2, 9)
    sw_vmr = constant_list('true', 9)
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
