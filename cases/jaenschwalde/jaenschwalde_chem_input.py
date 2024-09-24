"""
Simplified Jaenschwalde setup.
Profiles are slightly idealised, from 04:00 UTC ERA5 data.
Wind is rotated to be perfectly westerly.
"""

import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
import netCDF4 as nc

pl.close('all'); pl.ion()

float_type = 'f8'

Mair = 28.9664

MCO2 = 44.01
MNO2 = 46.0
MNO = 30.0
MCO = 28.0

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
with open('jaenschwalde_chem.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])
        if(line.split('=')[0]=='ysize'):
            ysize = float(line.split('=')[1])

# Read TUV output table.
columns = ['time', 'sza', 'jo31d', 'jh2o2', 'jno2', 'jno3', 'jn2o5', 'jch2or', 'jch2om', 'jch3o2h']
tuv = pd.read_table(
        'jaensw_tuv_output.txt',
        sep='\\s+',
        skiprows=12,
        skipfooter=1,
        engine='python',
        names=columns,
        index_col='time')

# Case runs (hardcoded!) from 04:00 to 16:00 UTC!
tuv = tuv.loc[4:16]

# Convert to seconds, and subtract starting time.
tuv.index *= 3600
tuv.index -= tuv.index.values[0]

# Emissions (?)
emi_no = np.zeros(tuv.index.size)
emi_isop = np.zeros(tuv.index.size)

# Concentrations, for now constant with height......
species = {
    'co': 1.1e-7,
    'no': 3e-11,
    'no2': 4e-10,
    'hno3': 2e-9,
    'h2o2': 1.1e-9,
    'o3': 9e-8,
    'hcho': 5e-10,
    'ho2': 1.5e-11,
    'oh': 1e-13,
    'no3': 5e-14,
    'n2o5': 8e-14,
    'rooh': 3.6e-10,
    'c3h6': 3.3e-9,
    'ro2': 7e-12}

# Enable resolved plume rise:
sw_plume_rise = False

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
td1 = 11.5*3600
td2 = 14*3600

time = np.linspace(t0, t1, 32)
wthl = 0.17   * np.sin(np.pi * (time-t0-1800) / td1)
wqt  = 8.3e-5 * np.sin(np.pi * (time-t0) / td2)

#wthl[:] = 0.1
#wqt[:] = 0.1e-3

# Write input NetCDF file
def add_nc_var(name, dims, nc, data):
    if dims is None:
        var = nc.createVariable(name, np.float64)
    else:
        var = nc.createVariable(name, np.float64, dims)
    var[:] = data

nc_file = nc.Dataset('jaenschwalde_chem_input.nc', mode='w', datamodel='NETCDF4', clobber=True)

nc_file.createDimension('z', kmax)
add_nc_var('z', ('z'), nc_file, z)

# Atmospheric input.
nc_group_init = nc_file.createGroup('init');

add_nc_var('u', ('z'), nc_group_init, u)
add_nc_var('v', ('z'), nc_group_init, v)
add_nc_var('thl', ('z'), nc_group_init, thl)
add_nc_var('qt', ('z'), nc_group_init, qt)
add_nc_var('co2', ('z'), nc_group_init, co2)
add_nc_var('co2_inflow', ('z'), nc_group_init, co2)

nc_tdep = nc_file.createGroup('timedep');
nc_tdep.createDimension("time_surface", time.size)

add_nc_var('time_surface', ('time_surface'), nc_tdep, time-time[0])
add_nc_var('thl_sbot', ('time_surface'), nc_tdep, wthl)
add_nc_var('qt_sbot', ('time_surface'), nc_tdep, wqt)

# Chemistry input.
nc_chem = nc_file.createGroup('timedep_chem');
nc_chem.createDimension("time_chem", tuv.index.size)

add_nc_var("time_chem", ('time_chem'), nc_chem, tuv.index)
add_nc_var("jo31d", ('time_chem'), nc_chem, tuv.jo31d)
add_nc_var("jh2o2", ('time_chem'), nc_chem, tuv.jh2o2)
add_nc_var("jno2", ('time_chem'), nc_chem, tuv.jno2)
add_nc_var("jno3", ('time_chem'), nc_chem, tuv.jno3)
add_nc_var("jn2o5", ('time_chem'), nc_chem, tuv.jn2o5)
add_nc_var("jch2or", ('time_chem'), nc_chem, tuv.jch2or)
add_nc_var("jch2om", ('time_chem'), nc_chem, tuv.jch2om)
add_nc_var("jch3o2h", ('time_chem'), nc_chem, tuv.jch3o2h)
add_nc_var("emi_isop", ('time_chem'), nc_chem, emi_isop)
add_nc_var("emi_no", ('time_chem'), nc_chem, emi_no)

for name, value in species.items():
    profile = np.ones(kmax, dtype=np.float64)*value
    add_nc_var(name, ('z'), nc_group_init, profile)
    add_nc_var('{}_inflow'.format(name), ('z'), nc_group_init, profile)

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
strength_co2 = 732.5  / 9. / MCO2   # kmol(CO2) s-1
strength_no2 = 0.0289 / 9. / MNO2   # kmol(NO2) s-1
strength_no  = 0.359  / 9. / MNO    # kmol(NO) s-1
strength_co  = 0.223  / 9. / MCO    # kmol(CO) s-1

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
    raise Exception('Plume rise not yet implemented...')
else:
    source_list = \
            constant_list('co2', 9) + ',' +\
            constant_list('no2', 9) + ',' +\
            constant_list('no', 9) + ',' +\
            constant_list('co', 9)

    strength = \
            constant_list(strength_co2, 9) + ',' +\
            constant_list(strength_no2, 9) + ',' +\
            constant_list(strength_no, 9) + ',' +\
            constant_list(strength_co, 9)

    sw_vmr = constant_list('true', 4*9)
    repeat = 1

print('sourcelist={}'.format(source_list))

print('source_x0={}'.format(to_cs_list(source_x0, repeat*4)))
print('source_y0={}'.format(to_cs_list(source_y0, repeat*4)))
print('source_z0={}'.format(to_cs_list(source_z0, repeat*4)))

print('sigma_x={}'.format(constant_list(sigma_x, repeat*4*9)))
print('sigma_y={}'.format(constant_list(sigma_y, repeat*4*9)))
print('sigma_z={}'.format(constant_list(sigma_z, repeat*4*9)))

print('strength={}'.format(strength))
print('swvmr={}'.format(sw_vmr))

print('line_x={}'.format(constant_list(0, repeat*4*9)))
print('line_y={}'.format(constant_list(0, repeat*4*9)))
print('line_z={}'.format(constant_list(0, repeat*4*9)))
