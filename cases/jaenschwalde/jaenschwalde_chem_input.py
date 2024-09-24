"""
Simplified Jaenschwalde setup.
Profiles are slightly idealised, from 04:00 UTC ERA5 data.
Wind is rotated to be perfectly westerly.
"""

import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
import netCDF4 as nc

import microhh_tools as mht

pl.close('all')


"""
Help functions.
"""
def calc_esat(T):
    """
    Calculate saturation vapor pressure (Pa).
    """
    a_ab = 611.21; b_ab = 18.678; c_ab = 234.5; d_ab = 257.14
    return a_ab * np.exp((b_ab - ((T-T0) / c_ab)) * ((T-T0) / (d_ab + (T-T0))))


def calc_qsat(T, p):
    """
    Calculate saturation specific humidity (kg kg-1)
    """
    esat = calc_esat(T)
    return ep*esat / (p-(1.-ep)*esat)


class Emissions:
    """
    Help class to gather emission settings.
    """
    def __init__(self):

        self.source_list = []
        self.strength = []
        self.sw_vmr = []

        self.x0 = []
        self.y0 = []
        self.z0 = []

        self.sigma_x = []
        self.sigma_y = []
        self.sigma_z = []

        self.line_x = []
        self.line_y = []
        self.line_z = []

    def add(
            self, name, strength, sw_vmr,
            x0, y0, z0,
            sigma_x, sigma_y, sigma_z,
            line_x=0, line_y=0, line_z=0):

        self.source_list.append(name)
        self.strength.append(strength)
        self.sw_vmr.append(sw_vmr)

        self.x0.append(x0)
        self.y0.append(y0)
        self.z0.append(z0)

        self.sigma_x.append(sigma_x)
        self.sigma_y.append(sigma_y)
        self.sigma_z.append(sigma_z)

        self.line_x.append(line_x)
        self.line_y.append(line_y)
        self.line_z.append(line_z)


"""
Settings.
"""
float_type = 'f8'

xsize = 6400
ysize = 3200
zsize = 4000

itot = 64
jtot = 32
ktot = 96

start_hour = 4
end_hour = 12

# Enable resolved plume rise:
sw_plume_rise = True


"""
Constants.
"""
Mair = 28.9664

MCO2 = 44.01
MNO2 = 46.0
MNO = 30.0
MCO = 28.0

T0 = 273.15
Rd = 287.04
Rv = 461.5
ep = Rd/Rv


# Read .ini file for case settings.
ini = mht.Read_namelist('jaenschwalde_chem.ini.base')

"""
Create case input.
"""
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

# NOTE: `.loc` is value based, not in index.
tuv = tuv.loc[start_hour:end_hour]

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

# Vertical grid LES
dz = zsize / ktot
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

v = np.zeros(ktot)
co2 = np.zeros(ktot)

# Surface fluxes, again idealised from ERA5.
t0 = start_hour*3600
t1 = end_hour*3600
td1 = 11.5*3600
td2 = 14*3600

time = np.linspace(t0, t1, 32)
wthl = 0.17   * np.sin(np.pi * (time-t0-1800) / td1)
wqt  = 8.3e-5 * np.sin(np.pi * (time-t0) / td2)


"""
Write input NetCDF file.
"""
def add_nc_var(name, dims, nc, data):
    if dims is None:
        var = nc.createVariable(name, np.float64)
    else:
        var = nc.createVariable(name, np.float64, dims)
    var[:] = data

nc_file = nc.Dataset('jaenschwalde_chem_input.nc', mode='w', datamodel='NETCDF4', clobber=True)

nc_file.createDimension('z', ktot)
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
    profile = np.ones(ktot, dtype=np.float64)*value
    add_nc_var(name, ('z'), nc_group_init, profile)
    add_nc_var('{}_inflow'.format(name), ('z'), nc_group_init, profile)

nc_file.close()


"""
Define emissions.
"""
# Coordinates of central cooling tower (m):
x0 = 1000
y0 = ysize/2.

# Std-dev of plume widths:
sigma_x = 25
sigma_y = 25

if sw_plume_rise:
    # Emissions from tower height.
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

# Emission input model:
emi = Emissions()

for j in range(-1,2):
    for i in range(-1,2):
        x = x0 + i*dx + j*ddx
        y = y0 + j*dy
        z = z0
            
        emi.add('co2', strength_co2, True, x, y, z, sigma_x, sigma_y, sigma_z)
        emi.add('no2', strength_no2, True, x, y, z, sigma_x, sigma_y, sigma_z)
        emi.add('no',  strength_no,  True, x, y, z, sigma_x, sigma_y, sigma_z)
        emi.add('co',  strength_co,  True, x, y, z, sigma_x, sigma_y, sigma_z)

        if sw_plume_rise:
            emi.add('thl', strength_T, False, x, y, z, sigma_x, sigma_y, sigma_z)
            emi.add('qt',  strength_q, False, x, y, z, sigma_x, sigma_y, sigma_z)



"""
Add settings to .ini file.
"""
ini['grid']['itot'] = itot
ini['grid']['jtot'] = jtot
ini['grid']['ktot'] = ktot

ini['grid']['xsize'] = xsize
ini['grid']['ysize'] = ysize
ini['grid']['zsize'] = zsize

ini['time']['endtime'] = (end_hour - start_hour) * 3600

ini['cross']['xz'] = ysize/2

ini['source']['sourcelist'] = emi.source_list

ini['source']['source_x0'] = emi.x0
ini['source']['source_y0'] = emi.y0
ini['source']['source_z0'] = emi.z0

ini['source']['sigma_x'] = emi.sigma_x
ini['source']['sigma_y'] = emi.sigma_y
ini['source']['sigma_z'] = emi.sigma_z

ini['source']['strength'] = emi.strength
ini['source']['swvmr'] = emi.sw_vmr

ini['source']['line_x'] = emi.line_x
ini['source']['line_y'] = emi.line_y
ini['source']['line_z'] = emi.line_z

ini.save('jaenschwalde_chem.ini', allow_overwrite=True)
