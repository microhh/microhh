import numpy as np
import netCDF4 as nc
import microhh_tools as mht
from ls2d import grid as ls2d_grid

float_type = 'f8'

# switch microphysics scheme
sw_micro = 'sb06'       # options: 'sb06', 'nsw6', 0
sw_ice=0                # switch for including/excluding ice (ony for swmicro = 'sb06' and swmicro=0)
sw_thldeep = 0          # switch between thlD and thlE (available for all micro options)

# with the default resolution of 500m it is not recommended to change the advection scheme and/or leave out the flux limiter
sw_advec = '2i5'
fluxlim_qx = True

# ***** Parameters for WK ********
# Tropopause parameters
T_tr = 213.  # Temperature  (K)
theta_tr = 343.  # Potential temperature (K)
z_tr = 12000.  # Height (m)

# Surface Parameters
theta_0 = 300.  # Potential temperature (K)

# Constants (from constants.h)
cp = 1005.  # Specific heat of air at constant pressure [J kg-1 K-1]
Rd = 287.04  # Gas constant for dry air [J K-1 kg-1]
Rv = 461.5  # Gas constant for water vapor [J K-1 kg-1]
p0 = 1.e5  # Reference pressure [Pa]
T0 = 273.15  # Freezing / melting temperature [K]
ep = Rd / Rv
rdcp = Rd / cp
cprd = cp / Rd
grav = 9.81

# Variations used
qv0 = 0.013  # Vapor Mixing ratio in BL (maximum value) (kg/kg) [ 11 / 14 / 16 ]
Us = 10  # Shear velocity (m/s) [ 10 17.5 25 ] [5, 15, 25, 35, 45)
z_sh = 2500.  # Height of maximum shear (m) [2500 5000]


# From thermo_moist_functions.h
def esat_liq(T):
    x = min(max(-100., T - 273.15), 50)
    return 611.21*np.exp(17.502*x / (240.97+x))


def qv_rh(p, T, rh):
    return ep * esat_liq(T) * rh / (p - (1. - ep) * esat_liq(T) * rh)


# vertical grid
dz0 = 50
kmax = 128
heights = [0, 200, 500, 1000, 2000, 4000, 8000, 14000, 1e12]
factors = [1.16704, 1.04472, 1.02662, 1.01704, 1.01136, 1.00807, 1.00680, 1.00676]
grid = ls2d_grid.Grid_stretched_manual(kmax=kmax, dz0=dz0, heights=heights, factors=factors)
z = grid.z
dz = grid.dz[1:]

thl = np.zeros(z.size)
qt = np.zeros(z.size)
u = np.zeros(z.size)
v = np.zeros(z.size)

# Temporary Arrays
theta = np.zeros(z.size)
rh = np.zeros(z.size)

# Search for the tropopause Ask Chiel
if z[kmax - 1] < z_tr:
    raise SystemExit('Domain is too small in z to fit the tropopause')

k = 0
while z[k] < z_tr:
    k_tr = k
    k = k + 1

# Pressure in tropopause
p_tr = p0 * pow(T_tr / theta_tr, cprd)

# Calculate theta and rh arrays and velocity
for k in range(k_tr):
    thl[k] = theta_0 + (theta_tr - theta_0) * pow(z[k] / z_tr, 5 / 4)
    rh[k] = 1. - 0.75 * (z[k] / z_tr) ** 1.25
    rh[k] = min(rh[k], 0.95)        # limit max rh to be in agreement with ICON

for k in range(k_tr, kmax):
    thl[k] = theta_tr * np.exp(grav / cp / T_tr * (z[k] - z_tr))
    rh[k] = 0.1

# wind profile following weismann et al. 1996
for k in range(kmax):
    u[k] = np.minimum(Us, Us/z_sh * z[k])

# Calculate values in above tropopause (constant temperature)
const_trop = -grav / Rd / T_tr
for k in range(k_tr, kmax):
    p_loc = p0 * pow(T_tr / thl[k], cprd)
    qt[k] = qv_rh(p_loc, T_tr, rh[k])

# Calculate values below tropopause
# cpres = grav * pow(p0, rdcp) * dz / cp
qfg = 0
p_up = p_tr
for k in range(k_tr - 1, -2, -1):
    # First guess no humidity
    cpres = grav * pow(p0, rdcp) * dz[k] / cp if k > -1 else grav * pow(p0, rdcp) * dz[0] / cp
    pfg = pow(pow(p_up, rdcp) + cpres / thl[k] / (1. + ep * qfg), cprd)
    Ttemp = thl[k] * pow(pfg / p0, rdcp) if k > -1 else theta_0 * pow(pfg / p0, rdcp)
    qfg = min(qv0, qv_rh(pfg, Ttemp, rh[k]))

    # Second guess with humidity
    p_up = pow(pow(p_up, rdcp) + cpres / thl[k] / (1. + ep * qfg), cprd)
    Ttemp = thl[k] * pow(p_up / p0, rdcp) if k > -1 else theta_0 * pow(p_up / p0, rdcp)
    if k > -1:
        qt[k] = min(qv0, qv_rh(p_up, Ttemp, rh[k]))


# create .ini file
ini = mht.Read_namelist('squall_line.ini.base')
ini['grid']['ktot'] = kmax
ini['grid']['zsize'] = grid.zsize
ini['buffer']['zstart'] = grid.zsize*5/6.

ini['thermo']['pbot'] = p_up

ini['micro']['swmicro'] = sw_micro

if sw_micro == 'nsw6':
    prog_species = ['qr', 'qs', 'qg']
    diag_species = ['ql', 'qi']
    bonus_cross = ['thl', 'qt', 'u', 'v', 'w', 'qlqi','T']
    precip_rates = ["rr_bot", "rs_bot", "rg_bot"]
    ini['thermo']['swsatadjust_ql'] = 1
    ini['thermo']['swsatadjust_qi'] = 1

elif sw_micro == 'sb06' and sw_ice == 0:
    prog_species = ['qr', 'nr']
    diag_species = ['ql']
    bonus_cross = ['thl', 'qt', 'u', 'v', 'w','T']
    precip_rates = ["rr_bot"]
    ini['thermo']['swsatadjust_ql'] = 1
    ini['thermo']['swsatadjust_qi'] = 0
    ini['micro']['swice'] = sw_ice

elif sw_micro == 'sb06' and sw_ice == 1:
    prog_species = ["qr","nr","qi","ni","qs","ns","qg","ng","qh","nh","ina"]
    diag_species = ['ql']
    bonus_cross = ['thl', 'qt', 'u', 'v', 'w','T']
    precip_rates = ["ri_bot","rs_bot","rg_bot","rh_bot","rr_bot"]
    ini['thermo']['swsatadjust_ql'] = 1
    ini['thermo']['swsatadjust_qi'] = 0
    ini['micro']['swice'] = sw_ice

elif sw_micro == 0:
    prog_species = []
    bonus_cross = ['thl', 'qt', 'u', 'v', 'w','T']
    precip_rates = []
    ini['thermo']['swsatadjust_ql'] = 1
    ini['thermo']['swsatadjust_qi'] = sw_ice
    if sw_ice:
        diag_species = ['ql', 'qi']
    else:
        diag_species = ['ql']

else:
    raise Exception('Unknown sw_micro!')

ini['cross']['crosslist'] = ['thl'] + precip_rates
ini['limiter']['limitlist'] = ['qt'] + prog_species
ini['limiter']['cliplist'] = ['qt'] + prog_species
ini['dump']['dumplist'] = prog_species + diag_species + bonus_cross

ini['advec']['swadvec'] = sw_advec
if fluxlim_qx:
    ini['advec']['fluxlimit_list'] = ['qt'] + prog_species

ini['thermo']['swthldeep'] = sw_thldeep

ini.save('squall_line.ini', allow_overwrite=True)

# write the data to a file
nc_file = nc.Dataset("squall_line_input.nc", mode="w", datamodel="NETCDF4", clobber=True)
nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, "z")

nc_group_init = nc_file.createGroup("init")
nc_thl = nc_group_init.createVariable("thl", float_type, ("z"))
nc_qt = nc_group_init.createVariable("qt", float_type, ("z"))
nc_u = nc_group_init.createVariable("u", float_type, ("z"))
nc_v = nc_group_init.createVariable("v", float_type, ("z"))

nc_z[:] = z[:]
nc_thl[:] = thl[:]
nc_qt[:] = qt[:]
nc_u[:] = u[:]
nc_v[:] = v[:]

nc_file.close()

