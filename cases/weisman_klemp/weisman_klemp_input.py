import numpy as np
import netCDF4 as nc

float_type = 'f8'

# ***** Parameters for WK ********
# Tropopause parameters
T_tr     = 213.         # Temperature  (K)
theta_tr = 343.         # Potential temperature (K)
z_tr     = 12000.       # Height (m)

# Surface Parameters
theta_0  = 300.         # Potential temperature (K)

# Shear Parameters
z_sh    =  3000.        # Heigth of maximum shear (m)

# Constants (from constants.h)
cp       = 1005.        # Specific heat of air at constant pressure [J kg-1 K-1]
Rd       = 287.04       # Gas constant for dry air [J K-1 kg-1]
Rv       = 461.5        # Gas constant for water vapor [J K-1 kg-1]
p0       = 1.e5         # Reference pressure [Pa]
T0       = 273.15       # Freezing / melting temperature [K]
ep       = Rd/Rv
rdcp     = Rd/cp
cprd     = cp/Rd
grav     = 9.81

# Variations used 
qv0      = 0.014       # Vapor Mixing ratio in BL (maximum value) (kg/kg) [ 11 / 14 / 16 ]
Us       = 25.         # Shear velocity (m/s) [ 5 15 25 35 45 ]

# My functions
# From thermo_moist_functions.h
#def esat_liq(T):
#  Tc = T - T0
#  return 611.21 * numpy.exp(17.502 * Tc / (240.97 + Tc))

#def qsat_liq(p,T):
#  return ep*esat_liq(T)/(p-(1.0-ep)*esat_liq(T))

#Chiel
# From thermo_moist_functions.h
def esat_liq(T):
    c0 = 0.6105851e+03; c1 = 0.4440316e+02; c2 = 0.1430341e+01; c3 = 0.2641412e-01
    c4 = 0.2995057e-03; c5 = 0.2031998e-05; c6 = 0.6936113e-08; c7 = 0.2564861e-11
    c8 = -.3704404e-13
    x  = max(-80.,T-273.15)
    return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

def qsat_liq(p, T):
    return ep*esat_liq(T)/(p-(1.-ep)*esat_liq(T))

def qv_rh(p, T, rh):
    return ep*esat_liq(T)*rh/(p-(1.-ep)*esat_liq(T)*rh)


# Get number of vertical levels and size from .ini file
with open('weisman_klemp.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

# set the height
z     = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl   = np.zeros(z.size)
qt    = np.zeros(z.size)
u     = np.zeros(z.size)
v     = np.zeros(z.size)

# Temporary Arrays
theta = np.zeros(z.size)
rh    = np.zeros(z.size)

# Search for the tropopause Ask Chiel
if ( z[kmax-1] < z_tr ):
    raise SystemExit('Domain is too small in z to fit the tropopause')

k = 0
while z[k] < z_tr:
    k_tr = k
    k = k + 1

# Pressure in tropopause
p_tr = p0 * pow( T_tr/ theta_tr, cprd)



# Calculate theta and rh arrays and velocity
for k in range(k_tr):
    thl[k] = theta_0 + (theta_tr - theta_0) *pow( z[k] / z_tr ,5/4)
    rh[k]  = max(1. - 0.75 * pow( z[k] / z_tr ,5/4), 0.95)

for k in range(k_tr,kmax):
    thl[k] = theta_tr * np.exp( grav/cp/T_tr * (z[k] - z_tr) )
    rh[k] = 0.25

for k in range(kmax):
    u[k] = Us*np.tanh( z[k]/ z_sh )

# Calculate values in above tropopause (constant temperature)
const_trop = -grav/Rd/T_tr
for k in range(k_tr,kmax):
    p_loc =  p0 * pow( T_tr/ thl[k], rdcp)
    qt[k] = rh[k] * qv_rh(p_loc,T_tr,rh[k])


# Calculate values below tropopause
cpres = grav*pow(p0,rdcp)*dz/cp
qfg = 0
p_up = p_tr
for k in range(k_tr-1, -2, -1):
    # First guess no humidity
    pfg   = pow(pow(p_up, rdcp) + cpres/thl[k]/(1. + ep*qfg), cprd)
    Ttemp = thl[k] * pow(pfg/p0, rdcp)
    qfg   = min(qv0, qv_rh(pfg, Ttemp,rh[k]))

    # Second guess with humidity
    p_up  = pow(pow(p_up,rdcp) + cpres/thl[k]/(1. + ep*qfg), cprd)
    Ttemp = thl[k] * pow(p_up/p0, rdcp)
    if ( k > -1):
        qt[k] = min(qv0, qv_rh(p_up, Ttemp,rh[k]))	
        #print(z[k], p_up/1.e5, Ttemp, qt[k])

print("Please set surface pressure in ini file to:", p_up)    

# write the data to a file
nc_file = nc.Dataset("weisman_klemp_input.nc", mode="w", datamodel="NETCDF4", clobber=True)
nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_thl = nc_group_init.createVariable("thl", float_type, ("z"))
nc_qt  = nc_group_init.createVariable("qt" , float_type, ("z"))
nc_u   = nc_group_init.createVariable("u"  , float_type, ("z"))
nc_v   = nc_group_init.createVariable("v"  , float_type, ("z"))

nc_z    [:] = z    [:]
nc_thl  [:] = thl  [:]
nc_qt   [:] = qt   [:]
nc_u    [:] = u    [:]
nc_v    [:] = v    [:]

nc_file.close()

