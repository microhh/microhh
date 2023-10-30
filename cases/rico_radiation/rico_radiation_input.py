import numpy as np
import netCDF4 as nc
import microhh_tools as mht
from scipy.interpolate import interp1d

float_type = 'f8'
eps = 287.04 / 461.5 #(Rd/Rv)
p00 = 1e5
Rd = 287.04
cp = 1005.

# If true, use Monte Carlo ray tracer for 3D radiative transfer. If false, use the two-stream solver.
use_rt = False
# If true: solar angles vary with time-of-day (by default, simulations then start at 00:00 LT).
diurnal_cycle = False

# Read ini file
ini = mht.Read_namelist('rico_radiation.ini')

# Get number of vertical levels and size
kmax  = ini['grid']['ktot']
zsize = ini['grid']['zsize']

dz = zsize / kmax

# Get surface values
thl_0 = ini['boundary']['sbot[thl]']
q_0 = ini['boundary']['sbot[qt]']
p_0 = ini['thermo']['pbot']

T_0 = thl_0*(p00/p_0)**(-Rd/cp)

# Spatially constant gas concentrations
vmr_gases_const = {'co2' : 348.e-6,
                   'ch4' : 1650.e-9,
                   'n2o' : 306.e-9,
                   'n2':  0.7808,
                   'o2':  0.2095,
                   'cfc11': 0,
                   'cfc12': 0,
                   'cfc22': 0,
                   'ccl4': 0}

# function for ozone profile
def calc_o3(z,p):
    g1 = 3.6478
    g2 = 0.83209
    g3 = 11.3515
    p_hpa = p/100.
    o3 = g1 * p_hpa**g2 * np.exp(-p_hpa/g3) * 1e-6
    return o3

# Create input file
nc_file = nc.Dataset("rico_radiation_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, ("z"))


#######################################
# Create init group: initial conditions and large-scale tendencies

# set the height
z     = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl   = np.zeros(z.size)
qt    = np.zeros(z.size)
u     = np.zeros(z.size)
ugeo  = np.zeros(z.size)
v     = np.zeros(z.size)
vgeo  = np.zeros(z.size)
wls   = np.zeros(z.size)
thlls = np.zeros(z.size)
qtls  = np.zeros(z.size)

for k in range(kmax):

    # Liquid water potential temperature: same in GCSS and SS08
    if(z[k] < 740.):
        thl[k] = 297.9
    else:
        thl[k] = 297.9 + (317.0 - 297.9)/(4000. - 740.) * (z[k] - 740.)

    if(z[k] < 740.):
        qt[k] = 16.0 + (13.8 - 16.0) / 740. * z[k]
    elif(z[k] < 3260.):
        qt[k] = 13.8 + (2.4 - 13.8) / (3260. - 740.) * (z[k] - 740.)
    else:
        qt[k] = 2.4 + (1.8 - 2.4)/(4000. - 3260.) * (z[k] - 3260.)

    # Subsidence
    if(z[k] < 2260):
        wls[k] = -0.005 * (z[k] / 2260.)
    else:
        wls[k] = -0.005

    # U and V component wind
    u[k] = -9.9 + 2.0e-3 * z[k]
    v[k] = -3.8
    ugeo[k] = u[k]
    vgeo[k] = v[k]

    # No large-scale tendency (comparable to irad-zero experiment in Seifert et al., 2015: http://dx.doi.org/10.1002/2015MS000489). To include advective cooling, set to -2.5 / 86400 K/s
    thlls[k] = 0

    # Advective tendency qt
    if(z[k] < 2980):
        qtls[k] = -1.0 / 86400. + (1.3456/ 86400.) * z[k] / 2980.
    else:
        qtls[k] = 4e-6

# normalize profiles to SI
qt  /= 1000.
qtls/= 1000.

# convert qt to h2o_vmr:
h2o = qt / (eps - eps * qt)

nc_group_init = nc_file.createGroup("init");
nc_thl   = nc_group_init.createVariable("thl"   , float_type, ("z"))
nc_qt    = nc_group_init.createVariable("qt"    , float_type, ("z"))
nc_h2o   = nc_group_init.createVariable("h2o"   , float_type, ("z"))
nc_o3    = nc_group_init.createVariable("o3"    , float_type, ("z"))
nc_u     = nc_group_init.createVariable("u"     , float_type, ("z"))
nc_ugeo  = nc_group_init.createVariable("u_geo" , float_type, ("z"))
nc_v     = nc_group_init.createVariable("v"     , float_type, ("z"))
nc_vgeo  = nc_group_init.createVariable("v_geo" , float_type, ("z"))
nc_wls   = nc_group_init.createVariable("w_ls"  , float_type, ("z"))
nc_thlls = nc_group_init.createVariable("thl_ls", float_type, ("z"))
nc_qtls  = nc_group_init.createVariable("qt_ls" , float_type, ("z"))

# domain-o3 is ommitted here. We compute this later on by interpolating o3 from the radiation profile to the domain grid
nc_z    [:] = z    [:]
nc_thl  [:] = thl  [:]
nc_qt   [:] = qt   [:]
nc_h2o  [:] = h2o  [:]
nc_u    [:] = u    [:]
nc_ugeo [:] = ugeo [:]
nc_v    [:] = v    [:]
nc_vgeo [:] = vgeo [:]
nc_wls  [:] = wls  [:]
nc_thlls[:] = thlls[:]
nc_qtls [:] = qtls [:]

for gas in vmr_gases_const.keys():
    nc_group_init.createVariable(gas, float_type)
    nc_group_init[gas][:] = vmr_gases_const[gas]

##########################################
# Create radiation group: background profiles and constant gas concentrations

# For the radiation background profile, we follow the RCEMIP specifications of Wing et al. (2018), but with T0 and q_0 of RICO
def calc_p_q_T_thl(z):
    z_q1 = 4.0e3
    z_q2 = 7.5e3
    z_t = 15.e3
    q_t = 1.e-14

    q = q_0 * np.exp(-z  /z_q1) * np.exp(-(z  /z_q2)**2)

    # CvH hack to remove moisture jump.
    q_tb = q_0 * np.exp(-z_t/z_q1) * np.exp(-(z_t/z_q2)**2)
    q -= q_tb + q_t

    i_above_zt = np.where(z > z_t)
    q[i_above_zt] = q_t

    gamma = 6.7e-3
    Tv_0 = (1. + 0.608*q_0)*T_0
    Tv = Tv_0 - gamma*z
    Tv_t = Tv_0 - gamma*z_t
    Tv[i_above_zt] = Tv_t
    T = Tv / (1. + 0.608*q)

    g = 9.79764

    p = p_0 * (Tv / Tv_0)**(g/(Rd*gamma))

    p_tmp = p_0 * (Tv_t/Tv_0)**(g/(Rd*gamma)) \
          * np.exp( -( (g*(z-z_t)) / (Rd*Tv_t) ) )

    p[i_above_zt] = p_tmp[i_above_zt]

    thl = T*(p00/p)**(Rd/cp)

    return p, q, T, thl

# Vertical layers and levels in background profile
z_rad_top = 70.e3
dz_rad = 500.
z_rad  = np.arange(dz_rad/2, z_rad_top, dz_rad)
zh_rad = np.arange(   0, z_rad_top-dz_rad/2, dz_rad)
zh_rad = np.append(zh_rad, z_rad_top)

p_lay, q_t, T_lay, _ = calc_p_q_T_thl( z_rad)
p_lev,   _, T_lev, _ = calc_p_q_T_thl(zh_rad)

# convert qt to h2o_vmr:
h2o = q_t/(eps - eps*q_t)
o3 = calc_o3(z_rad, p_lay)

nc_group_rad = nc_file.createGroup("radiation")

nc_group_rad.createDimension("lay", p_lay.size)
nc_group_rad.createDimension("lev", p_lev.size)

nc_rad_z_lay = nc_group_rad.createVariable("z_lay", float_type, ("lay"))
nc_rad_z_lev = nc_group_rad.createVariable("z_lev", float_type, ("lev"))
nc_rad_z_lay[:] = z_rad [:]
nc_rad_z_lev[:] = zh_rad[:]

nc_rad_p_lay = nc_group_rad.createVariable("p_lay", float_type, ("lay"))
nc_rad_p_lev = nc_group_rad.createVariable("p_lev", float_type, ("lev"))
nc_rad_p_lay[:] = p_lay[:]
nc_rad_p_lev[:] = p_lev[:]

nc_rad_T_lay = nc_group_rad.createVariable("t_lay", float_type, ("lay"))
nc_rad_T_lev = nc_group_rad.createVariable("t_lev", float_type, ("lev"))
nc_rad_T_lay[:] = T_lay[:]
nc_rad_T_lev[:] = T_lev[:]

nc_rad_O3  = nc_group_rad.createVariable("o3" , float_type, ("lay"))
nc_rad_H2O = nc_group_rad.createVariable("h2o", float_type, ("lay"))

nc_rad_O3 [:] = o3[:]
nc_rad_H2O[:] = h2o[:]

for gas in vmr_gases_const.keys():
    nc_group_rad.createVariable(gas, float_type)
    nc_group_rad[gas][:] = vmr_gases_const[gas]

# interpolate radiation-o3 to domain-o3
zzh_rad = np.append(z_rad, zh_rad)
pph_rad = np.append(p_lay, p_lev)
o3_rad = calc_o3(zzh_rad, pph_rad)
idx_sorted = zzh_rad.argsort()
o3_rad = o3_rad[idx_sorted]
zzh_rad = zzh_rad[idx_sorted]
print(zzh_rad,o3_rad)
f_o3 = interp1d(zzh_rad, o3_rad, fill_value="extrapolate")
domain_o3 = f_o3(z)
nc_o3[:] = domain_o3

##########################################
# Remaining radiation settings in ini file
if use_rt and not diurnal_cycle:
    raise Exception("Unfortunately, we do not support fixed azmiuth angles yet, hence ray tracing with fixed solar angles is not possible")

# two-stream or ray tracer
ini['radiation']['swradiation'] = 'rrtmgp' + '_rt'*use_rt

if diurnal_cycle:
    ini['radiation']['swfixedsza'] = False
    ini['radiation']['tsi_scaling'] = -999
    ini['radiation']['sza'] = -999

else:
    ini['radiation']['swfixedsza'] = True
    # tsi_scaling and fixed sza result in an insolation representing the daily average insolation during the period of RICO campaign on which the simulatin is based.
    ini['radiation']['tsi_scaling'] = 0.29081
    ini['radiation']['sza'] = np.deg2rad(42.05)

ini.save("rico_radiation.ini", allow_overwrite=True)

nc_file.close()

mht.copy_radfiles()
