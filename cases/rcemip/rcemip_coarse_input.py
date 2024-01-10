import numpy as np
import netCDF4 as nc

float_type = "f8"

#T_0 = 295.
#q_0 = 0.01200 # for 295 K SST.

T_0 = 300.
q_0 = 0.01864 # for 300 K SST.

#T_0 = 305.
#q_0 = 0.02400 # for 305 K SST.

eps = 18.01528 / 28.9647 # molar mass water / molar mass air

def q_sat(T, p):
    Tc = T - 273.15

    # Arden-Buck equation.
    e_sat = 611.21 * np.exp(17.502 * Tc / (240.97 + Tc))
    Rd, Rv = 287.04, 461.5
    return Rd/Rv * e_sat / (p - (1. - Rd/Rv)*e_sat)

def calc_p_q_T_thl_o3(z):
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
    Rd = 287.04
    cp = 1005.
    p0 = 101480.

    print("q_sat at T_0 = ", q_sat(T_0, p0))

    p = p0 * (Tv / Tv_0)**(g/(Rd*gamma))
    
    p_tmp = p0 * (Tv_t/Tv_0)**(g/(Rd*gamma)) \
          * np.exp( -( (g*(z-z_t)) / (Rd*Tv_t) ) )
    
    p[i_above_zt] = p_tmp[i_above_zt]

    p00 = 1e5
    thl = T*(p00/p)**(Rd/cp)

    g1 = 3.6478
    g2 = 0.83209
    g3 = 11.3515
    p_hpa = p/100.
    o3 = g1 * p_hpa**g2 * np.exp(-p_hpa/g3) * 1e-6

    return p, q, T, thl, o3

nc_file = nc.Dataset("rcemip_coarse_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

### RADIATION INIT ###
# Radiation profiles.
z_top = 100.e3
dz = 500.
z  = np.arange(dz/2, z_top, dz)
zh = np.arange(   0, z_top-dz/2, dz)
zh = np.append(zh, z_top)

p_lay, q, T_lay, _, o3 = calc_p_q_T_thl_o3( z)
p_lev, _, T_lev, _,  _ = calc_p_q_T_thl_o3(zh)

h2o = q / (eps - eps*q)

co2 =  348.e-6
ch4 = 1650.e-9
n2o =  306.e-9
n2 = 0.7808
o2 = 0.2095

g1 = 3.6478
g2 = 0.83209
g3 = 11.3515
p_hpa = p_lay/100.
o3 = g1 * p_hpa**g2 * np.exp(-p_hpa/g3) * 1e-6

nc_group_rad = nc_file.createGroup("radiation")

nc_group_rad.createDimension("lay", p_lay.size)
nc_group_rad.createDimension("lev", p_lev.size)

nc_z_lay = nc_group_rad.createVariable("z_lay", float_type, ("lay"))
nc_z_lev = nc_group_rad.createVariable("z_lev", float_type, ("lev"))
nc_z_lay[:] = z [:]
nc_z_lev[:] = zh[:]

nc_p_lay = nc_group_rad.createVariable("p_lay", float_type, ("lay"))
nc_p_lev = nc_group_rad.createVariable("p_lev", float_type, ("lev"))
nc_p_lay[:] = p_lay[:]
nc_p_lev[:] = p_lev[:]

nc_T_lay = nc_group_rad.createVariable("t_lay", float_type, ("lay"))
nc_T_lev = nc_group_rad.createVariable("t_lev", float_type, ("lev"))
nc_T_lay[:] = T_lay[:]
nc_T_lev[:] = T_lev[:]

nc_CO2 = nc_group_rad.createVariable("co2", float_type)
nc_CH4 = nc_group_rad.createVariable("ch4", float_type)
nc_N2O = nc_group_rad.createVariable("n2o", float_type)
nc_O3  = nc_group_rad.createVariable("o3" , float_type, ("lay"))
nc_H2O = nc_group_rad.createVariable("h2o", float_type, ("lay"))
nc_N2  = nc_group_rad.createVariable("n2" , float_type)
nc_O2  = nc_group_rad.createVariable("o2" , float_type)

nc_CFC11 = nc_group_rad.createVariable("cfc11", float_type)
nc_CFC12 = nc_group_rad.createVariable("cfc12", float_type)
nc_CFC22 = nc_group_rad.createVariable("cfc22", float_type)
nc_CCL4  = nc_group_rad.createVariable("ccl4" , float_type)

nc_CO2[:] = co2
nc_CH4[:] = ch4
nc_N2O[:] = n2o
nc_O3 [:] = o3 [:]
nc_H2O[:] = h2o[:]
nc_N2 [:] = n2
nc_O2 [:] = o2

nc_CFC11[:] = 0.
nc_CFC12[:] = 0.
nc_CFC22[:] = 0.
nc_CCL4 [:] = 0.

### INITIAL PROFILES ###
# Get number of vertical levels and size from .ini file
with open('rcemip_coarse.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

# set the height
# dz = zsize / kmax
# z = np.linspace(0.5*dz, zsize-0.5*dz, kmax)

# Official RCEMIP coarse vertical resolution.
z1 = np.array([37., 112., 194., 288., 395., 520., 667., 843., 1062., 1331., 1664., 2055., 2505., 3000.]);
z2 = np.arange(3500., 33001., 500.)
z = np.concatenate((z1, z2))
z = z[:-2]
zh = 0.5*(z[:-1] + z[1:])
zh = np.append(0., zh)
zh = np.append(zh, zsize)

if (z.size != kmax):
    raise RuntimeError("kmax does not match the RCEMIP profile")

_, qt, _, thl, o3 = calc_p_q_T_thl_o3(z)

h2o = qt / (eps - eps*qt)

nc_file.createDimension("z", kmax)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))
nc_z[:] = z[:]

# Initial profiles.
nc_group_init = nc_file.createGroup("init");
nc_thl = nc_group_init.createVariable("thl", float_type, ("z"))
nc_qt  = nc_group_init.createVariable("qt" , float_type, ("z"))

nc_CO2 = nc_group_init.createVariable("co2", float_type)
nc_CH4 = nc_group_init.createVariable("ch4", float_type)
nc_N2O = nc_group_init.createVariable("n2o", float_type)
nc_O3  = nc_group_init.createVariable("o3" , float_type, ("z"))
nc_H2O = nc_group_init.createVariable("h2o", float_type, ("z"))
nc_N2  = nc_group_init.createVariable("n2" , float_type)
nc_O2  = nc_group_init.createVariable("o2" , float_type)

nc_CFC11 = nc_group_init.createVariable("cfc11", float_type)
nc_CFC12 = nc_group_init.createVariable("cfc12", float_type)
nc_CFC22 = nc_group_init.createVariable("cfc22", float_type)
nc_CCL4  = nc_group_init.createVariable("ccl4" , float_type)

nc_thl[:] = thl[:]
nc_qt [:] = qt [:]

nc_CO2[:] = co2
nc_CH4[:] = ch4
nc_N2O[:] = n2o
nc_O3 [:] = o3
nc_H2O[:] = h2o[:]
nc_N2 [:] = n2
nc_O2 [:] = o2

nc_CFC11[:] = 0.
nc_CFC12[:] = 0.
nc_CFC22[:] = 0.
nc_CCL4 [:] = 0.

nc_file.close()
