#
#  MicroHH
#  Copyright (c) 2011-2023 Chiel van Heerwaarden
#  Copyright (c) 2011-2023 Thijs Heus
#  Copyright (c) 2014-2023 Bart van Stratum
#
#  This file is part of MicroHH
#
#  MicroHH is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MicroHH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
#

"""
Convert MicroHH input and output to a netCDF input file for the standalone version of the ray tracer (see /rte-rrtmgp-cpp/src_test)

How to use:
run 'python microhh_to_raytracing_input.py --name <simulation name> --time <time step to convert> --path <path to simulation files (defaults to "./")>'

Required input:
- <name>.ini
- <name>_input.nc
- <name>.default.0000000.nc
- grid.0000000
- T.<time> (binary 3D field of absolute temperature)
- qt.<time> (binary 3D field of specific humidity)

Optional input:
- ql.<time> (binary 3D field of liquid water specific humidity, if omitted, no liquid clouds are assumed to be present)
- qi.<time> (binary 3D field of ice specific humidity, if omitted, no ice clouds are assumed to be present)

To output these required 3D fields, add the following to <name>.ini
[dump]
swdump = 1
sampletime= <desired output time step>
dumplist = T,qt,ql,qi

"""

import microhh_tools as mht  # available in microhh/python directory
import numpy as np
import netCDF4 as nc
from scipy.interpolate import interp1d
import argparse

def solar_angles(lon, lat, day_of_year, year, seconds_since_midnight):
    #Based on: Paltridge, G. W. and Platt, C. M. R. (1976).
    #                 Radiative Processes in Meteorology and Climatology.
    #                 Elsevier, New York, 318 pp.

    if (year%4 == 0) and ((year%100 != 0) or (year%400 == 0)):
        days_per_year = 366
    else:
        days_per_year = 365

    doy = day_of_year-1
    radlat = lat * np.pi / 180.
    radlon = lon * np.pi / 180.

    doy_pi = 2.*np.pi*doy/days_per_year

    declination_angle = \
        0.006918 - 0.399912 * np.cos(doy_pi) + 0.070257 * np.sin(doy_pi)\
      - 0.006758 * np.cos(2*doy_pi) + 0.000907 * np.sin(2*doy_pi)\
      - 0.002697 * np.cos(3*doy_pi) + 0.00148  * np.sin(3*doy_pi)

    a1 = (1.00554 * doy -  6.28306) * np.pi/180
    a2 = (1.93946 * doy + 23.35089) * np.pi/180
    a3 = (7.67825 * np.sin(a1) + 10.09176 * np.sin(a2)) / 60

    hour_solar_time = (seconds_since_midnight/3600) - a3 + radlon * (180./np.pi/15.0)
    hour_angle = (hour_solar_time-12)*15.0*(np.pi/180)

    cos_zenith = np.sin(radlat) * np.sin(declination_angle)\
                    + np.cos(radlat) * np.cos(declination_angle) * np.cos(hour_angle)

    cos_elevation = np.cos(0.5*np.pi - np.arccos(cos_zenith))

    cos_azimuth = (\
          np.cos(radlat) * np.sin(declination_angle)\
        - np.sin(radlat) * np.cos(declination_angle) * np.cos(hour_angle) ) / cos_elevation

    azimuth = np.arccos(cos_azimuth) if hour_angle <= 0. else 2.*np.pi - np.arccos(cos_azimuth)

    return cos_zenith, azimuth

def calc_sun_distance_factor(day_of_year, seconds_since_midnight):
    # Based on: An Introduction to Atmospheric Radiation, Liou, Eq. 2.2.9.
    an = [1.000110, 0.034221, 0.000719]
    bn = [0,        0.001280, 0.000077]
    frac_doy = day_of_year + seconds_since_midnight / 86400.
    t = 2. * np.pi*(frac_doy - 1.)/ 365.

    factor = 0.
    for n in range(3):
        factor += an[n]*np.cos(n*t) + bn[n]*np.sin(n*t);

    return factor

def read_if_exists(var, t, dims):
    try:
        return np.fromfile("{}.{:07d}".format(var,t),TF).reshape(dims)
    except:
        return np.zeros(dims)

# functions to compute the relative humidity
Rd = 287.04
Rv = 461.5
ep = Rd/Rv
T0 = 273.15

def esat_liq(T):
    x = np.maximum(-75., T-T0)
    return 611.21*np.exp(17.502*x / (240.97+x))

def esat_ice(T):
    x = np.maximum(-100., T-T0)
    return 611.15*np.exp(22.452*x / (272.55+x))

def qsat_liq(p, T):
    return ep*esat_liq(T)/(p-(1.-ep)*esat_liq(T))

def qsat_ice(p, T):
    return ep*esat_ice(T)/(p-(1.-ep)*esat_ice(T))

def water_fraction(T):
    return np.maximum(0., np.minimum((T - 233.15) / (T0 - 233.15), 1.))

def qsat(p, T):
    alpha = water_fraction(T)
    return alpha*qsat_liq(p, T) + (1.-alpha)*qsat_ice(p, T)


parser = argparse.ArgumentParser()
parser.add_argument("-n","--name", type=str, help="simulating name, that is, name of the .ini file")
parser.add_argument("-t","--time", type=int, help="simulation time step to convert")
parser.add_argument("-p","--path", type=str, help="Path to simulation files, defaults to current directory",default="./")
args = parser.parse_args()

# some constants
xm_air = 28.97
xm_h2o = 18.01528
eps = xm_h2o / xm_air
TF = np.float32
g = 9.81

# size of null-collision grid:
ng_x = 48
ng_y = 48
ng_z = 32

# input: name of case, time step
name = args.name
time = args.time
path = args.path

### open all necessary files
# read namelist
nl = mht.Read_namelist(path+"{}.ini".format(name))

# precision of timestep in microhh output
iotimeprec = nl['time']['iotimeprec'] if 'iotimeprec' in nl['time'] else 0
iotime = int(round(time / 10**iotimeprec))

# convert simulation to local time
lon = nl['grid']['lon']
lat = nl['grid']['lat']
dt_utc = nl['time']['datetime_utc']
year = int(dt_utc[:4])
day_of_year = int((np.datetime64(dt_utc[:10]) - np.datetime64("{}-01-01".format(year))) / np.timedelta64(1, 'D')) + 1
seconds_since_midnight = int(dt_utc[11:13])*3600 + int(dt_utc[14:16])*60 + int(dt_utc[17:19]) + time

# compute solar angles
mu0, azi = solar_angles(lon, lat, day_of_year, year, seconds_since_midnight)

# scale top-of-atmosphere irradiance based on distance to sun at time of year
tsi_scaling = calc_sun_distance_factor(day_of_year, seconds_since_midnight)

# read stats
nc_stat = nc.Dataset(path+"{}.default.0000000.nc".format(name))
t_idx = np.where(nc_stat['time'][:] == time)[0][0]

# read input
nc_inp = nc.Dataset(path+"{}_input.nc".format(name))

# pressure and base state density profiles
play = nc_stat['thermo']['phydro'][t_idx]
plev = nc_stat['thermo']['phydroh'][t_idx]
rho = nc_stat['thermo']['rho'][t_idx]

### Read data
# switches
swtimedep = nl['radiation']['swtimedep_background'] if 'swtimedep_background' in nl['radiation'] else 0
gaslist = nl['radiation']['timedeplist_gas'] if 'timedeplist_gas' in nl['radiation'] else []
if 'aerosol' in nl.groups:
    swaerosol = nl['aerosol']['swaerosol'] if 'swaerosol' in nl['aerosol'] else 0
    swtimedep_aerosol = nl['aerosol']['swtimedep'] if 'swtimedep' in nl['aerosol'] else 0
else:
    swaerosol = 0
    swtimedep_aerosol = 0

# dimensions and grid
itot = nl['grid']['itot']
jtot = nl['grid']['jtot']
ktot = nl['grid']['ktot']
dims = (ktot, jtot, itot)

grid = mht.Read_grid(itot, jtot, ktot,filename=path+"grid.0000000")

dz = np.diff(grid.dim['zh'][:])

zlay = grid.dim['z']
zlev = grid.dim['zh']

# read temperature, humidity
qt = np.fromfile(path+"qt.{:07d}".format(iotime),TF).reshape(dims)
tlay = np.fromfile(path+"T.{:07d}".format(iotime),TF).reshape(dims)

# convert qt from kg/kg to vmr
h2o = qt / (ep - ep*qt)
h2o = np.maximum(0, h2o)

# cloud properties and effective radius
ql = read_if_exists(path+"ql", iotime, dims)
lwp = ql * (dz*rho)[:,np.newaxis,np.newaxis] # kg/m2

qi = read_if_exists(path+"qi", iotime, dims)
iwp = qi * (dz*rho)[:,np.newaxis,np.newaxis] # kg/m2

nc0 = nl['micro']['Nc0']
ftpnr_w = (4./3) * np.pi * nc0 * 1e3
ftpnr_i = (4./3) * np.pi * 1e5 * 7e2
sig_fac = np.exp(np.log(1.34)*np.log(1.34))
rel = np.where(lwp>0, 1e6 * sig_fac * (lwp / dz[:,np.newaxis,np.newaxis] / ftpnr_w)**(1./3), 0)
rei = np.where(iwp>0, 1e6 * sig_fac * (iwp / dz[:,np.newaxis,np.newaxis] / ftpnr_i)**(1./3), 0)

rel = np.maximum(2.5, np.minimum(rel, 21.5))
rei = np.maximum(10., np.minimum(rei, 180.))

lwp *= 1e3 # g/m2
iwp *= 1e3 # g/m2

nz = ktot
grid_z = grid.dim['z']
grid_zh = grid.dim['zh']

# read bg profile
if swtimedep:
    f_h2o_bg = interp1d(nc_inp['timedep']['time_rad'], nc_inp['timedep']['h2o_bg'], axis=0)
    h2o_bg = f_h2o_bg(time)
    f_zlay_bg = interp1d(nc_inp['timedep']['time_rad'], nc_inp['timedep']['z_lay'], axis=0)
    zlay_bg = f_zlay_bg(time)
    f_zlev_bg = interp1d(nc_inp['timedep']['time_rad'], nc_inp['timedep']['z_lev'], axis=0)
    zlev_bg = f_zlev_bg(time)
    f_play_bg = interp1d(nc_inp['timedep']['time_rad'], nc_inp['timedep']['p_lay'], axis=0)
    play_bg = f_play_bg(time)
    f_plev_bg = interp1d(nc_inp['timedep']['time_rad'], nc_inp['timedep']['p_lev'], axis=0)
    plev_bg = f_plev_bg(time)
    f_tlay_bg = interp1d(nc_inp['timedep']['time_rad'], nc_inp['timedep']['t_lay'], axis=0)
    tlay_bg = f_tlay_bg(time)
    f_tlev_bg = interp1d(nc_inp['timedep']['time_rad'], nc_inp['timedep']['t_lev'], axis=0)
    tlev_bg = f_tlev_bg(time)
else:
    h2o_bg = nc_inp['radiation']['h2o']
    zlay_bg = nc_inp['radiation']['z_lay']
    zlev_bg = nc_inp['radiation']['z_lev']
    play_bg = nc_inp['radiation']['p_lay']
    plev_bg = nc_inp['radiation']['p_lev']
    tlay_bg = nc_inp['radiation']['t_lay']
    tlev_bg = nc_inp['radiation']['t_lev']


# find lowest height in bg profile that is heigher than domain top
z_tod = grid.dim['zh'][-1]
zmin_idx = np.where(zlay_bg[:] > z_tod)[0][0]

zlay = np.append(zlay, zlay_bg[zmin_idx:])
zlev = np.append(zlev, zlev_bg[zmin_idx+1:])
zlev[nz] = nz*(grid_z[1]-grid_z[0])

# patch pressure profiles
play = np.append(play, play_bg[zmin_idx:])
plev = np.append(plev, plev_bg[zmin_idx+1:])

### Writing output
# create netcdf file
nc_out = nc.Dataset("rte_rrtmgp_input.nc", "w")

## create dimensions
nc_out.createDimension("band_sw", 14)
nc_out.createDimension("band_lw", 16)
nc_out.createDimension("lay", len(play))
nc_out.createDimension("lev", len(plev))
nc_out.createDimension("x", itot)
nc_out.createDimension("y", jtot)
nc_out.createDimension("z", nz)
nc_out.createDimension("xh", itot+1)
nc_out.createDimension("yh", jtot+1)
nc_out.createDimension("zh", nz+1)

# write raytracing grids
nc_x = nc_out.createVariable("x", TF, ("x",))
nc_x[:] = grid.dim['x'][:]
nc_y = nc_out.createVariable("y", TF, ("y",))
nc_y[:] = grid.dim['y'][:]
nc_z = nc_out.createVariable("z", TF, ("z",))
nc_z[:] = grid_z

grid_xh = np.append(grid.dim['xh'][:], grid.dim['xh'][-1] +  (grid.dim['xh'][1] -  grid.dim['xh'][0]))
grid_yh = np.append(grid.dim['yh'][:], grid.dim['yh'][-1] +  (grid.dim['yh'][1] -  grid.dim['yh'][0]))

nc_xh = nc_out.createVariable("xh", TF, ("xh",))
nc_xh[:] = grid_xh
nc_yh = nc_out.createVariable("yh", TF, ("yh",))
nc_yh[:] = grid_yh
nc_zh = nc_out.createVariable("zh", TF, ("zh",))
nc_zh[:] = grid_zh

nc_zlay = nc_out.createVariable("zlay", TF, ("lay"))
nc_zlev = nc_out.createVariable("z_lev", TF, ("lev"))
nc_zlay[:] = zlay
nc_zlev[:] = zlev

# write pressures
nc_play = nc_out.createVariable("p_lay", TF, ("lay","y","x"))
nc_play[:] = np.tile(play.reshape(len(play),1,1), (1, jtot, itot))
nc_plev = nc_out.createVariable("p_lev", TF, ("lev","y","x"))
nc_plev[:] = np.tile(plev.reshape(len(plev),1,1), (1, jtot, itot))

# write density profile
nc_rho = nc_out.createVariable("rho", TF, ("z",))
nc_rho[:] = rho

# remaining 3D variables
nc_h2o = nc_out.createVariable("vmr_h2o", TF, ("lay","y","x"))
nc_h2o[:] = np.append(h2o[:], np.tile(h2o_bg[zmin_idx:][:,None,None], (1, jtot, itot)), axis=0)
nc_tlay = nc_out.createVariable("t_lay", TF, ("lay","y","x"))
nc_tlay[:] = np.append(tlay[:], np.tile(tlay_bg[zmin_idx:][:,None,None], (1, jtot, itot)), axis=0)

# We do not bother about t_lev yet  because the ray tracer is shortwave-only, but we do need to supply it in the netcdf
nc_tlev = nc_out.createVariable("t_lev", TF, ("lev","y","x"))
nc_tlev[:] = 0 

# Liquid water path
nc_lwp = nc_out.createVariable("lwp" , TF, ("lay","y","x"))
nc_lwp[:] = 0
nc_lwp[:ktot] = lwp

# Liquid water effective radius
nc_rel = nc_out.createVariable("rel" , TF, ("lay","y","x"))
nc_rel[:] = 0
nc_rel[:ktot] = rel

# Ice water path
nc_iwp = nc_out.createVariable("iwp" , TF, ("lay","y","x"))
nc_iwp[:] = 0
nc_iwp[:ktot] = iwp

# Ice effective radius
nc_rei = nc_out.createVariable("rei" , TF, ("lay","y","x"))
nc_rei[:] = 0
nc_rei[:ktot] = rei

# surface properties
nc_alb_dir = nc_out.createVariable("sfc_alb_dir", TF, ("y","x","band_sw"))
nc_alb_dir[:] = nl['radiation']['sfc_alb_dir']
nc_alb_dif = nc_out.createVariable("sfc_alb_dif", TF, ("y","x","band_sw"))
nc_alb_dif[:] = nl['radiation']['sfc_alb_dif']
nc_emis = nc_out.createVariable("emis_sfc", TF, ("y","x","band_lw"))
nc_emis[:] = nl['radiation']['emis_sfc']
nc_tsfc = nc_out.createVariable("t_sfc", TF, ("y","x"))
nc_tsfc[:] = 0 # don't bother about longwave for now

# solar angles
nc_mu = nc_out.createVariable("mu0", TF, ("y","x"))
nc_mu[:] = mu0
nc_az = nc_out.createVariable("azi", TF, ("y","x"))
nc_az[:] = azi

# Scaling top-of-atmosphere irradiance
nc_ts = nc_out.createVariable("tsi_scaling", TF)
nc_ts[:] = tsi_scaling

# add gasses
possible_gases = ["co2", "o3", "n2o", "co", "ch4", "o2", "n2", "ccl4", "cfc11", "cfc12", "cfc22",
                  "hfc143a", "hfc125", "hfc23", "hfc32", "hfc134a", "cf4", "no2"]
for gas in possible_gases:
    if gas in nc_inp['radiation'].variables.keys():
        if len(nc_inp['radiation'][gas].dimensions) != 0 and len(nc_inp['init'][gas].dimensions) != 0:
            # in the domain
            if gas in gaslist:
                f_gas = interp1d(nc_inp['timedep']['time_rad'], nc_inp['timedep'][gas], axis=0)
                gas_prof = f_gas(time)
            else:
                gas_prof = nc_inp['init'][gas]
            # above the domain
            if gas in gaslist and swtimedep:
                f_gas_bg = interp1d(nc_inp['timedep']['time_rad'], nc_inp['timedep'][gas + '_bg'], axis=0)
                gas_bg = f_gas_bg(time)
            else:
                gas_bg = nc_inp['radiation'][gas]

            nc_out.createVariable("vmr_" + gas, TF, ("lay", "y", "x"))
            nc_out["vmr_" + gas][:] = np.tile(np.append(gas_prof[:], gas_bg[zmin_idx:])[:, None, None], (1, jtot, itot))

        elif len(nc_inp['radiation'][gas].dimensions) == 0 and len(nc_inp['init'][gas].dimensions) == 0:
            nc_out.createVariable("vmr_" + gas, TF)
            nc_out["vmr_" + gas][:] = nc_inp['radiation'][gas][:]

        else:
            print("not supported input for gas: " + gas)

# add aerosols + relative humidity
if swaerosol:
    aerosol_species = ["aermr01", "aermr02", "aermr03", "aermr04", "aermr05", "aermr06",
                         "aermr07", "aermr08", "aermr09", "aermr10", "aermr11"]

    # aerosol mixing ratios
    for aerosol in aerosol_species:
        # in the domain
        if swtimedep_aerosol:
            f_aermr = interp1d(nc_inp['timedep']['time_rad'], nc_inp['timedep'][aerosol], axis=0)
            aermr = np.maximum(f_aermr(time), 0)
        else:
            aermr = nc_inp['init'][aerosol]
        # above the domain
        if swtimedep_aerosol and swtimedep:
            f_aermr_bg = interp1d(nc_inp['timedep']['time_rad'], nc_inp['timedep'][aerosol + '_bg'], axis=0)
            aermr_bg = np.maximum(f_aermr_bg(time), 0)
        else:
            aermr_bg = nc_inp['radiation'][aerosol]

        nc_out.createVariable(aerosol, TF, ("lay", "y", "x"))
        nc_out[aerosol][:] = np.tile(np.append(aermr[:], aermr_bg[zmin_idx:])[:, None, None], (1, jtot, itot))

    # relative humidity
    q = nc_h2o[:] * eps / (1 + nc_h2o[:] * eps)
    qs = qsat(nc_play[:, :, :], nc_tlay[:])
    rh = np.maximum(np.minimum(q / qs, 1), 0)
    nc_out.createVariable("rh", TF, ("lay", "y", "x"))
    nc_out['rh'][:] = rh

# size of null-collision grid
nc_ng_x = nc_out.createVariable("ngrid_x", TF)
nc_ng_x[:] = ng_x
nc_ng_y = nc_out.createVariable("ngrid_y", TF)
nc_ng_y[:] = ng_y
nc_ng_z = nc_out.createVariable("ngrid_z", TF)
nc_ng_z[:] = ng_z

nc_out.close()
