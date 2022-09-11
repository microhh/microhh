import microhh_tools as mht  # available in microhh/python directory
import numpy as np
import netCDF4 as nc
from scipy.interpolate import interp1d
import argparse

def solar_angles(lon, lat, day_of_year, year, seconds_since_midnight):

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

    return cos_zenith, azimuth #np.rad2deg(np.arccos(cos_zenith)),np.rad2deg(azimuth))

def calc_sun_distance_factor(day_of_year, seconds_since_midnight):
    an = [1.000110, 0.034221, 0.000719]
    bn = [0,        0.001280, 0.000077]
    day_of_year = 228.3
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

parser = argparse.ArgumentParser()
parser.add_argument("--time", type=int)
parser.add_argument("--path", type=str, default="./")
parser.add_argument("--name", type=str, default="")
parser.add_argument("--merge_bg", action='store_true')
parser.add_argument("--ztilt", type=float, default=-1)
parser.add_argument("--regrid", action='store_true')
parser.add_argument("--zmax", type=float, default=4000)
parser.add_argument("--dz", type=float, default=25)
args = parser.parse_args()

### example usage:
# First: make sure that grid.0000000 is in folder
# For ray tracing:
# python mhh_to_raytracer.py --regrid --zmax 2500 --dz 25
# For tilted column:
# python mhh_to_raytracer --merge_bg --ztilt 4000
# constants
ep = 0.622
TF = np.float32
g = 9.81

#ngrid:
ng_x = 48
ng_y = 48
ng_z = 32

# input: name of case, time step
name = args.name
time = args.time
path = args.path
merge_bg = args.merge_bg
z_tilt = args.ztilt
f_regrid = args.regrid
regrid_zmax = args.zmax
regrid_dz = args.dz

if merge_bg & f_regrid:
    raise "nope, can't do both"

### open all necessary files
# read namelist
nl = mht.Read_namelist(path+"{}.ini".format(name))

lon = nl['radiation']['lon']
lat = nl['radiation']['lat']
dt_utc = nl['time']['datetime_utc']
year = int(dt_utc[:4])
day_of_year = int((np.datetime64(dt_utc[:10]) - np.datetime64("{}-01-01".format(year))) / np.timedelta64(1, 'D')) + 1
seconds_since_midnight = int(dt_utc[11:13])*3600 + int(dt_utc[14:16])*60 + int(dt_utc[17:19]) + time

mu0, azi = solar_angles(lon, lat, day_of_year, year, seconds_since_midnight)
tsi_scaling = calc_sun_distance_factor(day_of_year, seconds_since_midnight)

# read stats
nc_stat = nc.Dataset(path+"{}.default.0000000.nc".format(name))
t_idx = np.where(nc_stat['time'][:] == time)[0][0]

# read input
nc_inp = nc.Dataset(path+"{}_input.nc".format(name))

# pressure and base state density profiles
play = nc_stat['thermo']['phydro'][t_idx]
plev = nc_stat['thermo']['phydroh'][t_idx]
rhoh = nc_stat['thermo']['rhoref'][:]

# height diff

### Read data
# dimensions and grid
itot = nl['grid']['itot']
jtot = nl['grid']['jtot']
ktot = nl['grid']['ktot']
dims = (ktot, jtot, itot)
grid = mht.Read_grid(itot, jtot, ktot)
dz = np.diff(grid.dim['zh'][:])

zlay = grid.dim['z']
zlev = grid.dim['zh']

# read temperature, humidity
qt = np.fromfile(path+"qt.{:07d}".format(time),TF).reshape(dims)
tlay = np.fromfile(path+"T.{:07d}".format(time),TF).reshape(dims)

# convert qt from kg/kg to vmr
h2o = qt / (ep - ep*qt)

# cloud properties and effective radius
ql = read_if_exists(path+"ql", time, dims)
lwp = ql * (dz*rhoh)[:,np.newaxis,np.newaxis] # kg/m2

qi = read_if_exists(path+"qi", time, dims)
iwp = qi * (dz*rhoh)[:,np.newaxis,np.newaxis] # kg/m2

ftpnr_w = (4./3) * np.pi * 100e6 * 1e3
ftpnr_i = (4./3) * np.pi * 1e5 * 7e2
sig_fac = np.exp(np.log(1.34)*np.log(1.34))
rel = np.where(lwp>0, 1e6 * sig_fac * (lwp / dz[:,np.newaxis,np.newaxis] / ftpnr_w)**(1./3), 0)
rei = np.where(iwp>0, 1e6 * sig_fac * (iwp / dz[:,np.newaxis,np.newaxis] / ftpnr_i)**(1./3), 0)

rel = np.maximum(2.5, np.minimum(rel, 21.5))
rei = np.maximum(10., np.minimum(rei, 180.))
lwp *= 1e3 # g/m2
iwp *= 1e3 # g/m2

# ozone profile
o3 = nc_inp['init']['o3']

# read bg profile
h2o_bg = nc_inp['radiation']['h2o']
o3_bg = nc_inp['radiation']['o3']
zlay_bg = nc_inp['radiation']['z_lay']
zlev_bg = nc_inp['radiation']['z_lev']
play_bg = nc_inp['radiation']['p_lay']
plev_bg = nc_inp['radiation']['p_lev']
tlay_bg = nc_inp['radiation']['t_lay']
tlev_bg = nc_inp['radiation']['t_lev']


if f_regrid:
    zh_regrid = np.arange(0, regrid_zmax+.001, regrid_dz)
    z_regrid = (zh_regrid[1:] + zh_regrid[:-1]) / 2.
    grid_top_idx = np.where(grid.dim['z'] > regrid_zmax)[0][0]

    lwp_interp = interp1d(grid.dim['z'], lwp/dz[:,np.newaxis,np.newaxis], axis=0)
    lwp_grid = lwp_interp(z_regrid)*regrid_dz
    lwp = np.append(lwp_grid, np.tile(lwp[grid_top_idx:].mean(axis=(1,2))[:,None,None], (1, jtot, itot)), axis=0)

    iwp_interp = interp1d(grid.dim['z'], iwp/dz[:,np.newaxis,np.newaxis], axis=0)
    iwp_grid = iwp_interp(z_regrid)*regrid_dz
    iwp = np.append(iwp_grid, np.tile(iwp[grid_top_idx:].mean(axis=(1,2))[:,None,None], (1, jtot, itot)), axis=0)

    rel_interp = interp1d(grid.dim['z'], rel, axis=0)
    rel_grid = rel_interp(z_regrid)
    rel = np.append(rel_grid, np.tile(rel[grid_top_idx:].mean(axis=(1,2))[:,None,None], (1, jtot, itot)), axis=0)

    rei_interp = interp1d(grid.dim['z'], rei, axis=0)
    rei_grid = rei_interp(z_regrid)
    rei = np.append(rei_grid, np.tile(rei[grid_top_idx:].mean(axis=(1,2))[:,None,None], (1, jtot, itot)), axis=0)

    tlay_interp = interp1d(grid.dim['z'], tlay, axis=0)
    tlay_grid = tlay_interp(z_regrid)
    tlay = np.append(tlay_grid, np.tile(tlay[grid_top_idx:].mean(axis=(1,2))[:,None,None], (1, jtot, itot)), axis=0)

    h2o_interp = interp1d(grid.dim['z'], h2o, axis=0)
    h2o_grid = h2o_interp(z_regrid)
    h2o = np.append(h2o_grid, np.tile(h2o[grid_top_idx:].mean(axis=(1,2))[:,None,None], (1, jtot, itot)), axis=0)

    o3_interp = interp1d(grid.dim['z'], o3)
    o3_grid = o3_interp(z_regrid)
    o3 = np.append(o3_grid, o3[grid_top_idx:])

    p_hf = np.append(play,plev)
    z_hf = np.append(grid.dim['z'][:], grid.dim['zh'][:])
    idc = np.argsort(z_hf)
    p_hf = p_hf[idc]
    z_hf = z_hf[idc]

    play_interp = interp1d(z_hf, p_hf)
    play_grid = play_interp(z_regrid)
    play = np.append(play_grid, play[grid_top_idx:])

    plev_interp = interp1d(z_hf, p_hf)
    plev_grid = plev_interp(zh_regrid)
    plev = np.append(plev_grid, plev[grid_top_idx+1:])

    zlay = np.append(z_regrid, zlay[grid_top_idx:])
    ktot = len(play)

# find lowest height in bg profile that is heigher than domain top
z_tod = grid.dim['zh'][-1]
zmin_idx = np.where(zlay_bg[:] > z_tod)[0][0]

zlay = np.append(zlay, zlay_bg[zmin_idx:])
# patch pressure profiles
play = np.append(play, play_bg[zmin_idx:])
plev = np.append(plev, plev_bg[zmin_idx+1:])

if merge_bg:
    # for tilted column
    nz = len(play)
    grid_z = np.append(grid.dim['z'], zlay_bg[zmin_idx:])
    grid_zh= np.append(grid.dim['zh'], zlev_bg[zmin_idx+1:])
    if z_tilt == -1:
        z_tilt = grid_zh[-1]
elif f_regrid:
    nz = len(z_regrid)
    grid_z = z_regrid
    grid_zh = zh_regrid
else:
    nz = ktot
    grid_z = grid.dim['z']
    grid_zh = grid.dim['zh']

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
nc_x = nc_out.createVariable("x", "f8", ("x",))
nc_x[:] = grid.dim['x'][:]
nc_y = nc_out.createVariable("y", "f8", ("y",))
nc_y[:] = grid.dim['y'][:]
nc_z = nc_out.createVariable("z", "f8", ("z",))
nc_z[:] = grid_z

grid_xh = np.append(grid.dim['xh'][:], grid.dim['xh'][-1] +  (grid.dim['xh'][1] -  grid.dim['xh'][0]))
grid_yh = np.append(grid.dim['yh'][:], grid.dim['yh'][-1] +  (grid.dim['yh'][1] -  grid.dim['yh'][0]))

nc_xh = nc_out.createVariable("xh", "f8", ("xh",))
nc_xh[:] = grid_xh
nc_yh = nc_out.createVariable("yh", "f8", ("yh",))
nc_yh[:] = grid_yh
nc_zh = nc_out.createVariable("zh", "f8", ("zh",))
nc_zh[:] = grid_zh

nc_zlay = nc_out.createVariable("zlay", "f8", ("lay"))
print(nc_zlay[:].shape,zlay.shape)
nc_zlay[:] = zlay

# write pressures
nc_play = nc_out.createVariable("p_lay", "f8", ("lay","y","x"))
nc_play[:] = np.tile(play.reshape(len(play),1,1), (1, jtot, itot))
nc_plev = nc_out.createVariable("p_lev", "f8", ("lev","y","x"))
nc_plev[:] = np.tile(plev.reshape(len(plev),1,1), (1, jtot, itot))

# write ozone
nc_o3 = nc_out.createVariable("vmr_o3", "f8", ("lay","y","x"))
nc_o3[:] = np.tile(np.append(o3[:], o3_bg[zmin_idx:])[:,None,None], (1, jtot, itot))

# remaining 3D variables
nc_h2o = nc_out.createVariable("vmr_h2o", "f8", ("lay","y","x"))
nc_h2o[:] = np.append(h2o[:], np.tile(h2o_bg[zmin_idx:][:,None,None], (1, jtot, itot)), axis=0)
nc_tlay = nc_out.createVariable("t_lay", "f8", ("lay","y","x"))
nc_tlay[:] = np.append(tlay[:], np.tile(tlay_bg[zmin_idx:][:,None,None], (1, jtot, itot)), axis=0)

# don't bother about longwave :)
nc_tlev = nc_out.createVariable("t_lev", "f8", ("lev","y","x"))
nc_tlev[:] = 0 #np.append(tlev[:], np.tile(tlev_bg[zmin_idx+1:,None,None], (1, jtot, itot)), axis=0)

nc_lwp = nc_out.createVariable("lwp" , "f8", ("lev","y","x"))
nc_lwp[:] = 0
nc_lwp[:ktot] = lwp

nc_rel = nc_out.createVariable("rel" , "f8", ("lev","y","x"))
nc_rel[:] = 0
nc_rel[:ktot] = rel

nc_iwp = nc_out.createVariable("iwp" , "f8", ("lev","y","x"))
nc_iwp[:] = 0
nc_iwp[:ktot] = iwp

nc_rei = nc_out.createVariable("rei" , "f8", ("lev","y","x"))
nc_rei[:] = 0
nc_rei[:ktot] = rei

# surface properties
nc_alb_dir = nc_out.createVariable("sfc_alb_dir", "f8", ("y","x","band_sw"))
nc_alb_dir[:] = nl['radiation']['sfc_alb_dir']
nc_alb_dif = nc_out.createVariable("sfc_alb_dif", "f8", ("y","x","band_sw"))
nc_alb_dif[:] = nl['radiation']['sfc_alb_dif']
nc_emis = nc_out.createVariable("emis_sfc", "f8", ("y","x","band_lw"))
nc_emis[:] = nl['radiation']['emis_sfc']
nc_tsfc = nc_out.createVariable("t_sfc", "f8", ("y","x"))
nc_tsfc[:] = 0 # don't bother about longwave for now

# solar angles
nc_mu = nc_out.createVariable("mu0", "f8", ("y","x"))
nc_mu[:] = mu0
nc_az = nc_out.createVariable("azi", "f8", ("y","x"))
nc_az[:] = azi
nc_ts = nc_out.createVariable("tsi_scaling", "f8")
nc_ts[:] = tsi_scaling

#trace gases:
for var in nc_inp['radiation'].variables:
    if len(nc_inp['radiation'][var].dimensions) == 0:
        nc_gas = nc_out.createVariable("vmr_"+var, "f8")
        nc_gas[:] = nc_inp['radiation'][var][:]

nc_z_tilt = nc_out.createVariable("z_tilt", "f8")
nc_z_tilt[:] = z_tilt

nc_ng_x = nc_out.createVariable("ngrid_x", "f8")
nc_ng_x[:] = ng_x
nc_ng_y = nc_out.createVariable("ngrid_y", "f8")
nc_ng_y[:] = ng_y
nc_ng_z = nc_out.createVariable("ngrid_z", "f8")
nc_ng_z[:] = ng_z

nc_out.close()
