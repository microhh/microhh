import microhh_tools as mht  # available in microhh/python directory
import numpy as np
import netCDF4 as nc
def read_if_exists(var, t, dims):
    try:
        return np.fromfile("{}.{:07d}".format(var,t),TF).reshape(dims)
    except:
        return np.zeros(dims)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-t", type=int, nargs='?')

args = parser.parse_args()

# constants
ep = 0.622
TF = np.float32
g = 9.81
#ngrid:
ng_x = 48
ng_y = 48
ng_z = 32

# input: name of case, time step
name = "cabauw"
t = args.t

### open all necessary files
# read namelist
nl = mht.Read_namelist("{}.ini".format(name))
# read stats
nc_stat = nc.Dataset("{}.default.0000000.nc".format(name))
t_idx = np.where(nc_stat['time'][:] == t)[0][0]
# read input
nc_inp = nc.Dataset("{}_input.nc".format(name))

# pressure and base state density profiles
play = nc_stat['thermo']['phydro'][t_idx]
plev = nc_stat['thermo']['phydroh'][t_idx]
rhoh = nc_stat['thermo']['rhoref'][:]

# height diff
zh = nc_stat['zh'][:]
dz = np.diff(zh)

### Read data
# dimensions and grid
itot = nl['grid']['itot']
jtot = nl['grid']['jtot']
ktot = nl['grid']['ktot']
dims = (ktot, jtot, itot)
grid = mht.Read_grid(itot, jtot, ktot)

# read temperature, humidity
qt = np.fromfile("qt.{:07d}".format(t),TF).reshape(dims)
tlay = np.fromfile("T.{:07d}".format(t),TF).reshape(dims)

# convert qt from kg/kg to vmr
h2o = qt / (ep - ep*qt)

# cloud properties and effective radius
ql = read_if_exists("ql", t, dims)
lwp = ql * (dz*rhoh)[:,np.newaxis,np.newaxis] # kg/m2

qi = read_if_exists("qi", t, dims)
iwp = qi * (dz*rhoh)[:,np.newaxis,np.newaxis] # kg/m2

ftpnr_w = (4./3) * np.pi * 100e6 * 1e3
ftpnr_i = (4./3) * np.pi * 1e5 * 7e2
sig_fac = np.exp(np.log(1.34)*np.log(1.34))
rel = np.where(lwp>0, 1e6 * sig_fac * (lwp / dz[:,np.newaxis,np.newaxis] / ftpnr_w)**(1./3), 0)
rei = np.where(iwp>0, 1e6 * sig_fac * (iwp / dz[:,np.newaxis,np.newaxis] / ftpnr_i)**(1./3), 0)
rel = np.maximum(2.5, np.minimum(rel, 21.5))
rei = np.maximum(10., np.minimum(rel, 180.))
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

# find lowest height in bg profile that is heigher than domain top
zmin_idx = np.where(zlay_bg[:] > grid.dim['zh'][-1])[0][0]

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
nc_out.createDimension("z", ktot)

# write raytracing grids
nc_x = nc_out.createVariable("x", "f8", ("x",))
nc_x[:] = grid.dim['x'][:]
nc_y = nc_out.createVariable("y", "f8", ("y",))
nc_y[:] = grid.dim['y'][:]
nc_z = nc_out.createVariable("z", "f8", ("z",))
nc_z[:] = grid.dim['z'][:]

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
#print(np.tile(h2o_bg[zmin_idx:][:,None,None], (1, jtot, itot)),h2o)
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
nc_mu[:] = np.cos(nc_stat['radiation']['sza'][t_idx])
nc_az = nc_out.createVariable("azi", "f8", ("y","x"))
nc_az[:] = nc_stat['radiation']['saa'][t_idx]*1.1
nc_ts = nc_out.createVariable("tsi_scaling", "f8")
nc_ts[:] = nc_stat['radiation']['tsi_scaling'][t_idx]

#trace gases:
for var in nc_inp['radiation'].variables:
    if len(nc_inp['radiation'][var].dimensions) == 0:
        print(var)
        nc_gas = nc_out.createVariable("vmr_"+var, "f8")
        nc_gas[:] = nc_inp['radiation'][var][:]

nc_ng_x = nc_out.createVariable("ngrid_x", "f8")
nc_ng_x[:] = ng_x
nc_ng_y = nc_out.createVariable("ngrid_y", "f8")
nc_ng_y[:] = ng_y
nc_ng_z = nc_out.createVariable("ngrid_z", "f8")
nc_ng_z[:] = ng_z

nc_out.close()
