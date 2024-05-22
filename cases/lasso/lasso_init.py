import netCDF4 as nc
import numpy as np
import f90nml
import shutil
import os
from datetime import datetime

#import os.path
#import os
#import sys
#import glob
float_type = "f8"
largescale = True
# Define some constants
cp = 1004.
rd = 287.04
grav = 9.8
rho = 1.225
p0 = 1e5
Lv = 2.5e6
tau = 21600.

# Get number of vertical levels and size from .ini file
shutil.copy2('testbed.ini', 'testbed.tmp')
with open('testbed.tmp') as f:
    for line in f:
        if(line.split('=')[0] == 'ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0] == 'zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax
zstretch = 5800.
stretch = 1.04
# Read WRF Namelist
fnml = "config/namelist.input"
nml = f90nml.read(fnml)
runtime_in = nml["time_control"]["run_days"] * 86400 + nml["time_control"]["run_hours"] * \
    3600 + nml["time_control"]["run_minutes"] * 60 + nml["time_control"]["run_seconds"]
# Read Surface pressure
fname_in = "config/wrfinput_d01.nc"

f = nc.Dataset(fname_in, 'r+')
ps_in = f.variables['PB'][:][0, 0, 0, 0]
#wrfstattime = f.variables['XTIME'][:]*60.
f.close()
# Read WRF Surface Forcing
fname_in = "config/input_sfc_forcing.nc"
f = nc.Dataset(fname_in, 'r+')
H_in = f.variables['PRE_SH_FLX'][:]
LE_in = f.variables['PRE_LH_FLX'][:]
f.close()
# Read WRF LS Forcing & Nudging
fname_in = "config/input_ls_forcing.nc"
f = nc.Dataset(fname_in, 'r+')
timestr = f.variables['Times'][:]
z_in = f.variables['Z_LS'][:]
u_in = f.variables['U_LS'][:]
v_in = f.variables['V_LS'][:]
thl_in = f.variables['TH_RLX'][:]
qt_in = f.variables['QV_RLX'][:]
thlls_in = f.variables['TH_ADV'][:]
qtls_in = f.variables['QV_ADV'][:]
wls_in = f.variables['W_LS'][:]
f.close()

str = np.chararray(timestr.shape)
dt = np.empty(timestr.shape[0], dtype=datetime)
tnudge = np.zeros(timestr.shape[0])
for i in range(timestr.shape[0]):
    dt[i] = datetime.strptime(
        timestr[i].tostring().decode(),
        '%Y-%m-%d_%H:%M:%S')
    tnudge[i] = (dt[i] - dt[0]).total_seconds()
ntnudge = tnudge.size
#
# interpolate onto Microhh grid
ntnudge = timestr.shape[0]
# non-equidistant grid

z = np.zeros(kmax)
z[0] = 0.5 * dz
for k in range(1, kmax):
    z[k] = z[k - 1] + 0.5 * dz
    if z[k] > zstretch:
        dz *= stretch
    z[k] += 0.5 * dz

zh = np.zeros(kmax)
zh[:-1] = (z[1:] + z[:-1]) / 2
zh[-1] = 2 * z[-1] - zh[-2]
u = np.zeros((ntnudge, z.size))
v = np.zeros(np.shape(u))
thl = np.zeros(np.shape(u))
qt = np.zeros(np.shape(u))
thlls = np.zeros(np.shape(u))
qtls = np.zeros(np.shape(u))
wls = np.zeros(np.shape(u))
nudge_factor = np.zeros(np.shape(u))

p_sbot = np.zeros((ntnudge))
#p_sbot = np.interp(tnudge, wrfstattime,ps_in)
p_sbot[:] = ps_in
for t in range(tnudge.size):
    u[t, :] = np.interp(z, z_in[t, :], u_in[t, :])
    v[t, :] = np.interp(z, z_in[t, :], v_in[t, :])
    thl[t, :] = np.interp(z, z_in[t, :], thl_in[t, :])
    qt[t, :] = np.interp(z, z_in[t, :], qt_in[t, :])
    thlls[t, :] = np.interp(z, z_in[t, :], thlls_in[t, :])
    qtls[t, :] = np.interp(z, z_in[t, :], qtls_in[t, :])
    wls[t, :] = np.interp(zh, z_in[t, :], wls_in[t, :])

ug = u
vg = v
nudge_factor[:, :] = 1. / tau
# Surface fluxes
rhosurf = p_sbot / (rd * thl[:, 0] * (1. + 0.61 * qt[:, 0]))

sbotthl = H_in / (rhosurf * cp)
sbotqt = LE_in / (rhosurf * Lv)


# Modify .ini file: Add comments for case description; Alter lines where
# necessary.

inifile = open('testbed.ini', 'w')
inifile.write("#Converted from LASSO WRF" + "\n")
# inifile.write("#Start Date = " + timestr[0]+ "\n")
# inifile.write("#Stop Date = "  + timestr[-1]+"\n")
with open('testbed.tmp') as f:
    for line_in in f:
        if(line_in.split('=')[0] == 'zsize'):
            line = "zsize={0:f}\n".format(zh[-1])
        elif(line_in.split('=')[0] == 'pbot'):
            line = "pbot={0:f}\n".format(p_sbot[0])
        else:
            line = line_in
        inifile.write(line)
inifile.close()
os.remove('testbed.tmp')


# Save all the input data to NetCDF
nc_file = nc.Dataset(
    "testbed_input.nc",
    mode="w",
    datamodel="NETCDF4",
    clobber=True)

nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, ("z"))
nc_z[:] = z[:]

nc_file.createDimension("zh", kmax)
nc_zh = nc_file.createVariable("zh", float_type, ("zh"))
nc_zh[:] = zh[:]

# Create a group called "init" for the initial profiles.
nc_group_init = nc_file.createGroup("init")

nc_thl = nc_group_init.createVariable("thl", float_type, ("z"))
nc_qt = nc_group_init.createVariable("qt", float_type, ("z"))
nc_u = nc_group_init.createVariable("u", float_type, ("z"))
nc_v = nc_group_init.createVariable("v", float_type, ("z"))
nc_ug = nc_group_init.createVariable("u_geo", float_type, ("z"))
nc_vg = nc_group_init.createVariable("v_geo", float_type, ("z"))

nc_nudge_factor = nc_group_init.createVariable("nudgefac", float_type, ("z"))
nc_thl[:] = thl[0, :]
nc_qt[:] = qt[0, :]
nc_u[:] = u[0, :]
nc_v[:] = v[0, :]
nc_ug[:] = ug[0, :]
nc_vg[:] = vg[0, :]

nc_nudge_factor[:] = nudge_factor[0, :]

# Create a group called "timedep" for the timedep.
nc_group_timedep = nc_file.createGroup("timedep")

nc_group_timedep.createDimension("time", tnudge.size)
nc_time = nc_group_timedep.createVariable("time", float_type, ("time"))
nc_time[:] = tnudge[:]


nc_thl_sbot = nc_group_timedep.createVariable("thl_sbot", float_type, ("time"))
nc_qt_sbot = nc_group_timedep.createVariable("qt_sbot", float_type, ("time"))
nc_p_sbot = nc_group_timedep.createVariable("p_sbot", float_type, ("time"))
nc_thl_sbot[:] = sbotthl[:]
nc_qt_sbot[:] = sbotqt[:]
nc_p_sbot[:] = sbotqt[:]

nc_thl_ls = nc_group_timedep.createVariable(
    "thl_ls", float_type, ("time", "z"))
nc_qt_ls = nc_group_timedep.createVariable("qt_ls", float_type, ("time", "z"))
nc_w_ls = nc_group_timedep.createVariable("w_ls", float_type, ("time", "zh"))
nc_u_g = nc_group_timedep.createVariable("u_geo", float_type, ("time", "z"))
nc_v_g = nc_group_timedep.createVariable("v_geo", float_type, ("time", "z"))
nc_u_nudge = nc_group_timedep.createVariable(
    "u_nudge", float_type, ("time", "z"))
nc_v_nudge = nc_group_timedep.createVariable(
    "v_nudge", float_type, ("time", "z"))
nc_thl_nudge = nc_group_timedep.createVariable(
    "thl_nudge", float_type, ("time", "z"))
nc_qt_nudge = nc_group_timedep.createVariable(
    "qt_nudge", float_type, ("time", "z"))
nc_thl_ls[:, :] = thlls[:, :]
nc_qt_ls[:, :] = qtls[:, :]
nc_w_ls[:, :] = wls[:, :]
nc_u_g[:, :] = ug[:, :]
nc_v_g[:, :] = vg[:, :]
nc_u_nudge[:, :] = u[:, :]
nc_v_nudge[:, :] = v[:, :]
nc_thl_nudge[:, :] = thl[:, :]
nc_qt_nudge[:, :] = qt[:, :]

nc_file.close()

print("done")
