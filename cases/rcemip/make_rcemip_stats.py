import netCDF4 as nc
import numpy as np
import shutil
import subprocess

nc_file = nc.Dataset('rcemip_default_0000000.nc', 'r')
itot = 96
jtot = 96
ktot = 144
timetot = nc_file.variables['time'].shape[0]

case_name = "vert_300"

Lv = 2.45e6
cp = 1004.
Rd = 287.
p0 = 1e5


#########################
# Add the 0D variables. #
#########################
nc_0d = nc.Dataset("microhh_{}_0d.nc".format(case_name), "w")
nc_0d.createDimension("time", timetot)

time_var = nc_0d.createVariable("time", np.float32, ("time"))
time_var[:] = nc_file.variables["time"][:] / 86400.

def add_0d_variable(data_in, name, long_name, units):
    nc_var = nc_0d.createVariable(name, np.float32, ("time"))
    #nc_var[0] = data_in[0]
    #nc_var[1:] = (data_in[1:].reshape(100*24, 4)).mean(axis=-1)
    nc_var[:] = data_in[:]
    nc_var.units = units
    nc_var.long_name = long_name

pr_avg = nc_file.groups["thermo"].variables["rr"][:]
add_0d_variable(pr_avg, "pr_avg", "domain avg. surface precipitation rate", "kg m-2 s-1")
del(pr_avg)

hfls_avg = nc_file.groups["default"].variables["qt_flux"][:,0] * nc_file.groups["thermo"].variables["rhoh"][:,0] * Lv
add_0d_variable(hfls_avg, "hfls_avg", "domain avg. surface upward latent heat flux", "W m-2")

hfss_avg = nc_file.groups["default"].variables["thl_flux"][:,0] * nc_file.groups["thermo"].variables["rhoh"][:,0] * cp
add_0d_variable(hfss_avg, "hfss_avg", "domain avg. surface upward sensible heat flux", "W m-2")

clwvi_avg = nc_file.groups["thermo"].variables["ql_path"][:]
clivi_avg = nc_file.groups["thermo"].variables["qi_path"][:]
prw_avg   = nc_file.groups["default"].variables["qt_path"][:] - clwvi_avg - clivi_avg
add_0d_variable(clwvi_avg, "clwvi_avg", "domain avg. water condensed water path", "kg m-2")
add_0d_variable(clivi_avg, "clivi_avg", "domain avg. water ice water path", "kg m-2")
add_0d_variable(prw_avg  , "prw_avg"  , "domain avg. water vapor path", "kg m-2")
del(clwvi_avg, clivi_avg, prw_avg)

sprw_avg  = nc_file.groups["thermo"].variables["qsat_path"][:]
add_0d_variable(sprw_avg, "sprw_avg", "domain avg. saturated water vapor path", "kg m-2")
del(sprw_avg)

rlds_avg = nc_file.groups["radiation"].variables["lw_flux_dn"][:,0]
rlus_avg = nc_file.groups["radiation"].variables["lw_flux_up"][:,0]
rsds_avg = nc_file.groups["radiation"].variables["sw_flux_dn"][:,0]
rsus_avg = nc_file.groups["radiation"].variables["sw_flux_up"][:,0]
add_0d_variable(rlds_avg, "rlds_avg", "domain avg. surface downwelling longwave flux", "W m-2")
add_0d_variable(rlus_avg, "rlus_avg", "domain avg. surface upwelling longwave flux", "W m-2")
add_0d_variable(rsds_avg, "rsds_avg", "domain avg. surface downwelling shortwave flux", "W m-2")
add_0d_variable(rsus_avg, "rsus_avg", "domain avg. surface upwelling shortwave flux", "W m-2")
del(rlds_avg, rlus_avg, rsds_avg, rsus_avg)

rldscs_avg = nc_file.groups["radiation"].variables["lw_flux_dn_clear"][:,0]
rluscs_avg = nc_file.groups["radiation"].variables["lw_flux_up_clear"][:,0]
rsdscs_avg = nc_file.groups["radiation"].variables["sw_flux_dn_clear"][:,0]
rsuscs_avg = nc_file.groups["radiation"].variables["sw_flux_up_clear"][:,0]
add_0d_variable(rldscs_avg, "rldscs_avg", "domain avg. surface downwelling longwave flux - clear sky", "W m-2")
add_0d_variable(rluscs_avg, "rluscs_avg", "domain avg. surface upwelling longwave flux - clear sky", "W m-2")
add_0d_variable(rsdscs_avg, "rsdscs_avg", "domain avg. surface downwelling shortwave flux - clear sky", "W m-2")
add_0d_variable(rsuscs_avg, "rsuscs_avg", "domain avg. surface upwelling shortwave flux - clear sky", "W m-2")
del(rldscs_avg, rluscs_avg, rsdscs_avg, rsuscs_avg)

rlut_avg = nc_file.groups["radiation"].variables["lw_flux_up"][:,0]
rsdt_avg = nc_file.groups["radiation"].variables["sw_flux_dn"][:,0]
rsut_avg = nc_file.groups["radiation"].variables["sw_flux_up"][:,0]
add_0d_variable(rlut_avg, "rlut_avg", "domain avg. TOA upwelling longwave flux", "W m-2")
add_0d_variable(rsdt_avg, "rsdt_avg", "domain avg. TOA downwelling shortwave flux", "W m-2")
add_0d_variable(rsut_avg, "rsut_avg", "domain avg. TOA upwelling shortwave flux", "W m-2")
del(rlut_avg, rsdt_avg, rsut_avg)

rlutcs_avg = nc_file.groups["radiation"].variables["lw_flux_up_clear"][:,-2]
rsutcs_avg = nc_file.groups["radiation"].variables["sw_flux_up_clear"][:,-2]
add_0d_variable(rlutcs_avg, "rlutcs_avg", "domain avg. TOA upwelling longwave flux - clear sky", "W m-2")
add_0d_variable(rsutcs_avg, "rsutcs_avg", "domain avg. TOA upwelling shortwave flux - clear sky", "W m-2")
del(rlutcs_avg, rsutcs_avg)

nc_0d.close()


#########################
# Add the 1D variables. #
#########################
nc_1d = nc.Dataset("microhh_{}_1d.nc".format(case_name), "w")
nc_1d.createDimension("time", timetot)
nc_1d.createDimension("z" , ktot  )
nc_1d.createDimension("zh", ktot+1)

time_var = nc_1d.createVariable("time", np.float32, ("time"))
time_var[:] = nc_file.variables["time"][:] / 86400. # Convert to days.

z_var = nc_1d.createVariable("z", np.float32, ("z"))
z_var[:] = nc_file.variables["z"][:]

zh_var = nc_1d.createVariable("zh", np.float32, ("zh"))
zh_var[:] = nc_file.variables["zh"][:]

def add_1d_variable(data_in, name, long_name, units, z_name):
    nc_var = nc_1d.createVariable(name, np.float32, ("time", z_name))
    ktot_loc = ktot if z_name == "z" else ktot+1
    #nc_var[0,:] = data_in[0,:]
    #nc_var[1:] = (data_in[1:].reshape(100*24, 4, ktot_loc)).mean(axis=1)
    nc_var[:,:] = data_in[:,:]
    nc_var.units = units
    nc_var.long_name = long_name

ta_avg = nc_file.groups["thermo"].variables["T"][:,:]
add_1d_variable(ta_avg, "ta_avg", "domain avg. air temperature profile", "K", "z")
del(ta_avg)

ua_avg = nc_file.groups["default"].variables["u"][:,:]
va_avg = nc_file.groups["default"].variables["v"][:,:]
add_1d_variable(ua_avg, "ua_avg", "domain avg. eastward wind profile", "m s-1", "z")
add_1d_variable(va_avg, "va_avg", "domain avg. northward wind profile", "m s-1", "z")
del(ua_avg, va_avg)

hus_avg = nc_file.groups["default"].variables["qt"][:,:]
hur_avg = nc_file.groups["thermo"].variables["rh"][:,:] * 100.
add_1d_variable(hus_avg, "hus_avg", "domain avg. specific humidity profile", "kg kg-1", "z")
add_1d_variable(hur_avg, "hur_avg", "domain avg. relative humidity profile", "%", "z")
del(hus_avg, hur_avg)

clw_avg = nc_file.groups["thermo"].variables["ql"][:,:]
cli_avg = nc_file.groups["thermo"].variables["qi"][:,:]
plw_avg = nc_file.groups["default"].variables["qr"][:,:]
pli_avg = nc_file.groups["default"].variables["qs"][:,:] + nc_file.groups["default"].variables["qg"][:,:]
add_1d_variable(clw_avg, "clw_avg", "domain avg. mass fraction of the cloud liquid water profile", "kg kg-1", "z")
add_1d_variable(cli_avg, "cli_avg", "domain avg. mass fraction of the cloud ice profile", "kg kg-1", "z")
add_1d_variable(plw_avg, "plw_avg", "domain avg. mass fraction of the precipitating liquid water profile", "kg kg-1", "z")
add_1d_variable(pli_avg, "pli_avg", "domain avg. mass fraction of the precipitating ice profile", "kg kg-1", "z")
del(clw_avg, cli_avg, plw_avg, pli_avg)

exner = (nc_file.groups["thermo"].variables["phydro"][:,:] / p0)**(Rd/cp)
theta_avg = nc_file.groups["thermo"].variables["T"][:,:] / exner
qv_avg = nc_file.groups["default"].variables["qt"][:,:] - nc_file.groups["thermo"].variables["ql"][:,:] - nc_file.groups["thermo"].variables["qi"][:,:]
thetae_avg = (nc_file.groups["thermo"].variables["T"][:,:] + Lv/cp*qv_avg) / exner
add_1d_variable(theta_avg, "theta_avg", "domain avg. potential temperature profile", "K", "z")
add_1d_variable(thetae_avg, "thetae_avg", "domain avg. equivalent potential temperature profile", "K", "z")
del(theta_avg, thetae_avg, exner, qv_avg)

cldfrac_avg = np.maximum(nc_file.groups["thermo"].variables["ql_frac"][:,:], nc_file.groups["thermo"].variables["qi_frac"][:,:])
add_1d_variable(cldfrac_avg, "cldfrac_avg", "global cloud fraction profile", "-", "z")
del(cldfrac_avg)


sw_flux_up = nc_file.groups["radiation"].variables["sw_flux_up"][:,:]
sw_flux_dn = nc_file.groups["radiation"].variables["sw_flux_dn"][:,:]
lw_flux_up = nc_file.groups["radiation"].variables["lw_flux_up"][:,:]
lw_flux_dn = nc_file.groups["radiation"].variables["lw_flux_dn"][:,:]
sw_clear_flux_up = nc_file.groups["radiation"].variables["sw_flux_up_clear"][:,:]
sw_clear_flux_dn = nc_file.groups["radiation"].variables["sw_flux_dn_clear"][:,:]
lw_clear_flux_up = nc_file.groups["radiation"].variables["lw_flux_up_clear"][:,:]
lw_clear_flux_dn = nc_file.groups["radiation"].variables["lw_flux_dn_clear"][:,:]

dz = nc_file.variables["zh"][1:] - nc_file.variables["zh"][:-1]
fac = - 1./(nc_file.groups["thermo"].variables["rho"][:,:] * cp * dz[None, :])
sw_heating_rate = fac * (sw_flux_up[:,1:] - sw_flux_up[:,:-1] - sw_flux_dn[:,1:] + sw_flux_dn[:,:-1])
lw_heating_rate = fac * (lw_flux_up[:,1:] - lw_flux_up[:,:-1] - lw_flux_dn[:,1:] + lw_flux_dn[:,:-1])
sw_clear_heating_rate = fac * (sw_clear_flux_up[:,1:] - sw_clear_flux_up[:,:-1] - sw_clear_flux_dn[:,1:] + sw_clear_flux_dn[:,:-1])
lw_clear_heating_rate = fac * (lw_clear_flux_up[:,1:] - lw_clear_flux_up[:,:-1] - lw_clear_flux_dn[:,1:] + lw_clear_flux_dn[:,:-1])

add_1d_variable(sw_heating_rate, "tntrs_avg", "domain avg. shortwave radiative heating rate profile", "K s-1", "z")
add_1d_variable(lw_heating_rate, "tntrl_avg", "domain avg. longwave radiative heating rate profile", "K s-1", "z")
add_1d_variable(sw_clear_heating_rate, "tntrscs_avg", "domain avg. shortwave radiative heating rate profile", "K s-1", "z")
add_1d_variable(lw_clear_heating_rate, "tntrlcs_avg", "domain avg. longwave radiative heating rate profile", "K s-1", "z")
del(dz, fac, sw_heating_rate, lw_heating_rate, sw_clear_heating_rate, lw_clear_heating_rate)


#########################
# Add the 2D variables. #
#########################

rho = nc_file.groups["thermo"].variables["rhoh"][:,0]

var_in  = "rr_bot"
var_out = "pr"
nc_file_in = "{}.xy.nc".format(var_in)
nc_file_out = "microhh_{0}_2d_{1}.nc".format(case_name, var_out)
shutil.copy(nc_file_in, nc_file_out)
nc_2d = nc.Dataset(nc_file_out, "r+")
nc_2d.renameVariable(var_in, var_out)
nc_2d_var = nc_2d.variables[var_out]
nc_2d_var.units = "kg m-2 s-1"
nc_2d_var.long_name = "surface precipitation rate"
nc_2d.close()

var_in  = "qtfluxbot"
var_out = "hfls"
nc_file_in = "{}.xy.nc".format(var_in)
nc_file_out = "microhh_{0}_2d_{1}.nc".format(case_name, var_out)
shutil.copy(nc_file_in, nc_file_out)
nc_2d = nc.Dataset(nc_file_out, "r+")
nc_2d.renameVariable(var_in, var_out)
nc_2d_var = nc_2d.variables[var_out]
nc_2d_var.units = "W m-2"
nc_2d_var.long_name = "surface upward latent heat flux"
nc_2d_var[:,:,:] *= Lv * rho[:,None,None]
nc_2d.close()

var_in  = "thlfluxbot"
var_out = "hfss"
nc_file_in = "{}.xy.nc".format(var_in)
nc_file_out = "microhh_{0}_2d_{1}.nc".format(case_name, var_out)
shutil.copy(nc_file_in, nc_file_out)
nc_2d = nc.Dataset(nc_file_out, "r+")
nc_2d.renameVariable(var_in, var_out)
nc_2d_var = nc_2d.variables[var_out]
nc_2d_var.units = "W m-2"
nc_2d_var.long_name = "surface upward sensible heat flux"
nc_2d_var[:,:,:] *= cp * rho[:,None,None]
nc_2d.close()

#ncea -d level,6,6 -F hgt.mon.mean.nc hgt500.mon.mean.nc
var_in  = "lw_flux_dn"
var_out = "rlds"
nc_file_in = "{}.xy.nc".format(var_in)
nc_file_out = "microhh_{0}_2d_{1}.nc".format(case_name, var_out)
subprocess.run("ncks -O -h -d z,0,0 {0} {1}".format(nc_file_in, nc_file_out), shell=True)
nc_2d = nc.Dataset(nc_file_out, "r+")
nc_2d.renameVariable(var_in, var_out)
nc_2d_var = nc_2d.variables[var_out]
nc_2d_var.units = "W m-2"
nc_2d_var.long_name = "surface downwelling longwave flux"
nc_2d.close()


