import netCDF4 as nc
import numpy as np

nc_file = nc.Dataset("rcemip_default_vert_300.nc", "r")
itot = 96
jtot = 96
ktot = 144

"""
# Add the 0D variables.
nc_0d = nc.Dataset("microhh_vert_300_0d.nc", "w")
nc_0d.createDimension("time", 100*24+1)

time_var = nc_0d.createVariable("time", np.float32, ("time"))
time_var[:] = nc_file.variables["time"][::4] / 86400.

def add_0d_variable(data_in, name, long_name, units):
    nc_var = nc_0d.createVariable(name, np.float32, ("time"))
    nc_var[0] = data_in[0]
    nc_var[1:] = (data_in[1:].reshape(100*24, 4)).mean(axis=-1)
    nc_var.units = units
    nc_var.long_name = long_name

pr_avg = nc_file.groups["thermo"].variables["rr"][:]
add_0d_variable(pr_avg, "pr_avg", "domain avg. surface precipitation rate", "kg m-2 s-1")
del(pr_avg)

Lv = 2.45e6
hfls_avg = nc_file.groups["default"].variables["qt_flux"][:,0] * nc_file.groups["thermo"].variables["rhoh"][:,0] * Lv
add_0d_variable(hfls_avg, "hfls_avg", "domain avg. surface upward latent heat flux", "W m-2")

cp = 1004.
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
"""

#########################
# Add the 1D variables. #
#########################
nc_1d = nc.Dataset("microhh_vert_300_1d.nc", "w")
nc_1d.createDimension("time", 100*24+1)
nc_1d.createDimension("z" , ktot  )
nc_1d.createDimension("zh", ktot+1)

time_var = nc_1d.createVariable("time", np.float32, ("time"))
time_var[:] = nc_file.variables["time"][::4] / 86400.

z_var = nc_1d.createVariable("z", np.float32, ("z"))
z_var[:] = nc_file.variables["z"][:]

zh_var = nc_1d.createVariable("zh", np.float32, ("zh"))
zh_var[:] = nc_file.variables["zh"][:]

def add_1d_variable(data_in, name, long_name, units, z_name):
    nc_var = nc_1d.createVariable(name, np.float32, ("time", z_name))
    nc_var[0,:] = data_in[0,:]
    ktot_loc = ktot if z_name == "z" else ktot+1
    nc_var[1:] = (data_in[1:].reshape(100*24, 4, ktot_loc)).mean(axis=1)
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

