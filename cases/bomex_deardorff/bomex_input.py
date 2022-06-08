import numpy as np
import netCDF4 as nc

float_type = "f8"

# Get number of vertical levels and size from .ini file
with open('bomex.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

# set the height
z     = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl   = np.zeros(np.size(z))
qt    = np.zeros(np.size(z))
u     = np.zeros(np.size(z))
v     = np.zeros(np.size(z))
ugeo  = np.zeros(np.size(z))
vgeo  = np.zeros(np.size(z))
wls   = np.zeros(np.size(z))
thlls = np.zeros(np.size(z))
qtls  = np.zeros(np.size(z))
sgstke12 = np.zeros(np.size(z))

for k in range(kmax):
    # temperature
    if(z[k] <= 520.):
        thl[k] = 298.7
    elif(z[k] <= 1480):
        thl[k] = 298.7 + (z[k]-520.)*(302.4-298.7)/(1480.-520.)
    elif(z[k] <= 2000):
        thl[k] = 302.4 + (z[k]-1480.)*(308.2-302.4)/(2000.-1480.)
    else:
        thl[k] = 308.2 + (z[k]-2000.)*(311.85-308.2)/(3000.-2000.)

    # specific humidity
    if(z[k] <= 520.):
        qt[k] = 1e-3*(17.0 + z[k]*(16.3-17.0)/520.)
    elif(z[k] <= 1480):
        qt[k] = 1.e-3*(16.3 + (z[k]-520.)*(10.7-16.3)/(1480.-520.))
    elif(z[k] <= 2000):
        qt[k] = 1.e-3*(10.7 + (z[k]-1480.)*(4.2-10.7)/(2000.-1480.))
    else:
        qt[k] = 1.e-3*(4.2 + (z[k]-2000.)*(3.-4.2)/(3000.-2000.))


    # u-wind component
    if(z[k] <= 700.):
        u[k] = -8.75
    else:
        u[k] = -8.75 + (z[k]-700.)*(-4.61+8.75)/(3000.-700.)

    # ugeo-wind component
    ugeo[k] = -10. + 1.8e-3*z[k]

    # large scale vertical velocity
    if(z[k] <= 1500):
        wls[k] = z[k]*(-0.65)/1500.
    elif(z[k] <= 2100):
        wls[k] = -0.65 + (z[k]-1500)*(0.65)/(2100.-1500.)

    # large scale temperature tendency
    if(z[k] <= 1500):
        thlls[k] = (-2.)
    else:
        thlls[k] = (-2.) + (z[k]-1500)*(2.)/(3000.-1500.)

    # large scale moisture tendency
    if(z[k] <= 300):
        qtls[k] = -1.2
    elif(z[k] <= 500):
        qtls[k] = -1.2 + (z[k]-300)*(1.2)/(500.-300)

    # square root of SGS-TKE
    if(z[k] <= 3000.):
        sgstke12[k] = np.sqrt( (1 - z[k]/3000) )
    if(z[k] > 250.):
        sgstke12[k] = 0.

# normalize profiles to SI
#qtls /= 1000.  # from g/kg to kg/kg
wls  /= 100.   # from cm/s to m/s
thlls  /= 86400. # from K/d to K/s
qtls *= 1.e-8

nc_file = nc.Dataset("bomex_input.nc", mode="w", datamodel="NETCDF4", clobber=False)
nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_thl   = nc_group_init.createVariable("thl"   , float_type, ("z"))
nc_qt    = nc_group_init.createVariable("qt"    , float_type, ("z"))
nc_u     = nc_group_init.createVariable("u"     , float_type, ("z"))
nc_v     = nc_group_init.createVariable("v"     , float_type, ("z"))
nc_ugeo  = nc_group_init.createVariable("u_geo" , float_type, ("z"))
nc_vgeo  = nc_group_init.createVariable("v_geo" , float_type, ("z"))
nc_wls   = nc_group_init.createVariable("w_ls"  , float_type, ("z"))
nc_thlls = nc_group_init.createVariable("thl_ls", float_type, ("z"))
nc_qtls  = nc_group_init.createVariable("qt_ls" , float_type, ("z"))
nc_sgstke12 = nc_group_init.createVariable('sgstke12', float_type, ('z'))

nc_sgstke12 [:] = sgstke12[:] # square root of sgs tke

nc_z    [:] = z    [:]
nc_thl  [:] = thl  [:]
nc_qt   [:] = qt   [:]
nc_u    [:] = u    [:]
nc_ugeo [:] = ugeo [:]
nc_v    [:] = v    [:]
nc_vgeo [:] = vgeo [:]
nc_wls  [:] = wls  [:]
nc_thlls[:] = thlls[:]
nc_qtls [:] = qtls [:]

nc_file.close()
