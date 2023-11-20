import numpy as np
import netCDF4 as nc

float_type = "f8"

# Get number of vertical levels and size from .ini file
with open('arm.ini') as f:
    for line in f:
        if line.split('=')[0] == 'ktot':
            kmax = int(line.split('=')[1])
        if line.split('=')[0]=='zsize':
            zsize = float(line.split('=')[1])

dz = zsize / kmax

# set the height
z   = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl = np.zeros(np.size(z))
qt  = np.zeros(np.size(z))
u   = np.zeros(np.size(z))
ug  = np.zeros(np.size(z))

for k in range(kmax):
    # temperature
    if z[k] <= 50.:
        thl[k] = 299.0  + (z[k]     )*(301.5 -299.0 )/(50.)
        qt[k] = 15.20  + (z[k]     )*(15.17 -15.20 )/(50.)
    elif z[k] <=  350.:
        thl[k] = 301.5  + (z[k]-  50.)*(302.5 -301.5 )/(350.-50.)
        qt[k] = 15.17  + (z[k]-  50.)*(14.98 -15.17 )/(350.-50.)
    elif z[k] <=  650.:
        thl[k] = 302.5  + (z[k]- 350.)*(303.53-302.5 )/(650.-350.)
        qt[k] = 14.98  + (z[k]- 350.)*(14.80 -14.98 )/(650.-350.)
    elif z[k] <=  700.:
        thl[k] = 303.53 + (z[k]- 650.)*(303.7 -303.53)/(700.-650.)
        qt[k] = 14.80  + (z[k]- 650.)*(14.70 -14.80 )/(700.-650.)
    elif z[k] <= 1300.:
        thl[k] = 303.7  + (z[k]- 700.)*(307.13-303.7 )/(1300.-700.)
        qt[k] = 14.70  + (z[k]- 700.)*( 13.50-14.80 )/(1300.-700.)
    elif z[k] <= 2500.:
        thl[k] = 307.13 + (z[k]-1300.)*(314.0 -307.13)/(2500.-1300.)
        qt[k] = 13.50  + (z[k]-1300.)*( 3.00 - 13.50)/(2500.-1300.)
    elif z[k] <= 5500.:
        thl[k] = 314.0  + (z[k]-2500.)*(343.2 -314.0 )/(5500.-2500.)
        qt[k] =  3.00

    # u-wind component
    u[:] = 10.

    # ug-wind component
    ug[k] = 10.

# normalize profiles to SI
qt /= 1000.  # g to kg

# set the time series
time_surface = np.array([  0.,   4.,  6.5,  7.5,  10., 12.5, 14.5])

H  = np.array([-30.,  90., 140., 140., 100., -10.,  -10])
LE = np.array([  5., 250., 450., 500., 420., 180.,    0])

advthl = np.array([ 0.   , 0.  ,  0.  , -0.08, -0.16, -0.16])
radthl = np.array([-0.125, 0.  ,  0.  ,  0.  ,  0.   , -0.1])
advqt  = np.array([ 0.08 , 0.02, -0.04, -0.10, -0.16, -0.30])

time_ls = np.array([  0.,   3.,  6.,  9.,  12., 14.5])
thlls  = np.zeros((time_ls.size, kmax))
qtls   = np.zeros((time_ls.size, kmax))

# large scale forcings
for n in range(time_ls.size):
    tendthl = advthl[n] + radthl[n]
    tendqt  = advqt[n]
    for k in range(kmax):
        if z[k] <= 1000.:
            thlls[n,k] = tendthl
            qtls[n,k] = tendqt
        else:
            thlls[n,k] = tendthl - (z[k]-1000.)*(tendthl)/(5500.-1000.)
            qtls[n,k] = tendqt - (z[k]-1000.)*(tendqt)/(5500.-1000.)

time_ls *= 3600. # h to s
thlls /= 3600. # h to s
qtls  /= 3600. # h to s
qtls  /= 1000. # g to kg

# Calculate the surface fluxes in the correct units.
Rd  = 287.
cp  = 1005.
Lv  = 2.5e6
p0  = 97000.
rho = p0/(Rd*thl[0]*(1. + 0.61*qt[0]))
time_surface *= 3600. # h to s
sbotthl = H/(rho*cp)
sbotqt  = LE/(rho*Lv)

# Save all the input data to NetCDF
nc_file = nc.Dataset("arm_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, ("z"))
nc_z[:] = z[:]

# Create a group called "init" for the initial profiles.
nc_group_init = nc_file.createGroup("init")

nc_thl = nc_group_init.createVariable("thl",   float_type, ("z"))
nc_qt  = nc_group_init.createVariable("qt",    float_type, ("z"))
nc_u   = nc_group_init.createVariable("u",     float_type, ("z"))
nc_ug  = nc_group_init.createVariable("u_geo", float_type, ("z"))
nc_thl[:] = thl[:]
nc_qt [:] = qt [:]
nc_u  [:] = u  [:]
nc_ug [:] = ug [:]

# Create a group called "timedep" for the timedep.
nc_group_timedep = nc_file.createGroup("timedep")
nc_group_timedep.createDimension("time_surface", time_surface.size)
nc_group_timedep.createDimension("time_ls", time_ls.size)

nc_time_surface = nc_group_timedep.createVariable("time_surface", float_type, ("time_surface"))
nc_thl_sbot = nc_group_timedep.createVariable("thl_sbot", float_type, ("time_surface"))
nc_qt_sbot  = nc_group_timedep.createVariable("qt_sbot" , float_type, ("time_surface"))
nc_time_surface[:] = time_surface[:]
nc_thl_sbot    [:] = sbotthl     [:]
nc_qt_sbot     [:] = sbotqt      [:]

nc_time_ls = nc_group_timedep.createVariable("time_ls", float_type, ("time_ls"))
nc_thl_ls  = nc_group_timedep.createVariable("thl_ls" , float_type, ("time_ls", "z"))
nc_qt_ls   = nc_group_timedep.createVariable("qt_ls"  , float_type, ("time_ls", "z"))
nc_time_ls[:]   = time_ls[:]
nc_thl_ls [:,:] = thlls  [:,:]
nc_qt_ls  [:,:] = qtls   [:,:]

nc_file.close()
