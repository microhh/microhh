import numpy as np
import netCDF4 as nc

float_type = "f8"

# Get number of vertical levels and size from .ini file
with open('dycoms.ini') as f:
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
ug    = np.zeros(np.size(z))
v     = np.zeros(np.size(z))
vg    = np.zeros(np.size(z))
wls   = np.zeros(np.size(z))
thlls = np.zeros(np.size(z))
qtls  = np.zeros(np.size(z))

for k in range(kmax):
    # temperature
    if(z[k] <= 840.):
        thl[k] = 289.0
    else:
        thl[k] = 297.5 + (z[k]-840.)** (1./3.)

    # specific humidity
    if(z[k] <= 840.):
        qt[k] = 1e-3*9.0
    else:
        qt[k] = 1.e-3*1.5

    wls[k] = -3.75E-6 * z[k]

    # u-wind component
    u[k] = 6

    # ug-wind component
    ug[k] = 7

    # u-wind component
    v[k] = -4.25

    # ug-wind component
    vg[k] = -5.5

# write the data to a file
# proffile = open('dycoms.prof', 'w')
# proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s} {5:^20s} {6:^20s} {7:^20s}\n'.format('z', 'thl', 'qt', 'u', 'ug', 'v', 'vg', 'wls'))
# for k in range(kmax):
#     proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E} {6:1.14E} {7:1.14E}\n'.format(z[k], thl[k], qt[k], u[k], ug[k], v[k], vg[k], wls[k]))
# proffile.close()

nc_file = nc.Dataset("dycoms_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", kmax)
nc_z  = nc_file.createVariable("z", float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_thl   = nc_group_init.createVariable("thl"  , float_type, ("z"))
nc_qt    = nc_group_init.createVariable("qt"   , float_type, ("z"))
nc_u     = nc_group_init.createVariable("u"    , float_type, ("z"))
nc_u_geo = nc_group_init.createVariable("u_geo", float_type, ("z"))
nc_v     = nc_group_init.createVariable("v"    , float_type, ("z"))
nc_v_geo = nc_group_init.createVariable("v_geo", float_type, ("z"))
nc_w_ls  = nc_group_init.createVariable("w_ls" , float_type, ("z"))

nc_z    [:] = z  [:]
nc_thl  [:] = thl[:]
nc_qt   [:] = qt [:]
nc_u    [:] = u  [:]
nc_u_geo[:] = ug [:]
nc_v    [:] = v  [:]
nc_v_geo[:] = vg [:]
nc_w_ls [:] = wls[:]

nc_file.close()
