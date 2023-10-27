import numpy as np
import netCDF4 as nc

float_type = 'f8'

case = 'gcss' # Original RICO
#case = 'ss08' # Moist RICO from Stevens/Seifert & Seifert/Heus
#case = 'test' # More moist mixed-layer for testing

# Get number of vertical levels and size from .ini file
with open('rico.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

# set the height
z     = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl   = np.zeros(z.size)
qt    = np.zeros(z.size)
u     = np.zeros(z.size)
ugeo  = np.zeros(z.size)
v     = np.zeros(z.size)
vgeo  = np.zeros(z.size)
wls   = np.zeros(z.size)
thlls = np.zeros(z.size)
qtls  = np.zeros(z.size)

#print('Setup = %s'%case)
for k in range(kmax):

    # Liquid water potential temperature: same in GCSS and SS08
    if(z[k] < 740.):
        thl[k] = 297.9
    else:
        thl[k] = 297.9 + (317.0 - 297.9)/(4000. - 740.) * (z[k] - 740.) 

    if(case == 'gcss'):
        if(z[k] < 740.):
            qt[k] = 16.0 + (13.8 - 16.0) / 740. * z[k]
        elif(z[k] < 3260.):
            qt[k] = 13.8 + (2.4 - 13.8) / (3260. - 740.) * (z[k] - 740.) 
        else:
            qt[k] = 2.4 + (1.8 - 2.4)/(4000. - 3260.) * (z[k] - 3260.) 

    elif(case == 'ss08'):
        if(z[k] < 740.):
            qt[k] = 16.0 + (13.8 - 16.0) / 740. * z[k]
        elif(z[k] < 3260.):
            qt[k] = 13.8 + (4.4 - 13.8) / (3260. - 740.) * (z[k] - 740.) 
        else:
            qt[k] = 4.4 + (3.6 - 4.4)/(4000. - 3260.) * (z[k] - 3260.) 

    elif(case == 'test'):
        q0 = 18.
        q1 = 15.8
        if(z[k] < 740.):
            qt[k] = q0 + (q1 - q0) / 740. * z[k]
        elif(z[k] < 3260.):
            qt[k] = q1 + (2.4 - q1) / (3260. - 740.) * (z[k] - 740.) 
        else:
            qt[k] = 2.4 + (1.8 - 2.4)/(4000. - 3260.) * (z[k] - 3260.) 

    # Subsidence
    if(z[k] < 2260):
        wls[k] = -0.005 * (z[k] / 2260.)
    else:
        wls[k] = -0.005

    # U and V component wind
    u[k] = -9.9 + 2.0e-3 * z[k]
    v[k] = -3.8
    ugeo[k] = u[k]
    vgeo[k] = v[k]

    # Advective and radiative tendency thl
    thlls[k] = -2.5 / 86400.

    # Advective tendency qt
    if(z[k] < 2980):
        qtls[k] = -1.0 / 86400. + (1.3456/ 86400.) * z[k] / 2980.
    else:
        qtls[k] = 4e-6

# normalize profiles to SI
qt  /= 1000.
qtls/= 1000.

# write the data to a file
nc_file = nc.Dataset("rico_input.nc", mode="w", datamodel="NETCDF4", clobber=True)
nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_thl   = nc_group_init.createVariable("thl"   , float_type, ("z"))
nc_qt    = nc_group_init.createVariable("qt"    , float_type, ("z"))
nc_u     = nc_group_init.createVariable("u"     , float_type, ("z"))
nc_ugeo  = nc_group_init.createVariable("u_geo" , float_type, ("z"))
nc_v     = nc_group_init.createVariable("v"     , float_type, ("z"))
nc_vgeo  = nc_group_init.createVariable("v_geo" , float_type, ("z"))
nc_wls   = nc_group_init.createVariable("w_ls"  , float_type, ("z"))
nc_thlls = nc_group_init.createVariable("thl_ls", float_type, ("z"))
nc_qtls  = nc_group_init.createVariable("qt_ls" , float_type, ("z"))

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

ep = 287.04 / 461.5 

# Surface settings
def esat(T):
    c0 = 0.6105851e+03; c1 = 0.4440316e+02; c2 = 0.1430341e+01; c3 = 0.2641412e-01 
    c4 = 0.2995057e-03; c5 = 0.2031998e-05; c6 = 0.6936113e-08; c7 = 0.2564861e-11 
    c8 = -.3704404e-13 
    x  = max(-80.,T-273.15)
    return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

def qsat(p, T):
    return ep*esat(T)/(p-(1.-ep)*esat(T))

ps  = 101540.
SST = 299.8 
ths = SST / (ps/1.e5)**(287.04/1005.)
qs  = qsat(ps, SST) 
#print('sbot[thl]=%f, sbot[qt]=%f'%(ths, qs))
