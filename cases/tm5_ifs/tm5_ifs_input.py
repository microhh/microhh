import numpy as np
import pylab as pl
import netCDF4 as nc

float_type = "f8"
# float_type = "f4"

case = 'orlando' # According to the protocol.

# Get number of vertical levels and size from .ini file
with open('tm5_ifs.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax
u0 = 5.             # Zonal wind speed
v0 = 0.             # Meridional wind speed
co20 = 400e3        # ppb CO2


# set the height
z  = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl   = np.zeros(z.size)
qt    = np.zeros(z.size)
u     = np.zeros(z.size)
ug    = np.zeros(z.size)
v     = np.zeros(z.size)
vg    = np.zeros(z.size)
wls   = np.zeros(z.size)
thlls = np.zeros(z.size)
qtls  = np.zeros(z.size)
hno3 = np.zeros(np.size(z))
n2o5 = np.zeros(np.size(z))
h2o2 = np.zeros(np.size(z))
co = np.zeros(np.size(z))
pan = np.zeros(np.size(z))
c2o3 = np.zeros(np.size(z))
hcho = np.zeros(np.size(z))
rooh = np.zeros(np.size(z))
rh  = np.zeros(np.size(z))
ro2  = np.zeros(np.size(z))
ho2 = np.zeros(np.size(z))
o3 = np.zeros(np.size(z))
no = np.zeros(np.size(z))
no3 = np.zeros(np.size(z))
no2 = np.zeros(np.size(z))
oh = np.zeros(np.size(z))
co2 = np.zeros(np.size(z))

for k in range(kmax):

    # potential temperature: for now assume no liquid water:
    if (z[k] < 352.5):
        thl[k] = 296.6
        qt[k] = 11.8
    elif (z[k] < 442.5):
        thl[k] = 296.6 + 1.5*(z[k]-352.5)/(442.5-352.5)
        qt[k] = 11.8 - 4.0*(z[k]-352.5)/(442.5-352.5)
    else:
        thl[k] = 298.1 + (z[k] - 443.5)*0.003
        qt[k] = max(0.0,7.8 - (z[k]- 443.5 )*0.004)


    # Subsidence
    if(z[k] < 2260):
        wls[k] = -0.005 * (z[k] / 2260.)
    else:
        wls[k] = -0.005

    # U and V component wind
    u[k]  = u0
    ug[k] = u0
    v[k]  = v0
    vg[k] = v0
    co2[k] = co20

    thlls[k] = 0.0
    qtls[k] = 0.0

    # Advective tendency qt
#    qtls[k] = 1.5e-4
    # Advective and radiative tendency thl
#    thlls[k] = 6.4e-4   # K/s

    # Advective tendency qt
#    qtls[k] = 1.5e-4

# normalize profiles to SI
qt  /= 1000.
qtls/= 1000.

td = 19.5 - 6.0
time_surface = np.arange(13*4+1)*0.25
H = 0.1*np.sin(np.pi*(time_surface-1.0)/td)
idx = np.where(H < 0.0)
H[idx] = 0.0
td = 19.5 - 7.0
sbotqt = 0.15*np.sin(np.pi*(time_surface-2.0)/td)
idx = np.where(sbotqt < 0.0)
sbotqt[idx] = 0.0

# Calculate the surface fluxes in the correct units.
Rd  = 287.
cp  = 1005.
Lv  = 2.5e6
p0  = 1e5
rho = p0/(Rd*thl[0]*(1. + 0.61*qt[0]))
time_surface *= 3600. # h to s
sbotthl = H    #  /(rho*cp)
sbotqt  *= 1e-3    # kg/kg 



for k in range(kmax):

    # potential temperature: for now assume no liquid water:
    if (z[k] <= 1000):
        hno3[k] = 2.4
        n2o5[k] = 0.0
        h2o2[k] = 2.0
        co[k] = 213.0 
        pan[k] = 0.0
        c2o3[k] = 0.0
        hcho[k] = 4.2 
        rooh[k] = 2.7
        rh[k] = 1.0
        ro2[k] = 0.047
        ho2[k] = 0.071
        o3[k] = 69.5
        no[k] = 0.081
        no3[k] = 0.0
        no2[k] = 0.52
        oh[k]  = 0.0
    else:
        hno3[k] = 2.4
        n2o5[k] = 0.0
        h2o2[k] = 2.0
        co[k] = 213.0 
        pan[k] = 0.0
        c2o3[k] = 0.0
        hcho[k] = 4.2 
        rooh[k] = 2.7
        rh[k] = 1.0
        ro2[k] = 0.047
        ho2[k] = 0.071
        o3[k] = 69.5
        no[k] = 0.081
        no3[k] = 0.0
        no2[k] = 0.52
        oh[k]  = 0.0

# model time dependent parameters:
with open('input_chem_jvals') as f:
    line = f.readline()
    line = f.readline()
    line = f.readline()
    line = f.readline()
    header = line.split(' ')
    lines = f.readlines()
vals = []
for line in lines:
    vals.append([float(i) for i in line.split()])
vals = np.array(vals)
time_chem = vals[:,0]
# now the run starts at 5 AM: so we have to make sure that t = 0 corresponds to 5 AM.
tstart = 5.0
idx = abs(time_chem-tstart).argmin()
vals = vals[idx:,:]
time_chem = (vals[:,0]-tstart)*3600.0   # start at zero
jo31d = vals[:,2]
jh2o2= vals[:,3]
jno2= vals[:,4]
jno3a= vals[:,5]
jno3b= vals[:,6]
jch2or= vals[:,7]
jch2om= vals[:,8]
jch3o2h= vals[:,9]

# isoprene emissions:
td = 19. - 8.0
emi_isop = 1.1*np.sin(np.pi*(time_chem/3600.-3.0)/td)
idx = np.where(emi_isop < 0.0)
emi_isop[idx] = 0.0
# no emissions:
td = 19. - 8.0
emi_no = 0.03*np.sin(np.pi*(time_chem/3600.-3.0)/td)
idx = np.where(emi_no < 0.0)
emi_no[idx] = 0.0




"""
# well mixed profile with jump
h    = 1000.
dth  = 10.
dthz = 100.
dthetadz = 0.003

for k in range(kmax):
    if(z[k] <= h - 0.5*dthz):
        th[k] = 300.
    elif(z[k] <= h + 0.5*dthz):
        th[k] = 300. + dth/dthz * (z[k]-(h-0.5*dthz))
    else:
        th[k] = 300. + dth + dthetadz*(z[k]-(h+0.5*dthz))
    thls[k] = 2.*(z[k]/zsize - 0.5) / 3600.
    wls [k] = -0.01*(z[k]/zsize)
"""

nc_file = nc.Dataset("tm5_ifs_input.nc", mode="w", datamodel="NETCDF4", clobber=False)

nc_file.createDimension("z", kmax)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_thl   = nc_group_init.createVariable("thl"   , float_type, ("z"))
nc_qt    = nc_group_init.createVariable("qt"    , float_type, ("z"))
nc_u     = nc_group_init.createVariable("u"     , float_type, ("z"))
nc_ug    = nc_group_init.createVariable("u_geo"    , float_type, ("z"))
nc_v     = nc_group_init.createVariable("v"     , float_type, ("z"))
nc_vg    = nc_group_init.createVariable("v_geo"    , float_type, ("z"))
nc_wls   = nc_group_init.createVariable("w_ls"  , float_type, ("z"))
nc_thlls = nc_group_init.createVariable("thl_ls", float_type, ("z"))
nc_qtls  = nc_group_init.createVariable("qt_ls" , float_type, ("z"))

nc_z    [:] = z    [:]
nc_thl  [:] = thl  [:]
nc_qt   [:] = qt   [:]
nc_u    [:] = u    [:]
nc_ug   [:] = ug   [:]
nc_v    [:] = v    [:]
nc_vg   [:] = vg   [:]
nc_wls  [:] = wls  [:]
nc_thlls[:] = thlls[:]
nc_qtls [:] = qtls [:]

nc_hno3 = nc_group_init.createVariable("hno3", float_type, ("z"))
nc_n2o5 = nc_group_init.createVariable("n2o5", float_type, ("z"))
nc_h2o2 = nc_group_init.createVariable("h2o2", float_type, ("z"))
nc_co = nc_group_init.createVariable("co", float_type, ("z"))
nc_pan = nc_group_init.createVariable("pan", float_type, ("z"))
nc_c2o3 = nc_group_init.createVariable("c2o3", float_type, ("z"))
nc_hcho = nc_group_init.createVariable("hcho", float_type, ("z"))
nc_rooh = nc_group_init.createVariable("rooh", float_type, ("z"))
nc_rh = nc_group_init.createVariable("rh", float_type, ("z"))
nc_ro2 = nc_group_init.createVariable("ro2", float_type, ("z"))
nc_ho2 = nc_group_init.createVariable("ho2", float_type, ("z"))
nc_o3 = nc_group_init.createVariable("o3", float_type, ("z"))
nc_no2 = nc_group_init.createVariable("no2", float_type, ("z"))
nc_no3 = nc_group_init.createVariable("no3", float_type, ("z"))
nc_no = nc_group_init.createVariable("no", float_type, ("z"))
nc_oh = nc_group_init.createVariable("oh", float_type, ("z"))
nc_co2 = nc_group_init.createVariable("co2", float_type, ("z"))

nc_hno3[:] = hno3[:]
nc_n2o5[:] = n2o5[:] 
nc_h2o2[:] = h2o2[:] 
nc_co[:] = co[:] 
nc_pan[:] = pan[:] 
nc_c2o3[:] = c2o3[:] 
nc_hcho[:] = hcho[:] 
nc_rooh[:] = rooh[:] 
nc_rh[:] = rh[:] 
nc_ro2[:] = ro2[:] 
nc_ho2[:] = ho2[:] 
nc_o3[:] = o3[:] 
nc_no2[:] = no2[:] 
nc_no3[:] = no3[:] 
nc_no[:] = no[:] 
nc_oh[:] = oh[:] 
nc_co2[:] = co2[:] 


# Create a group called "timedep" for the timedep.
nc_group_timedep = nc_file.createGroup("timedep")
nc_group_timedep.createDimension("time_surface", time_surface.size)
nc_time_surface = nc_group_timedep.createVariable("time_surface", float_type, ("time_surface"))
nc_thl_sbot = nc_group_timedep.createVariable("thl_sbot", float_type, ("time_surface"))
nc_qt_sbot  = nc_group_timedep.createVariable("qt_sbot" , float_type, ("time_surface"))
nc_time_surface[:] = time_surface[:]
nc_thl_sbot    [:] = sbotthl     [:]
nc_qt_sbot     [:] = sbotqt      [:]

# Create a group called "timedep_chem" for the chemistry time dependent emmission/photolysis
nc_group_timedepchem = nc_file.createGroup("timedep_chem")
nc_group_timedepchem.createDimension("time_chem", time_chem.size)
nc_time_chem = nc_group_timedepchem.createVariable("time_chem", float_type, ("time_chem"))
nc_jo31d = nc_group_timedepchem.createVariable("jo31d", float_type, ("time_chem"))
nc_jh2o2 = nc_group_timedepchem.createVariable("jh2o2", float_type, ("time_chem"))
nc_jno2 = nc_group_timedepchem.createVariable("jno2", float_type, ("time_chem"))
nc_jno3a = nc_group_timedepchem.createVariable("jno3a", float_type, ("time_chem"))
nc_jno3b = nc_group_timedepchem.createVariable("jno3b", float_type, ("time_chem"))
nc_jch2or = nc_group_timedepchem.createVariable("jch2or", float_type, ("time_chem"))
nc_jch2om = nc_group_timedepchem.createVariable("jch2om", float_type, ("time_chem"))
nc_jch3o2h = nc_group_timedepchem.createVariable("jch3o2h", float_type, ("time_chem"))
nc_emi_isop = nc_group_timedepchem.createVariable("emi_isop", float_type, ("time_chem"))
nc_emi_no = nc_group_timedepchem.createVariable("emi_no", float_type, ("time_chem"))
nc_time_chem[:] = time_chem[:]
nc_jo31d[:] = jo31d [:]
nc_jh2o2[:] = jh2o2 [:]
nc_jno2[:] = jno2 [:]
nc_jno3a[:] = jno3a [:]
nc_jno3b[:] = jno3b [:]
nc_jch2or[:] = jch2or [:]
nc_jch2om[:] = jch2om [:]
nc_jch3o2h[:] = jch3o2h [:]
nc_emi_isop[:] = emi_isop[:]
nc_emi_no[:] = emi_no[:]

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


# some constants for thl:
cp = 1005.7   # [J kg-1 K]   specific heat dry air
Lv = 2.501e6    # [J kg-1]     latent heat vaporization

