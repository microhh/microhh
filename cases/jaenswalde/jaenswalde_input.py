import matplotlib.pyplot as pl
import numpy as np
import netCDF4 as nc

float_type = "f8"

Mair = 28.9664
Mco2 = 44.01

# Get number of vertical levels and size from .ini file
with open('jaenswalde.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])
        if(line.split('=')[0]=='ysize'):
            ysize = float(line.split('=')[1])

dz = zsize / kmax

z = np.arange(0.5*dz, zsize, dz)
u = np.zeros(z.size)
v = np.zeros(z.size)
th = np.zeros(z.size)
thl   = np.zeros(z.size)
qt    = np.zeros(z.size)
thlls = np.zeros(z.size)
qtls  = np.zeros(z.size)
co2 = np.zeros(z.size)
hno3 = np.zeros(z.size)
co = np.zeros(z.size)
hcho = np.zeros(z.size)
noy = np.zeros(z.size)
rooh =np.zeros(z.size)
h2o2 = np.zeros(z.size)
rh = np.zeros(z.size)
no2 = np.zeros(z.size)
no = np.zeros(z.size)
o3 = np.zeros(z.size)
ro2 = np.zeros(z.size)
ho2 = np.zeros(z.size)
oh = np.zeros(z.size)

# well mixed profile with jump
th0 = 280.          # Bulk potential temperature
dth = 10.            # Potential temperature inversion
dthz = 300          # Inversion thickness
dthetadz = 0.006    # Free troposphere lapse rate
zi = 1200.          # Boundary layer depth
u0 = 5.             # Zonal wind speed
v0 = 0.             # Meridional wind speed
co20 = 400.0         # CO2 concentration (ppm)

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

for k in range(kmax):
    if(z[k] <= zi - 0.5*dthz):
        th[k] = th0
        thl[k] = th0
        qt[k] = 6.0
    elif(z[k] <= zi + 0.5*dthz):
        th[k] = th0 + dth/dthz * (z[k]-(zi-0.5*dthz))
        thl[k] = th0 + dth/dthz * (z[k]-(zi-0.5*dthz))
        qt[k] = 6.0 - 2.0*(z[k]-352.5)/(442.5-352.5)
    else:
        th[k] = th0 + dth + dthetadz*(z[k]-(zi+0.5*dthz))
        thl[k] = th0 + dth + dthetadz*(z[k]-(zi+0.5*dthz))
        qt[k] = max(0.0,3.0 - (z[k]- 443.5 )*0.004)
    # potential temperature: for now assume no liquid water:

    if (z[k] <= 1000):
        hno3[k] = 2.4
        h2o2[k] = 8.8
        o3[k] = 69.5
        co[k] = 213.0 
        no[k] = 0.081
        no2[k] = 0.52
        hcho[k] = 4.2 
        rooh[k] = 2.7 
        rh[k]  = 6.8 
        ho2[k] = 0.071
        ro2[k] = 0.047
    else:
        hno3[k] = 2.4
        h2o2[k] = 8.8
        o3[k] = 69.5
        co[k] = 213.0 
        no[k] = 0.081
        no2[k] = 0.52
        hcho[k] = 4.2 
        rooh[k] = 2.7 
        rh[k]  = 6.8 
        ho2[k] = 0.071
        ro2[k] = 0.047

u[:] = u0
v[:] = v0
co2[:] = co20 * 1e3    # now in ppb
# Calculate the surface fluxes in the correct units.
Rd  = 287.
cp  = 1005.
Lv  = 2.5e6
p0  = 1e5
rho = p0/(Rd*thl[0]*(1. + 0.61*qt[0]))
time_surface *= 3600. # h to s
sbotthl = H    #  /(rho*cp)
sbotqt  *= 1e-3    # kg/kg 



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



# Write input NetCDF file
nc_file = nc.Dataset("jaenswalde_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z" , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_u = nc_group_init.createVariable("u" , float_type, ("z"))
nc_v = nc_group_init.createVariable("v" , float_type, ("z"))
nc_th = nc_group_init.createVariable("th", float_type, ("z"))
nc_thl   = nc_group_init.createVariable("thl"   , float_type, ("z"))
nc_qt    = nc_group_init.createVariable("qt"    , float_type, ("z"))
nc_co2 = nc_group_init.createVariable("co2", float_type, ("z"))
nc_hno3 = nc_group_init.createVariable("hno3", float_type, ("z"))
nc_co = nc_group_init.createVariable("co", float_type, ("z"))
nc_hcho = nc_group_init.createVariable("hcho", float_type, ("z"))
nc_noy = nc_group_init.createVariable("noy", float_type, ("z"))
nc_rooh = nc_group_init.createVariable("rooh", float_type, ("z"))
nc_h2o2 = nc_group_init.createVariable("h2o2", float_type, ("z"))
nc_rh = nc_group_init.createVariable("rh", float_type, ("z"))
nc_no2 = nc_group_init.createVariable("no2", float_type, ("z"))
nc_no = nc_group_init.createVariable("no", float_type, ("z"))
nc_o3 = nc_group_init.createVariable("o3", float_type, ("z"))
nc_ro2 = nc_group_init.createVariable("ro2", float_type, ("z"))
nc_ho2 = nc_group_init.createVariable("ho2", float_type, ("z"))
nc_oh = nc_group_init.createVariable("oh", float_type, ("z"))

nc_z [:] = z[:]
nc_u [:] = u[:]
nc_v [:] = v[:]
nc_th[:] = th[:]
nc_thl  [:] = thl  [:]
nc_qt   [:] = qt   [:]
nc_co2[:] = co2[:]
nc_hno3 [:] = hno3[:]
nc_co [:] = co[:]
nc_hcho [:] = hcho[:]
nc_rooh [:] = rooh[:]
nc_h2o2 [:] = h2o2[:]
nc_rh [:] = rh[:]
nc_no2 [:] = no2[:]
nc_no [:] = no[:]
nc_o3 [:] = o3[:]
nc_ro2 [:] = ro2[:]
nc_ho2 [:] = ho2[:]
nc_oh [:] = oh[:]

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
nc_emi = nc_group_timedepchem.createVariable("emi_isop", float_type, ("time_chem"))
nc_time_chem[:] = time_chem[:]
nc_jo31d[:] = jo31d [:]
nc_jh2o2[:] = jh2o2 [:]
nc_jno2[:] = jno2 [:]
nc_jno3a[:] = jno3a [:]
nc_jno3b[:] = jno3b [:]
nc_jch2or[:] = jch2or [:]
nc_jch2om[:] = jch2om [:]
nc_jch3o2h[:] = jch3o2h [:]
nc_emi[:] = emi_isop[:]

nc_file.close()

# Print .ini settings emissions:
# Coordinates of central cooling tower (m):
x0 = 1000
y0 = ysize/2.
z0 = 150

# x,y spacing towers:
dx = 290
dy = 120
ddx = 40

# Std-dev of plume widths:
sigma_x = 50
sigma_y = 50
sigma_z = 50

# Strength of plumes:
# CO2 emissions Jaenswalde, from:
# https://wwfeu.awsassets.panda.org/downloads/european_dirty_thirty_may_2007.pdf,
# which is probably not the best source.....

xmco2 = 12.0 + 2*16.0 # kg/kmol
xmno  = 14.0 + 2*16.0  #  is given in mass NO2. kg NO2/kmol
strength_co2 = 23.7e9 / 365.25 / 24. / 3600. / 9. / xmco2
strength_co2 = np.round(strength_co2, decimals=3)
strength_no = 19e6 / 365.25 / 24. / 3600. / 9. / xmno    # strength is in kmol/s.
strength_no = np.round(strength_no, decimals=6)

# Generate lists with model input:
source_x0 = []
source_y0 = []
source_z0 = []

for tracer in range(2):
    for j in range(-1,2):
        for i in range(-1,2):
            source_x0.append( x0 + i*dx + j*ddx )
            source_y0.append( y0 + j*dy )
            source_z0.append( z0 )

def to_cs_list(lst):
    """ From Python list to comma separated string """
    lst = [str(x) for x in lst]
    return ','.join(lst)

def constant_list(value, n):
    """ Create comma separated list with constant values """
    lst = n*[value]
    return to_cs_list(lst)

print('sourcelist={},{}'.format(constant_list('co2', 9),constant_list('no',9)))

print('source_x0={}'.format(to_cs_list(source_x0)))
print('source_y0={}'.format(to_cs_list(source_y0)))
print('source_z0={}'.format(to_cs_list(source_z0)))

print('sigma_x={}'.format(constant_list(sigma_x, 18)))
print('sigma_y={}'.format(constant_list(sigma_y, 18)))
print('sigma_z={}'.format(constant_list(sigma_z, 18)))

print('strength={},{}'.format(constant_list(strength_co2, 9),constant_list(strength_no,9)))

print('line_x={}'.format(constant_list(0, 18)))
print('line_y={}'.format(constant_list(0, 18)))
print('line_z={}'.format(constant_list(0, 18)))
