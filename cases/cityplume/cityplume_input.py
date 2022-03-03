import matplotlib.pyplot as pl
import numpy as np
import netCDF4 as nc

float_type = "f8"
# float_type = "f4"

np_dtype = np.float64 if float_type == "f8" else np.float32

# Get number of vertical levels and size from .ini file
with open('cityplume.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])
        if(line.split('=')[0]=='itot'):
            itot = int(line.split('=')[1])
        if(line.split('=')[0]=='jtot'):
            jtot = int(line.split('=')[1])

dz = zsize / kmax

dthetadz = 0.003

# set the height
z  = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
u  = np.ones(np.size(z))*8
v  = np.zeros(np.size(z))
thl = np.zeros(np.size(z))
qt    = np.zeros(z.size)

hno3 = np.zeros(z.size)
co = np.zeros(z.size)
hcho = np.zeros(z.size)
rooh = np.zeros(z.size)
h2o2 = np.zeros(z.size)
c3h8 = np.zeros(z.size)
no2 = np.zeros(z.size)
no = np.zeros(z.size)
o3 = np.zeros(z.size)
ro2 = np.zeros(z.size)
ho2 = np.zeros(z.size)
oh = np.zeros(z.size)
no3 = np.zeros(z.size)
n2o5 = np.zeros(z.size)


# linearly stratified profile
for k in range(kmax):
    co[k] = 90.0
    h2o2[k] = 1.0
    hno3[k] = 1.0
    hcho[k] = 1.0
    rooh[k] = 1.0
    h2o2[k] = 1.0
    c3h8[k] = 2.0
    no2[k] = 0.10
    no[k] = 0.05

#  from the orlando case:
    if (z[k] < 352.5):
        thl[k] = 296.6
        qt[k] = 11.8
        o3[k] = 30.0
    elif (z[k] < 442.5):
        c3h8[k] = 1.0
        thl[k] = 296.6 + 1.5*(z[k]-352.5)/(442.5-352.5)
        qt[k] = 11.8 - 4.0*(z[k]-352.5)/(442.5-352.5)
        o3[k] = 30.0 + 4.0*(z[k]-352.5)/(442.5-352.5)

    else:
        c3h8[k] = 0.5
        o3[k] = 34.0 + (z[k] - 443.5)*0.005  # roughly 50 aloft
        thl[k] = 298.1 + (z[k] - 443.5)*0.003
        qt[k] = max(0.0,7.8 - (z[k]- 443.5 )*0.004)

        
qt *= 1e-3  # kg/kg

nc_file = nc.Dataset("cityplume_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", kmax)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_u  = nc_group_init.createVariable("u" , float_type, ("z"))
nc_v  = nc_group_init.createVariable("v" , float_type, ("z"))
nc_thl = nc_group_init.createVariable("thl", float_type, ("z"))
nc_qt    = nc_group_init.createVariable("qt"    , float_type, ("z"))
nc_co = nc_group_init.createVariable("co", float_type, ("z"))
nc_hno3 = nc_group_init.createVariable("hno3", float_type, ("z"))
nc_hcho = nc_group_init.createVariable("hcho", float_type, ("z"))
nc_rooh = nc_group_init.createVariable("rooh", float_type, ("z"))
nc_h2o2 = nc_group_init.createVariable("h2o2", float_type, ("z"))
nc_c3h8 = nc_group_init.createVariable("c3h8", float_type, ("z"))
nc_no2 = nc_group_init.createVariable("no2", float_type, ("z"))
nc_no = nc_group_init.createVariable("no", float_type, ("z"))
nc_o3 = nc_group_init.createVariable("o3", float_type, ("z"))
nc_ro2 = nc_group_init.createVariable("ro2", float_type, ("z"))
nc_ho2 = nc_group_init.createVariable("ho2", float_type, ("z"))
nc_oh = nc_group_init.createVariable("oh", float_type, ("z"))
nc_no3 = nc_group_init.createVariable("no3", float_type, ("z"))
nc_n2o5 = nc_group_init.createVariable("n2o5", float_type, ("z"))

nc_z [:] = z [:]
nc_u [:] = u [:]
nc_v [:] = v [:]
nc_thl[:] = thl[:]
nc_qt[:] = qt[:]
# No loger in ppb ---> mixing ratios
nc_hno3[:] = hno3[:]*1e-9
nc_co[:] = co[:]*1e-9 
nc_hcho[:] = hcho[:]*1e-9 
nc_rooh[:] = rooh[:]*1e-9 
nc_h2o2[:] = h2o2[:]*1e-9 
nc_c3h8[:] = c3h8[:]*1e-9 
nc_no2[:] = no2[:]*1e-9 
nc_no[:] = no[:]*1e-9 
nc_o3[:] = o3[:]*1e-9 
nc_ro2[:] = ro2[:]*1e-9 
nc_ho2[:] = ho2[:]*1e-9 
nc_oh[:] = oh[:]*1e-9 
nc_no3[:] = no3[:]*1e-9 
nc_n2o5[:] = n2o5[:]*1e-9 

# model time dependent parameters: now from TUV calculation for Ryad...
# Create time varying surface fields.   Simulate 8 hours: from 6 AM --> 2 PM (1:30 TROPOMI overpass). 8*5
# latitude and longitude of Riyadh is 24.63° N, 46.71°E. 
with open('riyadh.txt') as f:
    # output is in order: O3, H2O2, NO2, NO3, N2O5, CH2OR, CH2OM, ROOH
    line = 'xxx'
    while not line.startswith('time'):
        line = f.readline()
    lines = f.readlines()
vals = []
for line in lines[:-1]:
    vals.append([float(i) for i in line.split()])
vals = np.array(vals)
time_chem = vals[:,0]
tstart = 6.0
idx = abs(time_chem-tstart).argmin()
vals = vals[idx:,:]
time_chem = (vals[:,0]-tstart)*3600.0   # start at zero
print(time_chem,' chemistry time')
jo31d = vals[:,2]
jh2o2= vals[:,3]
jno2= vals[:,4]
jno3= vals[:,5]
jn2o5= vals[:,6]
jch2or= vals[:,7]
jch2om= vals[:,8]
jch3o2h= vals[:,9]

# these emissions are set to zero: do not worry yet.
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
# Create a group called "timedep_chem" for the chemistry time dependent emmission/photolysis
nc_group_timedepchem = nc_file.createGroup("timedep_chem")
nc_group_timedepchem.createDimension("time_chem", time_chem.size)
nc_time_chem = nc_group_timedepchem.createVariable("time_chem", float_type, ("time_chem"))
nc_jo31d = nc_group_timedepchem.createVariable("jo31d", float_type, ("time_chem"))
nc_jh2o2 = nc_group_timedepchem.createVariable("jh2o2", float_type, ("time_chem"))
nc_jno2 = nc_group_timedepchem.createVariable("jno2", float_type, ("time_chem"))
nc_jno3 = nc_group_timedepchem.createVariable("jno3", float_type, ("time_chem"))
nc_jn2o5 = nc_group_timedepchem.createVariable("jn2o5", float_type, ("time_chem"))
nc_jch2or = nc_group_timedepchem.createVariable("jch2or", float_type, ("time_chem"))
nc_jch2om = nc_group_timedepchem.createVariable("jch2om", float_type, ("time_chem"))
nc_jch3o2h = nc_group_timedepchem.createVariable("jch3o2h", float_type, ("time_chem"))
nc_emi_isop = nc_group_timedepchem.createVariable("emi_isop", float_type, ("time_chem"))
nc_emi_no = nc_group_timedepchem.createVariable("emi_no", float_type, ("time_chem"))
nc_time_chem[:] = time_chem[:]
nc_jo31d[:] = jo31d [:]
nc_jh2o2[:] = jh2o2 [:]
nc_jno2[:] = jno2 [:]
nc_jno3[:] = jno3 [:]
nc_jn2o5[:] = jn2o5 [:]
nc_jch2or[:] = jch2or [:]
nc_jch2om[:] = jch2om [:]
nc_jch3o2h[:] = jch3o2h [:]
nc_emi_isop[:] = emi_isop[:]
nc_emi_no[:] = emi_no[:]

nc_file.close()

# Create time varying surface fields.   Simulate 8 hours: from 6 AM --> 2 PM (1:30 TROPOMI overpass). 8*5
# generate time series of surface heat flux from 6AM to 6PM
t_surface = 6.0 + np.arange(24+1)*0.5
H = 0.1*np.sin(np.pi*(t_surface-6.0)/12.)
sbotqt = 0.15*np.sin(np.pi*(t_surface-6.0)/12.)
sbotqt/=1e3
endtime = 12*3600
dt = 1800
nt = int((endtime / dt)+1)

thl_sbot = np.zeros((nt, jtot, itot), dtype=np_dtype)
co_sbot = np.zeros((nt, jtot, itot), dtype=np_dtype)
c3h8_sbot = np.zeros((nt, jtot, itot), dtype=np_dtype)
no_sbot = np.zeros((nt, jtot, itot), dtype=np_dtype)
qt_sbot = np.zeros((nt, jtot, itot), dtype=np_dtype)
# js = 
for ii,hh in enumerate(H):
    thl_sbot[ii,:,:] = hh
    thl_sbot[ii, 30:130, 5:105] = hh
# calculate emissions based on Srijana's numbers:
# co emissions  11.8 kg/s 
# c3h8 emissions: assume 10 kg/s
# no emissions: 9 kg/s (as NO2)
area = 30e3*30e3  # area of assumed city
e_co = 11.8/area  # kg/m2/s
e_no = 9.0/area
e_c3h8 = 0.0/area   # set to zero for this run
rhoz = 1.1662 # kg/m3
xmair = 28.9647 # kg/kmol
xmco = 28.0
xmno = 46.0   # as mass NO2
xmc3h8 = 42.0   # arbitrary, but propene here 
conv = rhoz/xmair # kmol/m3 

co_sbot[ii,:,:] = 0.0
co_sbot[:, 30:130, 5:105] = e_co/(conv*xmco)  # kg/m2s ---> kmolCO/m2s [/kmol-air/m3 = xmair ] --> mixing ratio change m/s
c3h8_sbot[ii,:,:] = 0.0
c3h8_sbot[:, 30:130, 5:105] = e_c3h8/(conv*xmc3h8)
no_sbot[ii,:,:] = 0.0
no_sbot[:, 30:130, 5:105] = e_no/(conv*xmno)

for ii,hh in enumerate(sbotqt):
    qt_sbot[:] = hh
    qt_sbot[ii, 30:130, 5:105] = hh

# Write as binary input files for MicroHH
for t in range(nt):
    print(thl_sbot[t,:].mean(), thl_sbot[t,:].max(), thl_sbot[t,:].min())
    thl_sbot[t,:].tofile('thl_flux.{0:07d}'.format(t*dt))
    co_sbot[t,:].tofile('co_flux.{0:07d}'.format(t*dt))
    c3h8_sbot[t,:].tofile('c3h8_flux.{0:07d}'.format(t*dt))
    no_sbot[t,:].tofile('no_flux.{0:07d}'.format(t*dt))
    qt_sbot[t,:].tofile('qt_flux.{0:07d}'.format(t*dt))
