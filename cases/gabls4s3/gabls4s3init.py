import numpy as np
import pylab as pl
import netCDF4 as nc4

# Read the stage 3 driver file
class read_driver:
    def __init__(self):
        nc = nc4.Dataset('GABLS4_SCM_LES_STAGE3.nc', 'r')
        self.t = nc.variables["time"][:]
        self.z = nc.variables["height"][:][::-1]
        self.p = nc.variables["pf"][:][::-1]

        # initial profiles
        self.th = nc.variables["theta"][:][::-1]
        self.T  = nc.variables["t"][:][::-1]
        self.q  = nc.variables["qv"][:][::-1] # = zero
        self.u  = nc.variables["u"][:][::-1]
        self.v  = nc.variables["v"][:][::-1]

        # time varying forcings
        self.ug   = nc.variables["Ug"][0,:][::-1] # u geo wind = constant in time
        self.vg   = nc.variables["Vg"][0,:][::-1] # v geo wind = consant in time
        self.advT = nc.variables["hadvT"][0,:][::-1] # = zero
        self.advq = nc.variables["hadvQ"][0,:][::-1] # = zero
        self.Ts   = nc.variables["Tg"][:]

        self.ps  = nc.variables["psurf"][0]
        self.z0m = nc.variables["z0m"][0]

        # Calculate theta_s
        self.ths = self.Ts / (self.ps / 1.e5)**(287.04/1005.)

s3 = read_driver()

# Get number of vertical levels and size from .ini file
with open('gabls4s3.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

print('kmax = %i, zsize = %.1f m'%(kmax, zsize))

# set the height
dz = zsize / kmax
z  = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
th = np.zeros(np.size(z))
u  = np.zeros(np.size(z))
ug = np.zeros(np.size(z))
v  = np.zeros(np.size(z))
vg = np.zeros(np.size(z))

th = np.interp(z, s3.z, s3.th)
u  = np.interp(z, s3.z, s3.u)
v  = np.interp(z, s3.z, s3.v)
ug = np.interp(z, s3.z, s3.ug)
vg = np.interp(z, s3.z, s3.vg)

# write the data to a file
proffile = open('gabls4s3.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s} {5:^20s}\n'.format('z','th','u','ug','v','vg'))
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E} \n'.format(z[k], th[k], u[k], ug[k], v[k], vg[k]))
proffile.close()

# write surface temperature
timefile = open('gabls4s3.time','w')
timefile.write('{0:^20s} {1:^20s} \n'.format('t','sbot[th]'))
for t in range(s3.t.size):
    timefile.write('{0:1.14E} {1:1.14E} \n'.format(s3.t[t], s3.ths[t]))
timefile.close()

"""
# Plot
pl.close('all')

pl.figure()
pl.subplot(221)
pl.plot(th, z, 'k-', label='mhh')
pl.plot(s3.th, s3.z, 'go', mfc='none', label='s3')
pl.ylim(0,1100)
pl.xlim(270,285)
pl.legend(frameon=False, loc=2)

pl.subplot(222)
pl.plot(u, z, 'k-', label='mhh')
pl.plot(s3.u, s3.z, 'go', mfc='none', label='s3')
pl.plot(ug, z, 'k--', label='mhh')
pl.plot(s3.ug, s3.z, 'bo', mfc='none', label='s3')
pl.ylim(0,1100)
pl.xlim(0,10)
pl.legend(frameon=False, loc=2)

pl.subplot(223)
pl.plot(v, z, 'k-', label='mhh')
pl.plot(s3.v, s3.z, 'go', mfc='none', label='s3')
pl.plot(vg, z, 'k--', label='mhh')
pl.plot(s3.vg, s3.z, 'bo', mfc='none', label='s3')
pl.ylim(0,1100)
pl.xlim(0,10)
pl.legend(frameon=False, loc=2)

pl.subplot(224)
pl.plot(s3.t, s3.ths, 'k-', label='s3')
pl.legend(frameon=False, loc=2)
"""
