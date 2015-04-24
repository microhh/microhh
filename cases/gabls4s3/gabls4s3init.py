import numpy as np
import pylab as pl
import netCDF4 as nc4

# load the init script to get variables like ug, vg, ts
from gabls4s3init import *

pl.close('all')

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

        self.ps  = nc.variables["psurf"].getValue()
        self.z0m = nc.variables["z0m"].getValue()

        # Calculate theta_s
        self.ths = self.Ts / (self.ps / 1.e5)**(287.04/1005.)

s3 = read_driver()

# ----------------
# non-stretched grid
# ----------------
if(False):
   # Get number of vertical levels and size from .ini file
   with open('gabls4s3.ini') as f:
       for line in f:
           if(line.split('=')[0]=='ktot'):
               kmax = int(line.split('=')[1])
           if(line.split('=')[0]=='zsize'):
               zsize = float(line.split('=')[1])
   
   print('kmax = %i, zsize = %.1f m'%(kmax, zsize))
   dz = zsize / kmax
   z  = np.linspace(0.5*dz, zsize-0.5*dz, kmax)

# ----------------
# stretched grid
# ----------------
if(True):
    kmax  = 500
    zsize = 1000
    dz1   = 2.
    dz2   = 40.
    z1    = 450
    z2    = 600
    z3    = 700
    fac   = (z3 - z1) / 10.
    
    dz0   = zsize / kmax
    z0    = np.linspace(0.5 * dz0, zsize - 0.5*dz0, kmax)
    dz    = dz1 + (dz2 - dz1) / (1. + np.exp(-(z0 - z2) / fac))
    
    z = [0.5 * dz1]
    for k in range(kmax):
        if(z[-1] + dz[k] >= zsize):
            break
        z.append(z[-1] + dz[k])
    
    zsize = z[-1] + 0.5*(z[-1] - z[-2])
    z = np.array(z)
    kmax  = np.size(z) 

print('kmax = %i, zsize = %f'%(kmax, zsize))

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

# Plot
if(True):
    zh_fleur = np.loadtxt('grille_stretche')
    z_fleur  = 0.5 * (zh_fleur[1:] + zh_fleur[:-1])
    dz_fleur = zh_fleur[1:] - zh_fleur[:-1]

    pl.figure()
    pl.plot(dz[:kmax], z, 'k-o', mfc='none')
    pl.plot(dz_fleur, z_fleur, 'g-x', mfc='none')
    pl.xlabel('dz [m]')
    pl.ylabel('z [m]')

if(False):
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
