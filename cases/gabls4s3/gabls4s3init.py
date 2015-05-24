import numpy as np
import pylab as pl
import netCDF4 as nc4

# load the init script to get variables like ug, vg, ts
from gabls4s3init import *

pl.close('all')

# Create stretched grid
class Grid:
    def __init__(self, kmax, nloc1, nbuf1,dz1, dz2):
        dn         = 1./kmax
        n          = np.linspace(dn, 1.-dn, kmax)
        nloc1     *= dn
        nbuf1     *= dn
        dzdn1      = dz1/dn
        dzdn2      = dz2/dn

        dzdn       = dzdn1 + 0.5*(dzdn2-dzdn1)*(1. + np.tanh((n-nloc1)/nbuf1))
        self.dz    = dzdn*dn

        self.kmax  = kmax
        self.z     = np.zeros(self.dz.size)
        stretch    = np.zeros(self.dz.size)

        self.z[0]  = 0.5*self.dz[0]
        stretch[0] = 1.

        for k in range(1, self.kmax):
              self.z [k] = self.z[k-1] + 0.5*(self.dz[k-1]+self.dz[k])
              stretch[k] = self.dz[k]/self.dz[k-1]

        zsize = self.z[kmax-1] + 0.5*self.dz[kmax-1]
        print('kmax=%i, zsize=%f\n'%(kmax,zsize))

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

outname = 'gabls4s3_restart'

s3 = read_driver()

#g = Grid(300, 250, 20, 2, 10) # 2 m 
#g = Grid(550, 500, 20, 1, 10) # 1 m
#g = Grid(80, 50, 10, 10, 20) # 10 m
g = Grid(144, 110, 10, 1, 3) # restart grid

th = np.zeros(g.z.size)
u  = np.zeros(g.z.size)
ug = np.zeros(g.z.size)
v  = np.zeros(g.z.size)
vg = np.zeros(g.z.size)

th = np.interp(g.z, s3.z, s3.th)
u  = np.interp(g.z, s3.z, s3.u)
v  = np.interp(g.z, s3.z, s3.v)
ug = np.interp(g.z, s3.z, s3.ug)
vg = np.interp(g.z, s3.z, s3.vg)

# write the data to a file
proffile = open(outname+'.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s} {5:^20s}\n'.format('z','th','u','ug','v','vg'))
for k in range(g.kmax):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E} \n'.format(g.z[k], th[k], u[k], ug[k], v[k], vg[k]))
proffile.close()

# write surface temperature
timefile = open(outname+'.time','w')
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
    pl.plot(g.dz[:g.kmax], g.z, 'k-o', mfc='none')
    pl.plot(dz_fleur, z_fleur, 'g-x', mfc='none')
    pl.xlabel('dz [m]')
    pl.ylabel('z [m]')

if(False):
    pl.figure()
    pl.subplot(221)
    pl.plot(th, g.z, 'k-', label='mhh')
    pl.plot(s3.th, s3.z, 'go', mfc='none', label='s3')
    pl.ylim(0,1100)
    pl.xlim(270,285)
    pl.legend(frameon=False, loc=2)
    
    pl.subplot(222)
    pl.plot(u, g.z, 'k-', label='mhh')
    pl.plot(s3.u, s3.z, 'go', mfc='none', label='s3')
    pl.plot(ug, g.z, 'k--', label='mhh')
    pl.plot(s3.ug, s3.z, 'bo', mfc='none', label='s3')
    pl.ylim(0,1100)
    pl.xlim(0,10)
    pl.legend(frameon=False, loc=2)
    
    pl.subplot(223)
    pl.plot(v, g.z, 'k-', label='mhh')
    pl.plot(s3.v, s3.z, 'go', mfc='none', label='s3')
    pl.plot(vg, g.z, 'k--', label='mhh')
    pl.plot(s3.vg, s3.z, 'bo', mfc='none', label='s3')
    pl.ylim(0,1100)
    pl.xlim(0,10)
    pl.legend(frameon=False, loc=2)
    
    pl.subplot(224)
    pl.plot(s3.t, s3.ths, 'k-', label='s3')
    pl.legend(frameon=False, loc=2)
