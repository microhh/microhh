from netCDF4 import Dataset
from pylab import *
import numpy as np
import sys

slist = 'hno3,co,hcho,rooh,h2o2,rh,no2,no,o3,ro2,ho2,oh,n2o5,no3,c2o3,pan,co2'
#slist = 'hcho,isopooh,isop,mvkmacr,no,o3,oh'
specs = slist.split(',')
profs = []
error = []
flux = []
with Dataset("tm5_ifs.default.0032400.nc") as x:
    zax = x.variables['z'][:]
    zaxh = x.variables['zh'][:]
    z = x.groups['default']
    for spec in specs:
        profs.append(z.variables[spec][:])
        spece = spec+'_2'
        error.append(z.variables[spece][:])
        if spec != specs[-1]:
            spece = spec+'_flux'
            flux.append(z.variables[spece][:])
    z = x.groups['thermo']
    thl = z.variables['thl'][:]
    thl_2 = z.variables['thl_2'][:]

f,ax = subplots()
for th in thl:
    ax.plot(th[:], zax)
ax.plot(thl[-1,:],zax, 'ro')
ax.set_ylabel('z (m)')
f.savefig('theta_ed_48_1e1.png')
profs = np.array(profs)

error = np.array(error)
flux = np.array(flux)

#for i,spec in enumerate(specs[:-1]):
#    f,ax = subplots()
#    for j in range(132):
#        ax.plot(flux[i,j,:], zaxh)
#    ax.set_xlabel(spec + '(m/s)')
#    ax.set_ylabel('z (m)')
np = profs.shape[1]
for i,spec in enumerate(specs):
    f,ax = subplots()
    for j in range(0,np,1):
        ax.plot(profs[i,j,:], zax)
    ax.errorbar(profs[i,-1,:],zax, xerr=sqrt(error[i,-1,:]), fmt='.')
    ax.set_xlabel(spec + '(ppb)')
    ax.set_ylabel('z (m)')
    ax.set_ylim((0,3000))
    f.savefig(spec+'ed_small_ed100.png')
