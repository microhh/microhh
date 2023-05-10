from netCDF4 import Dataset
from pylab import *
import numpy as np

slist = 'hcho,isopooh,isop,mvkmacr,no,o3'
specs = slist.split(',')
profs = []
flux = []
error = []
with Dataset("test.nc") as x:
    zax = x.variables['z'][:]
    zaxh = x.variables['zh'][:]
    z = x.groups['default']
    for spec in specs:
        profs.append(z.variables[spec][:])
        spece = spec+'_2'
        error.append(z.variables[spece][:])
        spece = spec+'_flux'
        flux.append(z.variables[spece][:])
    z = x.groups['thermo']
    thl = z.variables['thl'][:]
    thl_2 = z.variables['thl_2'][:]

profs = np.array(profs)
error = np.array(error)
flux = np.array(flux)

for i,spec in enumerate(specs):
    f,ax = subplots()
    for j in range(50):
        ax.plot(flux[i,j,:], zaxh)
    ax.set_xlabel(spec + '(m/s)')
    ax.set_ylabel('z (m)')
for i,spec in enumerate(specs):
    f,ax = subplots()
    for j in range(50):
        ax.plot(profs[i,j,:], zax)
    ax.errorbar(profs[i,-1,:],zax, xerr=sqrt(error[i,-1,:]), fmt='.')
    ax.set_xlabel(spec + '(ppb)')
    ax.set_ylabel('z (m)')
    f.savefig(spec+'ed_32_4e2.png')
