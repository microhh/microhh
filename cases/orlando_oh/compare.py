from netCDF4 import Dataset
from pylab import *
import numpy as np
import sys

slist = 'ch4,h2o2,n2o5,hald,co,hcho,isopooh,isop,xo2,mvkmacr,isopao2,no2,no3,ch3o2,isopbo2,no,ho2,o3,oh'
#slist = 'hcho,isopooh,isop,mvkmacr,no,o3,oh'
specs = slist.split(',')
profs = []
error = []
flux = []
with Dataset("orlando.default.0010800.nc") as x:
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

specs = slist.split(',')
profs_1 = []
error_1 = []
flux_1 = []
with Dataset("switch_dt_1e2/orlando.default.0010800.nc") as x:
    zax = x.variables['z'][:]
    zaxh = x.variables['zh'][:]
    z = x.groups['default']
    for spec in specs:
        profs_1.append(z.variables[spec][:])
        spece = spec+'_2'
        error_1.append(z.variables[spece][:])
        if spec != specs[-1]:
            spece = spec+'_flux'
            flux_1.append(z.variables[spece][:])
    z = x.groups['thermo']

profs = array(profs)
profs_1 = array(profs_1)
for i,spec in enumerate(specs):
    f,ax = subplots()
    ax.plot(profs[i,-1,:], zax)
    ax.plot(profs_1[i,-1,:], zax)
    ax.set_xlabel(spec + '(ppb)')
    ax.set_ylabel('z (m)')
    ax.set_ylim((0,4000))
    f.savefig(spec+'ed_64_1e2.png')
