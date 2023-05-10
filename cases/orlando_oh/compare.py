from netCDF4 import Dataset
from pylab import *
import numpy as np
import sys

xtime = '0021600'
xtime = '0021600'
projs = ['dt_100','./','dt_1e4']
slist = 'ch4,h2o2,n2o5,hald,co,hcho,isopooh,isop,xo2,mvkmacr,isopao2,no2,no3,ch3o2,isopbo2,no,ho2,o3,oh'
#slist = 'hcho,isopooh,isop,mvkmacr,no,o3,oh'
specs = slist.split(',')
aprofs = []
th = []
for iproj,proj in enumerate(projs):
    with Dataset(proj+'/orlando.default.'+xtime+'.nc') as x:
        profs = []
        error = []
        flux = []
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
        aprofs.append(array(profs))
        th.append(array(thl))

f,ax = subplots()
for ip, xprofs in enumerate(th):
    ax.plot(xprofs[-1,:], zax, label = projs[ip])
    ax.set_xlabel('thl (K)')
    ax.legend(loc='best')
for i,spec in enumerate(specs):
    f,ax = subplots()
    for ip, xprofs in enumerate(aprofs):
        ax.plot(xprofs[i,-1,:], zax, label = projs[ip])
    ax.set_xlabel(spec + '(ppb)')
    ax.set_ylabel('z (m)')
    ax.set_ylim((0,4000))
    ax.legend(loc='best')
