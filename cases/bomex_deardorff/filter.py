import numpy
import struct
import netCDF4

from pylab import *

stats0 = netCDF4.Dataset("bomex.default.0000000.nc","r")
t0   = stats0.variables["t"][:]
z0   = stats0.variables["z"][:]
zh0  = stats0.variables["zh"][:]
vt0  = stats0.variables["v"][:,:]
v2t0 = stats0.variables["v2"][:,:]
bfluxt0 = stats0.variables["bflux"][:,:]
area0 = stats0.variables["area"][:,:]
areah0 = stats0.variables["areah"][:,:]

stats1 = netCDF4.Dataset("bomex.wplus.0000000.nc","r")
t1   = stats1.variables["t"][:]
z1   = stats1.variables["z"][:]
zh1  = stats1.variables["zh"][:]
vt1 = stats1.variables["v"][:,:]
v2t1 = stats1.variables["v2"][:,:]
bfluxt1 = stats1.variables["bflux"][:,:]
area1 = stats1.variables["area"][:,:]
areah1 = stats1.variables["areah"][:,:]

stats2 = netCDF4.Dataset("bomex.wmin.0000000.nc","r")
t2   = stats2.variables["t"][:]
z2   = stats2.variables["z"][:]
zh2  = stats2.variables["zh"][:]
vt2 = stats2.variables["v"][:,:]
v2t2 = stats2.variables["v2"][:,:]
bfluxt2 = stats2.variables["bflux"][:,:]
area2 = stats2.variables["area"][:,:]
areah2 = stats2.variables["areah"][:,:]

#print(areah1[-1,:]+areah2[-1,:])
#print((area1[-1,:]*area2[-1,:]*(vt1[-1,:]-vt2[-1,:])**2. + v2t1[-1,:]*area1[-1,:]+v2t2[-1,:]*area2[-1,:])/(v2t0[-1,:]*area0[-1,:]))
print((bfluxt1[-1,:]*areah1[-1,:]+bfluxt2[-1,:]*areah2[-1,:])/(bfluxt0[-1,:]*areah0[-1,:]))
# enable LaTeX plotting
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

close('all')
#figure()
#plot(vt0[-1,:], z0)
#plot(vt1[-1,:], z1)
#plot(vt2[-1,:], z2)
#plot(area1[-1,:]*vt1[-1,:] + area2[-1,:]*vt2[-1,:], z0)

figure()
plot(v2t0[-1,:], z0)
plot(v2t1[-1,:], z1)
plot(v2t2[-1,:], z2)
plot(area1[-1,:]*v2t1[-1,:] + area2[-1,:]*v2t2[-1,:], z0)

figure()
plot(bfluxt0[-1,:], zh0)
plot(bfluxt1[-1,:], zh1)
plot(bfluxt2[-1,:], zh2)
plot(areah1[-1,:]*bfluxt1[-1,:] + areah2[-1,:]*bfluxt2[-1,:], zh0)

#figure()
#plot(areah0[-1,:], zh0)
#plot(areah1[-1,:], zh1)
#plot(areah2[-1,:], zh2)
#plot(areah1[-1,:]+areah2[-1,:], zh0)

