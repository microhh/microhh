import numpy
import struct
import netCDF4

from pylab import *

stats0 = netCDF4.Dataset("bomex.default.0000000.nc","r")
t0   = stats0.variables["t"][:]
z0   = stats0.variables["z"][:]
qtt0 = stats0.variables["qt"][:,:]
evisct0 = stats0.variables["evisc"][:,:]
area0 = stats0.variables["area"][:,:]

stats1 = netCDF4.Dataset("bomex.wplus.0000000.nc","r")
t1   = stats1.variables["t"][:]
z1   = stats1.variables["z"][:]
qtt1 = stats1.variables["qt"][:,:]
evisct1 = stats1.variables["evisc"][:,:]
area1 = stats1.variables["area"][:,:]

stats2 = netCDF4.Dataset("bomex.wmin.0000000.nc","r")
t2   = stats2.variables["t"][:]
z2   = stats2.variables["z"][:]
qtt2 = stats2.variables["qt"][:,:]
evisct2 = stats2.variables["evisc"][:,:]
area2 = stats2.variables["area"][:,:]

print(area1+area2)
print(qtt0 - (qtt1*area1+qtt2*area2))
# enable LaTeX plotting
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

close('all')
figure()
plot(qtt0[-1,:], z0)
plot(qtt1[-1,:], z1)
plot(qtt2[-1,:], z2)
plot(area1[-1,:]*qtt1[-1,:] + area2[-1,:]*qtt2[-1,:], z0)

figure()
plot(area0[-1,:], z0)
plot(area1[-1,:], z1)
plot(area2[-1,:], z2)
plot(area1[-1,:]+area2[-1,:], z0)

figure()
plot(evisct0[-1,:], z0)
plot(evisct1[-1,:], z1)
plot(evisct2[-1,:], z2)
plot(area1[-1,:]*evisct1[-1,:] + area2[-1,:]*evisct2[-1,:], z0)

