import numpy
import struct
import netCDF4

from pylab import *

stats = netCDF4.Dataset("bomex.0000000.nc","r")
t   = stats.variables["t"][:]
z   = stats.variables["z"][:]
zh  = stats.variables["zh"][:]

st  = stats.variables["s"][:,:]
qtt = stats.variables["qt"][:,:]*1000.
ut  = stats.variables["u"][:,:]
vt  = stats.variables["v"][:,:]
qlt = stats.variables["ql"][:,:]*1000.
cft = stats.variables["cfrac"][:,:]

sfluxt = stats.variables["sflux"][:,:]
ufluxt = stats.variables["uflux"][:,:]

end   = t.size
start = t.size - 36

s  = numpy.mean(st [start:end,:], 0)
qt = numpy.mean(qtt[start:end,:], 0)
u  = numpy.mean(ut [start:end,:], 0)
v  = numpy.mean(vt [start:end,:], 0)
ql = numpy.mean(qlt[start:end,:], 0)
cf = numpy.mean(cft[start:end,:], 0)

sflux = numpy.mean(sfluxt[start:end,:], 0)
uflux = numpy.mean(ufluxt[start:end,:], 0)

# enable LaTeX plotting
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

close('all')
figure()
for n in range(start,end):
  plot(st[n,:], z, color='#eeeeee')
plot(s, z)
plot(st[0,:], z, 'k-')
xlabel(r'$\theta$ [K]')
ylabel(r'z [m]')

figure()
for n in range(start,end):
  plot(qtt[n,:], z, color='#eeeeee')
plot(qt, z)
plot(qtt[0,:], z, 'k-')
xlabel(r'q$_t$ [g~kg$^{-1}$]')
ylabel(r'z [m]')

figure()
for n in range(start,end):
  plot(ut[n,:], z, color='#eeeeee')
plot(u, z)
plot(ut[0,:], z, 'k-')
xlabel(r'u [m~s$^{-1}$]')
ylabel(r'z [m]')

figure()
for n in range(start,end):
  plot(vt[n,:], z, color='#eeeeee')
plot(v, z)
plot(vt[0,:], z, 'k-')
xlabel(r'v [m~s$^{-1}$]')
ylabel(r'z [m]')

figure()
for n in range(start,end):
  plot(qlt[n,:], z, color='#eeeeee')
plot(ql, z)
xlabel(r'q$_l$ [g~kg$^{-1}$]')
ylabel(r'z [m]')

figure()
for n in range(start,end):
  plot(cft[n,:], z, color='#eeeeee')
plot(cf, z)
xlabel(r'cloud fraction [-]')
ylabel(r'z [m]')

figure()
for n in range(start,end):
  plot(sfluxt[n,:], zh, color='#eeeeee')
plot(sflux, zh)
xlabel(r'w`$\theta_l$` [K~m~s$^{-1}$]')
ylabel(r'z [m]')

figure()
for n in range(start,end):
  plot(ufluxt[n,:], zh, color='#eeeeee')
plot(uflux, zh)
xlabel(r'u`w` [m$^2$~s$^{-2}$]')
ylabel(r'z [m]')

