import numpy
import struct
import netCDF4

from pylab import *

start = 37
end   = 73

stats = netCDF4.Dataset("bomex.0000000.nc","r")
t   = stats.variables["t"][:]
z   = stats.variables["z"][:]
zh  = stats.variables["zh"][:]

st  = stats.variables["s"][:,:]
qtt = stats.variables["qt"][:,:]*1000.
qlt = stats.variables["ql"][:,:]*1000.
cft = stats.variables["cfrac"][:,:]

sfluxt = stats.variables["sflux"][:,:]

s  = numpy.mean(st [start:end,:],0)
qt = numpy.mean(qtt[start:end,:],0)
ql = numpy.mean(qlt[start:end,:],0)
cf = numpy.mean(cft[start:end,:],0)

sflux = numpy.mean(sfluxt[start:end,:],0)

# enable LaTeX plotting
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

close('all')
figure()
for n in range(start,end):
  plot(st[n,:], z, color='#eeeeee')
plot(s, z)
xlabel(r'$\theta [K]$')
ylabel(r'$z [m]$')

figure()
for n in range(start,end):
  plot(qtt[n,:], z, color='#eeeeee')
plot(qt, z)
xlabel(r'$q_t [g/kg]$')
ylabel(r'$z [m]$')

figure()
for n in range(start,end):
  plot(qlt[n,:], z, color='#eeeeee')
plot(ql, z)
xlabel(r'$q_l [g/kg]$')
ylabel(r'$z [m]$')

figure()
for n in range(start,end):
  plot(cft[n,:], z, color='#eeeeee')
plot(cf, z)
xlabel(r'$q_l [g/kg]$')
ylabel(r'$z [m]$')

figure()
for n in range(start,end):
  plot(sfluxt[n,:], zh, color='#eeeeee')
plot(sflux, zh)
xlabel(r'$w`\theta_l` [K m/s]$')
ylabel(r'$z [m]$')

