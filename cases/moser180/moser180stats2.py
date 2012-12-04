import numpy
import struct
import netCDF4

from pylab import *

start = 50
end   = 101

# read Moser's data
Mosermean = numpy.loadtxt("chan180.means", skiprows=25)
Moserrey  = numpy.loadtxt("chan180.reystress", skiprows=25)
Moseru2   = numpy.loadtxt("chan180.uubal", skiprows=25)
Moserv2   = numpy.loadtxt("chan180.wwbal", skiprows=25)
Moserw2   = numpy.loadtxt("chan180.vvbal", skiprows=25)
Mosertke  = numpy.loadtxt("chan180.kbal" , skiprows=25)

yplusMoser = Mosermean[:,1]
uavgMoser  = Mosermean[:,2]
uvarMoser  = Moserrey[:,2]
vvarMoser  = Moserrey[:,3]
wvarMoser  = Moserrey[:,4]

u2_shearMoser = Moseru2[:,3]
u2_turbMoser  = Moseru2[:,6]

v2_shearMoser  = Moserv2 [:,3]
tke_shearMoser = Mosertke[:,3]

stats = netCDF4.Dataset("moser180.0000000.nc","r")
t  = stats.variables["t"] [start:end]
z  = stats.variables["z"] [:]
zh = stats.variables["zh"][:]
uavgt = stats.variables["u"] [start:end,:]
vavgt = stats.variables["v"] [start:end,:]
uvart = stats.variables["u2"][start:end,:]
vvart = stats.variables["v2"][start:end,:]
wvart = stats.variables["w2"][start:end,:]

wut    = stats.variables["wu"]   [start:end,:]
wudt   = stats.variables["udiff"][start:end,:]
ufluxt = stats.variables["uflux"][start:end,:]

u2_sheart = stats.variables["u2_shear"][start:end,:]
u2_turbt  = stats.variables["u2_turb"] [start:end,:]

utotavgt = (uavgt**2. + vavgt**2.)**.5
visc   = 1.0e-5
#ustart = (visc * utotavgt[:,0] / z[0])**0.5
ustart = (ufluxt[:,0]**2.)**.25

uavg = numpy.mean(uavgt,0)
vavg = numpy.mean(vavgt,0)
uvar = numpy.mean(uvart,0)
vvar = numpy.mean(vvart,0)
wvar = numpy.mean(wvart,0)

wu    = numpy.mean(wut,0)
wud   = numpy.mean(wudt,0)
uflux = numpy.mean(ufluxt,0)

u2_shear = numpy.mean(u2_sheart,0)
u2_turb  = numpy.mean(u2_turbt ,0)

utotavg = numpy.mean(utotavgt,0)

ustar = numpy.mean(ustart)

print('Re_tau = %.2f' % (ustar / visc))

# create the theoretical lines
ypluslin = arange(0.5,15., 0.1)
ypluslog = arange(5.,800., 1.)
ulin     = ypluslin
ulog     = 2.5 * numpy.log( ypluslog ) + 5.

yplus  = z  * ustar / visc
yplush = zh * ustar / visc

starty = 0
endy   = z.size / 2

close('all')
figure()
for n in range(end-start):
  semilogx(yplus[starty:endy], utotavgt[n,starty:endy] / ustar, color='#cccccc')
semilogx(yplus[starty:endy], utotavg[starty:endy] / ustar, 'bo-', label='u')
semilogx(yplusMoser, uavgMoser, 'k--', label="Moser")
semilogx(ypluslin, ulin, 'k:')
semilogx(ypluslog, ulog, 'k:')
xlabel('y+')
ylabel('u+')
legend(loc=2, frameon=False)
grid()
axis([0.3, 200, 0, 22])

figure()
for n in range(end-start):
  plot(yplus [starty:endy], (uvart[n,starty:endy] / ustar**2.)**0.5, color='#cccccc')
  plot(yplus [starty:endy], (vvart[n,starty:endy] / ustar**2.)**0.5, color='#cccccc')
  plot(yplush[starty:endy], (wvart[n,starty:endy] / ustar**2.)**0.5, color='#cccccc')
plot(yplus [starty:endy], (uvar[starty:endy] / ustar**2.)**0.5, 'b-', label='u')
plot(yplus [starty:endy], (vvar[starty:endy] / ustar**2.)**0.5, 'g-', label='v')
plot(yplush[starty:endy], (wvar[starty:endy] / ustar**2.)**0.5, 'r-', label='w')
plot(yplusMoser, sqrt(uvarMoser), 'k--', label="Moser")
plot(yplusMoser, sqrt(vvarMoser), 'k--')
plot(yplusMoser, sqrt(wvarMoser), 'k--')
xlabel('y+')
ylabel('rms')
legend(loc=0, frameon=False)
grid()
axis([0, 200, 0, 3.5])

figure()
for n in range(end-start):
  plot(yplus[starty:endy], u2_sheart[n,starty:endy] * visc / ustar**4., color='#cccccc')
  plot(yplus[starty:endy], u2_turbt [n,starty:endy] * visc / ustar**4., color='#cccccc')
plot(yplus[starty:endy], u2_shear[starty:endy] * visc / ustar**4., 'b-', label='S')
plot(yplus[starty:endy], u2_turb [starty:endy] * visc / ustar**4., 'g-', label='Tt')
plot(yplusMoser, u2_shearMoser, 'k--', label="Moser")
plot(yplusMoser, u2_turbMoser , 'k--')
xlabel('y+')
ylabel('Rxx')
legend(loc=0, frameon=False)
grid()
axis([0, 200, -0.6, 0.6])

figure()
#for n in range(end-start):
#  plot(ufluxt[n,:] / ustar**2., zh, color='#cccccc')
plot(uflux / ustar**2., zh, 'b-' , label='total flux')
plot(wu    / ustar**2., zh, 'b--', label='turbulent flux')
plot(wud   / ustar**2., zh, 'b:' , label='diffusive flux')
xlabel('uflux')
ylabel('y+')
legend(loc=0, frameon=False)
grid()
axis([-1.1, 1.1, 0., 2.])
 
