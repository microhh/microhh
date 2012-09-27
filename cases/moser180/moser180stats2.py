import numpy
import struct
import netCDF4

from pylab import *

start = 2
end   = 4

# read Moser's data
Mosermean = numpy.loadtxt("chan180.means", skiprows=25)
Moserrey  = numpy.loadtxt("chan180.reystress", skiprows=25)

yplusMoser = Mosermean[:,1]
uavgMoser  = Mosermean[:,2]
uvarMoser  = Moserrey[:,2]
vvarMoser  = Moserrey[:,3]
wvarMoser  = Moserrey[:,4]

stats = netCDF4.Dataset("test.nc","r")
z  = stats.variables["z"] [:]
zh = stats.variables["zh"][:]
uavgt = stats.variables["u"] [start:end,:]
vavgt = stats.variables["v"] [start:end,:]
uvart = stats.variables["u2"][start:end,:]
vvart = stats.variables["v2"][start:end,:]
wvart = stats.variables["w2"][start:end,:]

utotavgt = (uavgt**2. + vavgt**2.)**.5
visc     = 1.0e-5
ustart   = (visc * utotavgt[:,0] / z[0])**0.5

uavg = numpy.mean(uavgt,0)
vavg = numpy.mean(vavgt,0)
uvar = numpy.mean(uvart,0)
vvar = numpy.mean(vvart,0)
wvar = numpy.mean(wvart,0)

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
for t in range(end-start):
  semilogx(yplus[starty:endy], utotavgt[t,starty:endy] / ustar, color='#cccccc')
semilogx(yplus[starty:endy], utotavg[starty:endy] / ustar, 'b-')
semilogx(yplusMoser, uavgMoser, 'k--', label="Moser")
semilogx(ypluslin, ulin, 'k:')
semilogx(ypluslog, ulog, 'k:')
xlabel('y+')
ylabel('u+')
legend(loc=2, frameon=False)
grid()
axis([0.3, 200, 0, 22])

figure()
for t in range(end-start):
  plot(yplus [starty:endy], (uvart[t,starty:endy] / ustar**2.)**0.5, color='#cccccc')
  plot(yplus [starty:endy], (vvart[t,starty:endy] / ustar**2.)**0.5, color='#cccccc')
  plot(yplush[starty:endy], (wvart[t,starty:endy] / ustar**2.)**0.5, color='#cccccc')
plot(yplus [starty:endy], (uvar[starty:endy] / ustar**2.)**0.5, 'b-')
plot(yplus [starty:endy], (vvar[starty:endy] / ustar**2.)**0.5, 'g-')
plot(yplush[starty:endy], (wvar[starty:endy] / ustar**2.)**0.5, 'r-')
plot(yplusMoser, sqrt(uvarMoser), 'k--', label="Moser")
plot(yplusMoser, sqrt(vvarMoser), 'k--')
plot(yplusMoser, sqrt(wvarMoser), 'k--')
xlabel('y+')
ylabel('rms')
legend(loc=0, frameon=False)
grid()
axis([0, 200, 0, 3.5])
 
