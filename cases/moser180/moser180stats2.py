import numpy
import struct
import netCDF4

from pylab import *

start = 20
end   = 30

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
u2_viscMoser  = Moseru2[:,7]
u2_dissMoser  = Moseru2[:,2]
u2_rdstrMoser = Moseru2[:,4]

v2_turbMoser  = Moserv2[:,6]
v2_viscMoser  = Moserv2[:,7]
v2_dissMoser  = Moserv2[:,2]
v2_rdstrMoser = Moserv2[:,4]

w2_turbMoser  = Moserw2[:,6]
w2_viscMoser  = Moserw2[:,7]
w2_dissMoser  = Moserw2[:,2]
w2_presMoser  = Moserw2[:,5]
w2_rdstrMoser = Moserw2[:,4]

tke_shearMoser = Mosertke[:,3]
tke_turbMoser  = Mosertke[:,6]
tke_viscMoser  = Mosertke[:,7]
tke_dissMoser  = Mosertke[:,2]
tke_presMoser  = Mosertke[:,5]

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
udifft = stats.variables["udiff"][start:end,:]
ufluxt = stats.variables["uflux"][start:end,:]

# variance budgets
u2_sheart = stats.variables["u2_shear"][start:end,:]
u2_turbt  = stats.variables["u2_turb"] [start:end,:]
u2_visct  = stats.variables["u2_visc"] [start:end,:]
u2_disst  = stats.variables["u2_diss"] [start:end,:]
u2_rdstrt = stats.variables["u2_rdstr"][start:end,:]

v2_turbt  = stats.variables["v2_turb"] [start:end,:]
v2_visct  = stats.variables["v2_visc"] [start:end,:]
v2_disst  = stats.variables["v2_diss"] [start:end,:]
v2_rdstrt = stats.variables["v2_rdstr"][start:end,:]

w2_turbt  = stats.variables["w2_turb"] [start:end,:]
w2_visct  = stats.variables["w2_visc"] [start:end,:]
w2_disst  = stats.variables["w2_diss"] [start:end,:]
w2_prest  = stats.variables["w2_pres"] [start:end,:]
w2_rdstrt = stats.variables["w2_rdstr"][start:end,:]

# tke budget
tke_sheart = stats.variables["tke_shear"][start:end,:]
tke_turbt  = stats.variables["tke_turb"] [start:end,:]
tke_visct  = stats.variables["tke_visc"] [start:end,:]
tke_disst  = stats.variables["tke_diss"] [start:end,:]
tke_prest  = stats.variables["tke_pres"] [start:end,:]

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
udiff = numpy.mean(udifft,0)
uflux = numpy.mean(ufluxt,0)

u2_shear = numpy.mean(u2_sheart,0)
u2_turb  = numpy.mean(u2_turbt ,0)
u2_visc  = numpy.mean(u2_visct ,0)
u2_diss  = numpy.mean(u2_disst ,0)
u2_rdstr = numpy.mean(u2_rdstrt,0)

v2_turb  = numpy.mean(v2_turbt ,0)
v2_visc  = numpy.mean(v2_visct ,0)
v2_diss  = numpy.mean(v2_disst ,0)
v2_rdstr = numpy.mean(v2_rdstrt,0)

w2_turb  = numpy.mean(w2_turbt ,0)
w2_visc  = numpy.mean(w2_visct ,0)
w2_diss  = numpy.mean(w2_disst ,0)
w2_pres  = numpy.mean(w2_prest ,0)
w2_rdstr = numpy.mean(w2_rdstrt,0)

tke_shear = numpy.mean(tke_sheart,0)
tke_turb  = numpy.mean(tke_turbt ,0)
tke_visc  = numpy.mean(tke_visct ,0)
tke_diss  = numpy.mean(tke_disst ,0)
tke_pres  = numpy.mean(tke_prest ,0)

utotavg = numpy.mean(utotavgt,0)
ustar   = numpy.mean(ustart)

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
semilogx(yplus[starty:endy], utotavg[starty:endy] / ustar, 'b-', label='u')
semilogx(yplusMoser, uavgMoser, 'k--', label="Moser")
semilogx(ypluslin, ulin, 'k:')
semilogx(ypluslog, ulog, 'k:')
xlabel('y+')
ylabel('u+')
legend(loc=2, frameon=False)
grid()
axis([0.1, 200, 0, 20])

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
  plot(yplus[starty:endy], u2_visct [n,starty:endy] * visc / ustar**4., color='#cccccc')
  plot(yplus[starty:endy], u2_disst [n,starty:endy] * visc / ustar**4., color='#cccccc')
  plot(yplus[starty:endy], u2_rdstrt[n,starty:endy] * visc / ustar**4., color='#cccccc')
plot(yplus[starty:endy], u2_shear[starty:endy] * visc / ustar**4., 'b-', label='S')
plot(yplus[starty:endy], u2_turb [starty:endy] * visc / ustar**4., 'g-', label='Tt')
plot(yplus[starty:endy], u2_visc [starty:endy] * visc / ustar**4., 'c-', label='Tv')
plot(yplus[starty:endy], u2_diss [starty:endy] * visc / ustar**4., 'r-', label='D')
plot(yplus[starty:endy], u2_rdstr[starty:endy] * visc / ustar**4., 'm-', label='P')
plot(yplusMoser, u2_shearMoser, 'k--', label="Moser")
plot(yplusMoser, u2_turbMoser , 'k--')
plot(yplusMoser, u2_viscMoser , 'k--')
plot(yplusMoser, u2_dissMoser , 'k--')
plot(yplusMoser, u2_rdstrMoser, 'k--')
xlabel('y+')
ylabel('Rxx')
legend(loc=0, frameon=False)
grid()
axis([0, 200, -0.6, 0.6])

figure()
for n in range(end-start):
  plot(yplus[starty:endy], v2_turbt [n,starty:endy] * visc / ustar**4., color='#cccccc')
  plot(yplus[starty:endy], v2_visct [n,starty:endy] * visc / ustar**4., color='#cccccc')
  plot(yplus[starty:endy], v2_disst [n,starty:endy] * visc / ustar**4., color='#cccccc')
  plot(yplus[starty:endy], v2_rdstrt[n,starty:endy] * visc / ustar**4., color='#cccccc')
plot(yplus[starty:endy], v2_turb [starty:endy] * visc / ustar**4., 'g-', label='Tt')
plot(yplus[starty:endy], v2_visc [starty:endy] * visc / ustar**4., 'c-', label='Tv')
plot(yplus[starty:endy], v2_diss [starty:endy] * visc / ustar**4., 'r-', label='D')
plot(yplus[starty:endy], v2_rdstr[starty:endy] * visc / ustar**4., 'm-', label='P')
plot(yplusMoser, v2_turbMoser , 'k--', label="Moser")
plot(yplusMoser, v2_viscMoser , 'k--')
plot(yplusMoser, v2_dissMoser , 'k--')
plot(yplusMoser, v2_rdstrMoser, 'k--')
xlabel('y+')
ylabel('Ryy')
legend(loc=0, frameon=False)
grid()
axis([0, 200, -0.06, 0.06])

figure()
for n in range(end-start):
  plot(yplush[starty:endy], w2_turbt [n,starty:endy] * visc / ustar**4., color='#cccccc')
  plot(yplush[starty:endy], w2_visct [n,starty:endy] * visc / ustar**4., color='#cccccc')
  plot(yplush[starty:endy], w2_disst [n,starty:endy] * visc / ustar**4., color='#cccccc')
  plot(yplush[starty:endy], w2_prest [n,starty:endy] * visc / ustar**4., color='#cccccc')
  plot(yplush[starty:endy], w2_rdstrt[n,starty:endy] * visc / ustar**4., color='#cccccc')
plot(yplush[starty:endy], w2_turb [starty:endy] * visc / ustar**4., 'g-', label='Tt')
plot(yplush[starty:endy], w2_visc [starty:endy] * visc / ustar**4., 'c-', label='Tv')
plot(yplush[starty:endy], w2_diss [starty:endy] * visc / ustar**4., 'r-', label='D')
plot(yplush[starty:endy], w2_pres [starty:endy] * visc / ustar**4., 'y-', label='Tp')
plot(yplush[starty:endy], w2_rdstr[starty:endy] * visc / ustar**4., 'm-', label='P')
plot(yplusMoser, w2_turbMoser , 'k--', label="Moser")
plot(yplusMoser, w2_viscMoser , 'k--')
plot(yplusMoser, w2_dissMoser , 'k--')
plot(yplusMoser, w2_presMoser , 'k--')
plot(yplusMoser, w2_rdstrMoser, 'k--')
xlabel('y+')
ylabel('Rzz')
legend(loc=0, frameon=False)
grid()
axis([0, 200, -0.06, 0.06])

figure()
for n in range(end-start):
  plot(yplus[starty:endy], tke_sheart[n,starty:endy] * visc / ustar**4., color='#cccccc')
  plot(yplus[starty:endy], tke_turbt [n,starty:endy] * visc / ustar**4., color='#cccccc')
  plot(yplus[starty:endy], tke_visct [n,starty:endy] * visc / ustar**4., color='#cccccc')
  plot(yplus[starty:endy], tke_disst [n,starty:endy] * visc / ustar**4., color='#cccccc')
  plot(yplus[starty:endy], tke_prest [n,starty:endy] * visc / ustar**4., color='#cccccc')
plot(yplus[starty:endy], tke_shear[starty:endy] * visc / ustar**4., 'b-', label='S')
plot(yplus[starty:endy], tke_turb [starty:endy] * visc / ustar**4., 'g-', label='Tt')
plot(yplus[starty:endy], tke_visc [starty:endy] * visc / ustar**4., 'c-', label='Tv')
plot(yplus[starty:endy], tke_diss [starty:endy] * visc / ustar**4., 'r-', label='D')
plot(yplus[starty:endy], tke_pres [starty:endy] * visc / ustar**4., 'y-', label='Tp')
plot(yplusMoser, tke_shearMoser, 'k--', label="Moser")
plot(yplusMoser, tke_turbMoser , 'k--')
plot(yplusMoser, tke_viscMoser , 'k--')
plot(yplusMoser, tke_dissMoser , 'k--')
plot(yplusMoser, tke_presMoser , 'k--')
xlabel('y+')
ylabel('tke')
legend(loc=0, frameon=False)
grid()
#axis([0, 200, -0.6, 0.6])

figure()
for n in range(end-start):
  plot(zh, ufluxt[n,:] / ustar**2., color='#cccccc')
  plot(zh, wut   [n,:] / ustar**2., color='#cccccc')
  plot(zh, udifft[n,:] / ustar**2., color='#cccccc')
plot(zh, uflux / ustar**2., 'b-' , label='total flux')
plot(zh, wu    / ustar**2., 'b--', label='turbulent flux')
plot(zh, udiff / ustar**2., 'b:' , label='diffusive flux')
xlabel('y+')
ylabel('uflux')
legend(loc=0, frameon=False)
grid()
axis([0., 2., -1.1, 1.1])
 
