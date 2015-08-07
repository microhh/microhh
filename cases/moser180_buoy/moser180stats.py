import numpy
import struct
import netCDF4

from pylab import *

start = 10
end   = 12
plotens = True

stats = netCDF4.Dataset("moser180.default.0000000.nc","r")
t  = stats.variables["t"] [start:end]
z  = stats.variables["z"] [:]
zh = stats.variables["zh"][:]
uavgt = stats.variables["u"] [start:end,:]
vavgt = stats.variables["v"] [start:end,:]
uvart = stats.variables["u2"][start:end,:]
vvart = stats.variables["v2"][start:end,:]
wvart = stats.variables["w2"][start:end,:]
bvart = stats.variables["b2"][start:end,:]

uwt    = stats.variables["uw"]   [start:end,:]
udifft = stats.variables["udiff"][start:end,:]
ufluxt = stats.variables["uflux"][start:end,:]

bwt    = stats.variables["bw"]   [start:end,:]
bdifft = stats.variables["bdiff"][start:end,:]
bfluxt = stats.variables["bflux"][start:end,:]

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

w2_buoyt  = stats.variables["w2_buoy"] [start:end,:]
w2_turbt  = stats.variables["w2_turb"] [start:end,:]
w2_visct  = stats.variables["w2_visc"] [start:end,:]
w2_disst  = stats.variables["w2_diss"] [start:end,:]
w2_prest  = stats.variables["w2_pres"] [start:end,:]
w2_rdstrt = stats.variables["w2_rdstr"][start:end,:]

# tke budget
tke_sheart = stats.variables["tke_shear"][start:end,:]
tke_buoyt  = stats.variables["tke_buoy"] [start:end,:]
tke_turbt  = stats.variables["tke_turb"] [start:end,:]
tke_visct  = stats.variables["tke_visc"] [start:end,:]
tke_disst  = stats.variables["tke_diss"] [start:end,:]
tke_prest  = stats.variables["tke_pres"] [start:end,:]

# stress budgets
uw_sheart = stats.variables["uw_shear"][start:end,:]
uw_visct  = stats.variables["uw_visc"] [start:end,:]
uw_turbt  = stats.variables["uw_turb"] [start:end,:]
uw_disst  = stats.variables["uw_diss"] [start:end,:]
uw_prest  = stats.variables["uw_pres"] [start:end,:]
uw_buoyt  = stats.variables["uw_buoy"] [start:end,:]
uw_rdstrt = stats.variables["uw_rdstr"][start:end,:]

# buoyancy varianc budget
b2_sheart = stats.variables["b2_shear"][start:end,:]
b2_turbt  = stats.variables["b2_turb"] [start:end,:]
b2_visct  = stats.variables["b2_visc"] [start:end,:]
b2_disst  = stats.variables["b2_diss"] [start:end,:]

utotavgt = (uavgt**2. + vavgt**2.)**.5
visc   = 1.0e-5
#ustart = (visc * utotavgt[:,0] / z[0])**0.5
ustart = (ufluxt[:,0]**2.)**.25

uavg = numpy.mean(uavgt,0)
vavg = numpy.mean(vavgt,0)
uvar = numpy.mean(uvart,0)
vvar = numpy.mean(vvart,0)
wvar = numpy.mean(wvart,0)
bvar = numpy.mean(bvart,0)

uw    = numpy.mean(uwt,0)
udiff = numpy.mean(udifft,0)
uflux = numpy.mean(ufluxt,0)

bw    = numpy.mean(bwt,0)
bdiff = numpy.mean(bdifft,0)
bflux = numpy.mean(bfluxt,0)

u2_shear = numpy.mean(u2_sheart,0)
u2_turb  = numpy.mean(u2_turbt ,0)
u2_visc  = numpy.mean(u2_visct ,0)
u2_diss  = numpy.mean(u2_disst ,0)
u2_rdstr = numpy.mean(u2_rdstrt,0)
u2_resid = u2_shear + u2_turb + u2_visc + u2_diss + u2_rdstr

v2_turb  = numpy.mean(v2_turbt ,0)
v2_visc  = numpy.mean(v2_visct ,0)
v2_diss  = numpy.mean(v2_disst ,0)
v2_rdstr = numpy.mean(v2_rdstrt,0)
v2_resid = v2_turb + v2_visc + v2_diss + v2_rdstr

w2_turb  = numpy.mean(w2_turbt ,0)
w2_buoy  = numpy.mean(w2_buoyt ,0)
w2_visc  = numpy.mean(w2_visct ,0)
w2_diss  = numpy.mean(w2_disst ,0)
w2_pres  = numpy.mean(w2_prest ,0)
w2_rdstr = numpy.mean(w2_rdstrt,0)
w2_resid = w2_turb + w2_visc + w2_diss + w2_pres + w2_rdstr + w2_buoy

tke_shear = numpy.mean(tke_sheart,0)
tke_buoy  = numpy.mean(tke_buoyt ,0)
tke_turb  = numpy.mean(tke_turbt ,0)
tke_visc  = numpy.mean(tke_visct ,0)
tke_diss  = numpy.mean(tke_disst ,0)
tke_pres  = numpy.mean(tke_prest ,0)
tke_resid = tke_shear + tke_turb + tke_visc + tke_diss + tke_pres + tke_buoy
 
uw_shear = numpy.mean(uw_sheart,0)
uw_turb  = numpy.mean(uw_turbt ,0)
uw_buoy  = numpy.mean(uw_buoyt ,0)
uw_visc  = numpy.mean(uw_visct ,0)
uw_diss  = numpy.mean(uw_disst ,0)
uw_pres  = numpy.mean(uw_prest ,0)
uw_rdstr = numpy.mean(uw_rdstrt,0)
uw_resid = uw_shear + uw_turb + uw_buoy + uw_visc + uw_diss + uw_pres + uw_rdstr

b2_shear = numpy.mean(b2_sheart,0)
b2_turb  = numpy.mean(b2_turbt ,0)
b2_visc  = numpy.mean(b2_visct ,0)
b2_diss  = numpy.mean(b2_disst ,0)
b2_resid = b2_shear + b2_turb + b2_visc + b2_diss

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
if(plotens):
 for n in range(end-start):
   semilogx(yplus[starty:endy], utotavgt[n,starty:endy] / ustar, color='#cccccc')
semilogx(yplus[starty:endy], utotavg[starty:endy] / ustar, 'b-', label='u')
semilogx(ypluslin, ulin, 'k:')
semilogx(ypluslog, ulog, 'k:')
xlabel('y+')
ylabel('u+')
legend(loc=2, frameon=False)
grid()
xlim(0, 100)

figure()
if(plotens):
  for n in range(end-start):
    plot(yplus [starty:endy], (uvart[n,starty:endy] / ustar**2.)**0.5, color='#cccccc')
    plot(yplus [starty:endy], (vvart[n,starty:endy] / ustar**2.)**0.5, color='#cccccc')
    plot(yplush[starty:endy], (wvart[n,starty:endy] / ustar**2.)**0.5, color='#cccccc')
    plot(yplus [starty:endy], (bvart[n,starty:endy] / (bflux[0]/ustar)**2.)**0.5, color='#cccccc')
plot(yplus [starty:endy], (uvar[starty:endy] / ustar**2.)**0.5, 'b-', label='u')
plot(yplus [starty:endy], (vvar[starty:endy] / ustar**2.)**0.5, 'g-', label='v')
plot(yplush[starty:endy], (wvar[starty:endy] / ustar**2.)**0.5, 'r-', label='w')
plot(yplus [starty:endy], (bvar[starty:endy] / (bflux[0]/ustar)**2.)**.5, 'm-', label='b')
xlabel('y+')
ylabel('rms')
legend(loc=0, frameon=False)
grid()
xlim(0, 100)

figure()
if(plotens):
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
plot(yplus[starty:endy], u2_resid[starty:endy] * visc / ustar**4., 'k--', label='resid')
xlabel('y+')
ylabel('Rxx')
legend(loc=0, ncol=2, frameon=False)
grid()
xlim(0, 100)

figure()
if(plotens):
  for n in range(end-start):
    plot(yplus[starty:endy], v2_turbt [n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplus[starty:endy], v2_visct [n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplus[starty:endy], v2_disst [n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplus[starty:endy], v2_rdstrt[n,starty:endy] * visc / ustar**4., color='#cccccc')
plot(yplus[starty:endy], v2_turb [starty:endy] * visc / ustar**4., 'g-', label='Tt')
plot(yplus[starty:endy], v2_visc [starty:endy] * visc / ustar**4., 'c-', label='Tv')
plot(yplus[starty:endy], v2_diss [starty:endy] * visc / ustar**4., 'r-', label='D')
plot(yplus[starty:endy], v2_rdstr[starty:endy] * visc / ustar**4., 'm-', label='P')
plot(yplus[starty:endy], v2_resid[starty:endy] * visc / ustar**4., 'k--', label='resid')
xlabel('y+')
ylabel('Ryy')
legend(loc=0, ncol=2, frameon=False)
grid()
xlim(0, 100)

figure()
if(plotens):
  for n in range(end-start):
    plot(yplush[starty:endy], w2_turbt [n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplush[starty:endy], w2_visct [n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplush[starty:endy], w2_disst [n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplush[starty:endy], w2_prest [n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplush[starty:endy], w2_rdstrt[n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplush[starty:endy], w2_buoyt [n,starty:endy] * visc / ustar**4., color='#cccccc')
plot(yplush[starty:endy], w2_turb [starty:endy] * visc / ustar**4., 'g-', label='Tt')
plot(yplush[starty:endy], w2_visc [starty:endy] * visc / ustar**4., 'c-', label='Tv')
plot(yplush[starty:endy], w2_diss [starty:endy] * visc / ustar**4., 'r-', label='D')
plot(yplush[starty:endy], w2_pres [starty:endy] * visc / ustar**4., 'y-', label='Tp')
plot(yplush[starty:endy], w2_rdstr[starty:endy] * visc / ustar**4., 'm-', label='P')
plot(yplush[starty:endy], w2_buoy [starty:endy] * visc / ustar**4., 'k-', label='B')
plot(yplush[starty:endy], w2_resid[starty:endy] * visc / ustar**4., 'k--', label='resid')
xlabel('y+')
ylabel('Rzz')
legend(loc=0, ncol=2, frameon=False)
grid()
xlim(0, 100)

figure()
if(plotens):
  for n in range(end-start):
    plot(yplus[starty:endy], tke_sheart[n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplus[starty:endy], tke_turbt [n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplus[starty:endy], tke_buoyt [n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplus[starty:endy], tke_visct [n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplus[starty:endy], tke_disst [n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplus[starty:endy], tke_prest [n,starty:endy] * visc / ustar**4., color='#cccccc')
plot(yplus[starty:endy], tke_shear[starty:endy] * visc / ustar**4., 'b-', label='S')
plot(yplus[starty:endy], tke_turb [starty:endy] * visc / ustar**4., 'g-', label='Tt')
plot(yplus[starty:endy], tke_visc [starty:endy] * visc / ustar**4., 'c-', label='Tv')
plot(yplus[starty:endy], tke_diss [starty:endy] * visc / ustar**4., 'r-', label='D')
plot(yplus[starty:endy], tke_pres [starty:endy] * visc / ustar**4., 'y-', label='Tp')
plot(yplus[starty:endy], tke_buoy [starty:endy] * visc / ustar**4., 'k-', label='B')
plot(yplus[starty:endy], tke_resid[starty:endy] * visc / ustar**4., 'k--', label='resid')
xlabel('y+')
ylabel('tke')
legend(loc=0, ncol=2, frameon=False)
grid()
xlim(0, 100)

figure()
if(plotens):
  for n in range(end-start):
    plot(yplus[starty:endy], uw_sheart[n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplus[starty:endy], uw_turbt [n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplus[starty:endy], uw_visct [n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplus[starty:endy], uw_disst [n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplus[starty:endy], uw_prest [n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplus[starty:endy], uw_rdstrt[n,starty:endy] * visc / ustar**4., color='#cccccc')
    plot(yplus[starty:endy], uw_buoyt [n,starty:endy] * visc / ustar**4., color='#cccccc')
plot(yplus[starty:endy], uw_shear[starty:endy] * visc / ustar**4., 'b-', label='S')
plot(yplus[starty:endy], uw_turb [starty:endy] * visc / ustar**4., 'g-', label='Tt')
plot(yplus[starty:endy], uw_visc [starty:endy] * visc / ustar**4., 'c-', label='Tv')
plot(yplus[starty:endy], uw_diss [starty:endy] * visc / ustar**4., 'r-', label='D')
plot(yplus[starty:endy], uw_pres [starty:endy] * visc / ustar**4., 'y-', label='Tp')
plot(yplus[starty:endy], uw_rdstr[starty:endy] * visc / ustar**4., 'm-', label='P')
plot(yplus[starty:endy], uw_buoy [starty:endy] * visc / ustar**4., 'k-', label='B')
plot(yplus[starty:endy], uw_resid[starty:endy] * visc / ustar**4., 'k--', label='resid')
xlabel('y+')
ylabel('Rxz')
legend(loc=0, ncol=2, frameon=False)
grid()
xlim(0, 100)

figure()
if(plotens):
  for n in range(end-start):
    plot(yplus[starty:endy], b2_sheart[n,starty:endy] / (bflux[0]/ustar)**2, color='#cccccc')
    plot(yplus[starty:endy], b2_turbt [n,starty:endy] / (bflux[0]/ustar)**2, color='#cccccc')
    plot(yplus[starty:endy], b2_visct [n,starty:endy] / (bflux[0]/ustar)**2, color='#cccccc')
    plot(yplus[starty:endy], b2_disst [n,starty:endy] / (bflux[0]/ustar)**2, color='#cccccc')
plot(yplus[starty:endy], b2_shear[starty:endy] / (bflux[0]/ustar)**2, 'b-', label='S')
plot(yplus[starty:endy], b2_turb [starty:endy] / (bflux[0]/ustar)**2, 'g-', label='Tt')
plot(yplus[starty:endy], b2_visc [starty:endy] / (bflux[0]/ustar)**2, 'c-', label='Tv')
plot(yplus[starty:endy], b2_diss [starty:endy] / (bflux[0]/ustar)**2, 'r-', label='D')
plot(yplus[starty:endy], b2_resid[starty:endy] / (bflux[0]/ustar)**2, 'k--', label='resid')
xlabel('y+')
ylabel('B2')
legend(loc=0, ncol=2, frameon=False)
grid()
xlim(0, 100)


figure()
if(plotens):
  for n in range(end-start):
    plot(zh, ufluxt[n,:] / ustar**2., color='#cccccc')
    plot(zh, uwt   [n,:] / ustar**2., color='#cccccc')
    plot(zh, udifft[n,:] / ustar**2., color='#cccccc')
plot(zh, uflux / ustar**2., 'b-' , label='total flux')
plot(zh, uw    / ustar**2., 'b--', label='turbulent flux')
plot(zh, udiff / ustar**2., 'b:' , label='diffusive flux')
xlabel('y+')
ylabel('uflux')
legend(loc=0, frameon=False)
grid()
axis([0., 2., -1.1, 1.1])
 
figure()
if(plotens):
  for n in range(end-start):
    plot(zh, bfluxt[n,:] / bfluxt[n,0], color='#cccccc')
    plot(zh, bwt   [n,:] / bfluxt[n,0], color='#cccccc')
    plot(zh, bdifft[n,:] / bfluxt[n,0], color='#cccccc')
plot(zh, bflux / bflux[0], 'b-' , label='total flux')
plot(zh, bw    / bflux[0], 'b--', label='turbulent flux')
plot(zh, bdiff / bflux[0], 'b:' , label='diffusive flux')
xlabel('y+')
ylabel('bflux')
legend(loc=0, frameon=False)
grid()
axis([0., 2., -0.1, 1.1])
