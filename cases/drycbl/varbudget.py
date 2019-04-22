import numpy
import struct
import netCDF4

from pylab import *

stats = netCDF4.Dataset("drycbl_default_0000000.nc","r")

t = stats.variables["time"][:]
end   = t.size
start = t.size-5
dt = t[1] - t[0]

z  = stats.variables["z"][:]
zh = stats.variables["zh"][:]
dz = zh[1::] - zh[0:-1]

B0 = 0.0032
N2 = 3.
L0 = (B0/N2**1.5)**.5
henct = (2.*B0/N2 * t)**.5
benct = henct*N2
wenct = (B0*henct)**(1./3.)
tenct = henct/(B0/N2**1.5)**.5

henc  = numpy.mean(henct[start:end])
benc  = numpy.mean(benct[start:end])
wenc  = numpy.mean(wenct[start:end])

u2_turb  = average(stats.variables["u2_turb" ][start:end,:], 0)
u2_visc  = average(stats.variables["u2_visc" ][start:end,:], 0)
u2_rdstr = average(stats.variables["u2_rdstr"][start:end,:], 0)
u2_diss  = average(stats.variables["u2_diss" ][start:end,:], 0)
u2_sum = u2_turb + u2_visc + u2_rdstr + u2_diss

v2_turb  = average(stats.variables["v2_turb" ][start:end,:], 0)
v2_visc  = average(stats.variables["v2_visc" ][start:end,:], 0)
v2_rdstr = average(stats.variables["v2_rdstr"][start:end,:], 0)
v2_diss  = average(stats.variables["v2_diss" ][start:end,:], 0)
v2_sum = v2_turb + v2_visc + v2_rdstr + v2_diss

w2_turb  = average(stats.variables["w2_turb" ][start:end,:], 0)
w2_visc  = average(stats.variables["w2_visc" ][start:end,:], 0)
w2_pres  = average(stats.variables["w2_pres" ][start:end,:], 0)
w2_rdstr = average(stats.variables["w2_rdstr"][start:end,:], 0)
w2_diss  = average(stats.variables["w2_diss" ][start:end,:], 0)
w2_buoy  = average(stats.variables["w2_buoy" ][start:end,:], 0)
w2_sum = w2_turb + w2_visc + w2_pres + w2_rdstr + w2_diss + w2_buoy

tke_turb = average(stats.variables["tke_turb"][start:end,:], 0)
tke_visc = average(stats.variables["tke_visc"][start:end,:], 0)
tke_pres = average(stats.variables["tke_pres"][start:end,:], 0)
tke_diss = average(stats.variables["tke_diss"][start:end,:], 0)
tke_buoy = average(stats.variables["tke_buoy"][start:end,:], 0)
tke_sum = tke_turb + tke_visc + tke_pres + tke_diss + tke_buoy

b2_shear = average(stats.variables["b2_shear"][start:end,:], 0)
b2_turb  = average(stats.variables["b2_turb" ][start:end,:], 0)
b2_visc  = average(stats.variables["b2_visc" ][start:end,:], 0)
b2_diss  = average(stats.variables["b2_diss" ][start:end,:], 0)
b2_sum = b2_shear + b2_turb + b2_visc + b2_diss

bw_shear = average(stats.variables["bw_shear"][start:end,:], 0)
bw_turb  = average(stats.variables["bw_turb" ][start:end,:], 0)
bw_visc  = average(stats.variables["bw_visc" ][start:end,:], 0)
bw_buoy  = average(stats.variables["bw_buoy" ][start:end,:], 0)
bw_rdstr = average(stats.variables["bw_rdstr"][start:end,:], 0)
bw_diss  = average(stats.variables["bw_diss" ][start:end,:], 0)
bw_pres  = average(stats.variables["bw_pres" ][start:end,:], 0)
bw_sum = bw_shear + bw_turb + bw_visc + bw_buoy + bw_rdstr + bw_diss + bw_pres


# enable LaTeX plotting
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

zlim = 1.8

close('all')
figure()
plot(u2_turb /B0, z/henc, 'b-', label='turb' )
plot(u2_visc /B0, z/henc, 'g-', label='visc' )
plot(u2_diss /B0, z/henc, 'r-', label='diss' )
plot(u2_rdstr/B0, z/henc, 'y-', label='rdstr')
plot(u2_sum  /B0, z/henc, 'k:', label='sum'  )
xlabel(r'$du^2/dt / B_0$')
ylabel(r'$z/h_{enc}$')
legend(loc=0, frameon=False)
ylim(0,zlim)

figure()
plot(v2_turb /B0, z/henc, 'b-', label='turb' )
plot(v2_visc /B0, z/henc, 'g-', label='visc' )
plot(v2_diss /B0, z/henc, 'r-', label='diss' )
plot(v2_rdstr/B0, z/henc, 'y-', label='rdstr')
plot(v2_sum  /B0, z/henc, 'k:', label='sum'  )
xlabel(r'$dv^2/dt / B_0$')
ylabel(r'$z/h_{enc}$')
legend(loc=0, frameon=False)
ylim(0,zlim)

figure()
plot(w2_turb /B0, zh/henc, 'b-', label='turb' )
plot(w2_visc /B0, zh/henc, 'g-', label='visc' )
plot(w2_diss /B0, zh/henc, 'r-', label='diss' )
plot(w2_pres /B0, zh/henc, 'm-', label='pres' )
plot(w2_rdstr/B0, zh/henc, 'y-', label='rdstr')
plot(w2_buoy /B0, zh/henc, 'c-', label='buoy' )
plot(w2_sum  /B0, zh/henc, 'k:', label='sum'  )
xlabel(r'$dw^2/dt / B_0$')
ylabel(r'$z/h_{enc}$')
legend(loc=0, frameon=False)
ylim(0,zlim)

figure()
plot(tke_turb/B0, z/henc, 'b-', label='turb')
plot(tke_visc/B0, z/henc, 'g-', label='visc')
plot(tke_diss/B0, z/henc, 'r-', label='diss')
plot(tke_pres/B0, z/henc, 'm-', label='pres')
plot(tke_buoy/B0, z/henc, 'c-', label='buoy')
plot(tke_sum /B0, z/henc, 'k:', label='sum' )
xlabel(r'$de_k/dt / B_0$')
ylabel(r'$z/h_{enc}$')
legend(loc=0, frameon=False)
ylim(0,zlim)

figure()
plot(b2_shear/B0, z/henc, 'c-', label='shear')
plot(b2_turb /B0, z/henc, 'b-', label='turb' )
plot(b2_visc /B0, z/henc, 'g-', label='visc' )
plot(b2_diss /B0, z/henc, 'r-', label='diss' )
plot(b2_sum  /B0, z/henc, 'k:', label='sum'  )
xlabel(r'$db^2/dt / B_0$')
ylabel(r'$z/h_{enc}$')
legend(loc=0, frameon=False)
ylim(0,zlim)

figure()
plot(bw_shear/B0, zh/henc, 'c-', label='shear' )
plot(bw_turb /B0, zh/henc, 'b-', label='turb' )
plot(bw_visc /B0, zh/henc, 'g-', label='visc' )
plot(bw_buoy /B0, zh/henc, 'k-', label='buoy' )
plot(bw_rdstr/B0, zh/henc, 'y-', label='rdstr')
plot(bw_diss /B0, zh/henc, 'r-', label='diss' )
plot(bw_pres /B0, zh/henc, 'm-', label='pres' )
plot(bw_sum  /B0, zh/henc, 'k:', label='sum'  )
xlabel(r'$dw^\prime b^\prime /dt / B_0$')
ylabel(r'$z/h_{enc}$')
legend(loc=0, frameon=False)
ylim(0,zlim)

