import numpy
import struct
import netCDF4

from pylab import *

stats = netCDF4.Dataset("drycbl.default.0000000.nc","r")

t = stats.variables["t"][:]
end   = t.size
start = t.size-20
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

tke = numpy.dot(stats.variables["tke"][:,:], dz)
pe  = numpy.dot(stats.variables["pe" ][:,:], dz)
ape = numpy.dot(stats.variables["ape"][:,:], dz)
bpe = numpy.dot(stats.variables["bpe"][:,:], dz)

bdzstardt = zeros(stats.variables["zsort"][:,:].shape)
zsort = stats.variables["zsort"][:,:]
b     = stats.variables["b"]    [:,:]
for k in range(z.size):
  bdzstardt[:,k] = -b[:,k] * numpy.gradient(zsort[:,k], dt)
del(b, zsort)

tke_turb = average(stats.variables["tke_turb"][start:end,:], 0)
tke_visc = average(stats.variables["tke_visc"][start:end,:], 0)
tke_pres = average(stats.variables["tke_pres"][start:end,:], 0)
tke_diss = average(stats.variables["tke_diss"][start:end,:], 0)
tke_buoy = average(stats.variables["tke_buoy"][start:end,:], 0)
tke_sum = tke_turb + tke_visc + tke_pres + tke_diss + tke_buoy

pe_turb = average(stats.variables["pe_turb" ][start:end,:], 0)
pe_visc = average(stats.variables["pe_visc" ][start:end,:], 0)
pe_bous = average(stats.variables["pe_bous" ][start:end,:], 0)
pe_buoy = average(stats.variables["tke_buoy"][start:end,:], 0)*-1.
pe_sum  = pe_turb + pe_visc + pe_bous + pe_buoy

bpe_turb = average(stats.variables["bpe_turb"][start:end,:], 0)
bpe_visc = average(stats.variables["bpe_visc"][start:end,:], 0)
bpe_diss = average(stats.variables["bpe_diss"][start:end,:], 0)
bpe_dzdt = average(bdzstardt[start:end,:], 0)
bpe_sum  = bpe_turb + bpe_visc + bpe_diss + bpe_dzdt

ape_turb = pe_turb - bpe_turb
ape_visc = pe_visc - bpe_visc
ape_diss = -bpe_diss
ape_bous = pe_bous
ape_buoy = pe_buoy
ape_dzdt = -bpe_dzdt
ape_sum  = ape_turb + ape_visc + ape_diss + ape_bous + ape_buoy + ape_dzdt

tke_turb_time = numpy.dot(stats.variables["tke_turb"][:,:], dz)
tke_visc_time = numpy.dot(stats.variables["tke_visc"][:,:], dz)
tke_pres_time = numpy.dot(stats.variables["tke_pres"][:,:], dz)
tke_diss_time = numpy.dot(stats.variables["tke_diss"][:,:], dz)
tke_buoy_time = numpy.dot(stats.variables["tke_buoy"][:,:], dz)
tke_sum_time  = tke_turb_time + tke_visc_time + tke_pres_time + tke_diss_time + tke_buoy_time

pe_turb_time = numpy.dot(stats.variables["pe_turb" ][:,:], dz)
pe_visc_time = numpy.dot(stats.variables["pe_visc" ][:,:], dz)
pe_buoy_time = numpy.dot(stats.variables["tke_buoy"][:,:], dz) * -1.
pe_bous_time = numpy.dot(stats.variables["pe_bous" ][:,:], dz)
pe_sum_time  = pe_turb_time + pe_visc_time + pe_buoy_time + pe_bous_time

ape_turb_time = numpy.dot(stats.variables["pe_turb" ][:,:] - stats.variables["bpe_turb"][:,:], dz)
ape_visc_time = numpy.dot(stats.variables["pe_visc" ][:,:] - stats.variables["bpe_visc"][:,:], dz)
#ape_diss_time = numpy.dot(stats.variables["bpe_diss"][:,:], dz) * -1.
# CvH HACK bug to fix
ape_diss_time = numpy.dot(stats.variables["bpe_diss"][:,0:-1], dz[0:-1]) * -1.
ape_buoy_time = numpy.dot(stats.variables["tke_buoy"][:,:], dz) * -1.
ape_bous_time = numpy.dot(stats.variables["pe_bous" ][:,:], dz)
ape_sum_time  = ape_turb_time + ape_visc_time + ape_diss_time + ape_buoy_time + ape_bous_time

bpe_turb_time = numpy.dot(stats.variables["bpe_turb"][:,:], dz)
bpe_visc_time = numpy.dot(stats.variables["bpe_visc"][:,:], dz)
bpe_diss_time = numpy.dot(stats.variables["bpe_diss"][:,0:-1], dz[0:-1])
bpe_sum_time  = bpe_turb_time + bpe_visc_time + bpe_diss_time

eff0 = (ape_diss_time - ape_bous_time) / (ape_diss_time - ape_bous_time + tke_diss_time)
eff1 = 1 - ape_bous_time / ape_diss_time
eff2 = ape_diss_time / (ape_diss_time + ape_bous_time)

# enable LaTeX plotting
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

zlim = 1.8

close('all')
figure()
plot(tke_turb/B0, z/henc, 'b-', label='turb')
plot(tke_visc/B0, z/henc, 'g-', label='visc')
plot(tke_diss/B0, z/henc, 'r-', label='diss')
plot(tke_pres/B0, z/henc, 'm-', label='pres')
plot(tke_buoy/B0, z/henc, 'c-', label='buoy')
plot(tke_sum /B0, z/henc, 'k:', label='sum' )
xlabel(r'de$_k$/dt / B$_0$')
ylabel(r'z/h$_{enc}$')
legend(loc=0, frameon=False)
ylim(0,zlim)

figure()
plot(pe_turb/B0, z/henc, 'b-', label='turb')
plot(pe_visc/B0, z/henc, 'g-', label='visc')
plot(pe_bous/B0, z/henc, 'm-', label='bous')
plot(pe_buoy/B0, z/henc, 'c-', label='buoy')
plot(pe_sum /B0, z/henc, 'k:', label='sum' )
xlabel(r'de$_p$/dt / B$_0$')
ylabel(r'z/h$_{enc}$')
legend(loc=0, frameon=False)
ylim(0,zlim)

figure()
plot(bpe_turb/B0, z/henc, 'b-', label='turb')
plot(bpe_visc/B0, z/henc, 'g-', label='visc')
plot(bpe_diss/B0, z/henc, 'r-', label='diss')
plot(bpe_dzdt/B0, z/henc, 'c-', label='dzdt')
plot(bpe_sum /B0, z/henc, 'k:', label='sum' )
xlabel(r'de$_b$/dt / B$_0$')
ylabel(r'z/h$_{enc}$')
legend(loc=0, frameon=False)
ylim(0,zlim)

figure()
plot(ape_turb/B0, z/henc, 'b-', label='turb')
plot(ape_visc/B0, z/henc, 'g-', label='visc')
plot(ape_bous/B0, z/henc, 'm-', label='bous')
plot(ape_buoy/B0, z/henc, 'c-', label='buoy')
plot(ape_dzdt/B0, z/henc, 'y-', label='dzdt')
plot(ape_diss/B0, z/henc, 'r-', label='diss')
plot(ape_sum /B0, z/henc, 'k:', label='sum' )
xlabel(r'de$_a$/dt / B$_0$')
ylabel(r'z/h$_{enc}$')
legend(loc=0, frameon=False)
ylim(0,zlim)

figure()
plot(tenct, tke_turb_time / (henct*B0), 'b-', label='turb')
plot(tenct, tke_visc_time / (henct*B0), 'g-', label='visc')
plot(tenct, tke_diss_time / (henct*B0), 'r-', label='diss')
plot(tenct, tke_buoy_time / (henct*B0), 'c-', label='buoy')
plot(tenct, tke_pres_time / (henct*B0), 'm-', label='pres')
plot(tenct, tke_sum_time  / (henct*B0), 'k:', label='sum' )
xlabel(r'h$_{enc}$/L$_0$')
ylabel(r'de$_k$/dt / B$_0$')
legend(loc=0, frameon=False)

figure()
plot(tenct, pe_turb_time / (henct*B0), 'b-', label='turb')
plot(tenct, pe_visc_time / (henct*B0), 'g-', label='visc')
plot(tenct, pe_buoy_time / (henct*B0), 'c-', label='buoy')
plot(tenct, pe_bous_time / (henct*B0), 'm-', label='bous')
plot(tenct, pe_sum_time  / (henct*B0), 'k:', label='sum' )
xlabel(r'h$_{enc}$/L$_0$')
ylabel(r'de$_p$/dt / B$_0$')
legend(loc=0, frameon=False)

figure()
plot(tenct, ape_turb_time / (henct*B0), 'b-', label='turb')
plot(tenct, ape_visc_time / (henct*B0), 'g-', label='visc')
plot(tenct, ape_diss_time / (henct*B0), 'r-', label='diss')
plot(tenct, ape_buoy_time / (henct*B0), 'c-', label='buoy')
plot(tenct, ape_bous_time / (henct*B0), 'm-', label='bous')
plot(tenct, ape_sum_time  / (henct*B0), 'k:', label='sum' )
xlabel(r'h$_{enc}$/L$_0$')
ylabel(r'de$_a$/dt / B$_0$')
legend(loc=0, frameon=False)

figure()
plot(tenct, eff0, 'b-', label='$\eta$ (Gayen)')
plot(tenct, eff1, 'g-', label='$\eta$ (Gayen2)')
plot(tenct, eff2, 'r-', label='$\eta$ (Scotti)')
xlabel(r'h$_{enc}$/L$_0$')
ylabel(r'$\eta$ [-]')
legend(loc=0, frameon=False)

figure()
plot(tenct, tke_sum_time / (henct*B0), 'b-', label='e$_k$')
plot(tenct, pe_sum_time  / (henct*B0), 'g-', label='e$_p$')
plot(tenct, (tke_sum_time + pe_sum_time) / (henct*B0), 'k:', label='e$_k$ + e$_p$')
xlabel(r'h$_{enc}$/L$_0$')
ylabel(r'de/dt / B$_0$')
legend(loc=0, frameon=False)

figure()
plot(tenct, tke_sum_time / (henct*B0), 'b-', label='e$_k$')
plot(tenct, ape_sum_time / (henct*B0), 'r-', label='e$_a$')
plot(tenct, bpe_sum_time / (henct*B0), 'c-', label='e$_b$')
plot(tenct, (tke_sum_time + ape_sum_time + bpe_sum_time) / (henct*B0), 'k:', label='e$_k$ + e$_a$ + e$_b$')
xlabel(r'h$_{enc}$/L$_0$')
ylabel(r'de/dt / B$_0$')
legend(loc=0, frameon=False)

figure()
plot(tenct, (tke       ) / (henct*wenct**2.), 'b-', label='tke')
plot(tenct, (pe-pe[0]  ) / (henct*wenct**2.), 'g-', label='pe' )
plot(tenct, (ape       ) / (henct*wenct**2.), 'r-', label='ape')
plot(tenct, (bpe-bpe[0]) / (henct*wenct**2.), 'm-', label='bpe')
plot(tenct, tke/ape                         , 'k:', label='tke/ape')
xlabel(r'h$_{enc}$/L$_0$')
ylabel(r'E/(h$_{enc}$ w$_{enc}^2)$')
legend(loc=0, frameon=False)

