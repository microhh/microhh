import numpy
import struct
import netCDF4

from pylab import *

stats = netCDF4.Dataset("drycbl.default.0000000.nc","r")

t = stats.variables["t"][:]
end   = t.size
start = t.size-1

z  = stats.variables["z"][:]
zh = stats.variables["zh"][:]

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
bpe_sum  = bpe_turb + bpe_visc + bpe_diss

ape_turb = pe_turb - bpe_turb
ape_visc = pe_visc - bpe_visc
ape_diss = -bpe_diss
ape_bous = pe_bous
ape_buoy = pe_buoy
ape_sum  = ape_turb + ape_visc + ape_diss + ape_bous + ape_buoy

# enable LaTeX plotting
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

close('all')
figure()
plot(tke_turb, z, 'b-', label='turb')
plot(tke_visc, z, 'g-', label='visc')
plot(tke_diss, z, 'r-', label='diss')
plot(tke_pres, z, 'm-', label='pres')
plot(tke_buoy, z, 'c-', label='buoy')
plot(tke_sum , z, 'k:', label='sum' )
xlabel(r'de$_k$/dt [m$^2$~s$^{-3}$]')
ylabel(r'z [m]')
legend(loc=0, frameon=False)
ylim(0,0.5)

figure()
plot(pe_turb, z, 'b-', label='turb')
plot(pe_visc, z, 'g-', label='visc')
plot(pe_bous, z, 'm-', label='bous')
plot(pe_buoy, z, 'c-', label='buoy')
plot(pe_sum , z, 'k:', label='sum' )
xlabel(r'de$_p$/dt [m$^2$~s$^{-3}$]')
ylabel(r'z [m]')
legend(loc=0, frameon=False)
ylim(0,0.5)

figure()
plot(bpe_turb, z, 'b-', label='turb')
plot(bpe_visc, z, 'g-', label='visc')
plot(bpe_diss, z, 'r-', label='diss')
plot(bpe_sum , z, 'k:', label='sum' )
xlabel(r'de$_p$/dt [m$^2$~s$^{-3}$]')
ylabel(r'z [m]')
legend(loc=0, frameon=False)
ylim(0,0.5)

figure()
plot(ape_turb, z, 'b-', label='turb')
plot(ape_visc, z, 'g-', label='visc')
plot(ape_bous, z, 'm-', label='bous')
plot(ape_buoy, z, 'c-', label='buoy')
plot(ape_diss, z, 'r-', label='diss')
plot(ape_sum , z, 'k:', label='sum' )
xlabel(r'de$_p$/dt [m$^2$~s$^{-3}$]')
ylabel(r'z [m]')
legend(loc=0, frameon=False)
ylim(0,0.5)

