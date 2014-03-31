import numpy
import struct
import netCDF4

from pylab import *

start = 0
end   = 17
step  = 8

stats = netCDF4.Dataset("drycbl.default.0000000.nc","r")

t = stats.variables["t"][:]

z  = stats.variables["z"][:]
zh = stats.variables["zh"][:]

b     = stats.variables["b"    ][:,:]
bsort = stats.variables["bsort"][:,:]

zsort = stats.variables["zsort"][:,:]

pe_total = stats.variables["pe_total"][:,:]
pe_bg    = stats.variables["pe_bg"   ][:,:]
pe_avail = stats.variables["pe_avail"][:,:]

tke = stats.variables["tke"][:,:]
#tkeref = 0.5*(stats.variables["u2"][start:end,:] + stats.variables["v2"][start:end,:] + 0.5*(stats.variables["w2"][start:end,0:-1]+stats.variables["w2"][start:end,1::]))
#ke  = stats.variables["ke" ][start:end,:]

dz = zh[1::] - zh[0:-1]

pe_total_sum = zeros(t.size)
tke_sum      = zeros(t.size)

for n in range(t.size):
  pe_total_sum [n] = numpy.sum(dz[:]*pe_total[n,:])
  tke_sum[n]       = numpy.sum(dz[:]*tke     [n,:])

pe_botflux = stats.variables["bflux"][:,0] * t
pe_topflux = stats.variables["bflux"][:,-1] * t

# enable LaTeX plotting
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

close('all')
figure()
for n in range(start,end,step):
  plot(b    [n,:], z, 'b-')
  plot(bsort[n,:], z, 'r-')
xlabel(r'b [m~s$^{-2}$]')
ylabel(r'z [m]')

figure()
for n in range(start,end,step):
  plot(zsort[n,:], z, 'b-')
xlabel(r'zsort [m]')
ylabel(r'z [m]')

figure()
for n in range(start,end,step):
  plot(pe_total[n,:], z, 'b-')
  plot(pe_bg   [n,:], z, 'g-')
  plot(pe_avail[n,:], z, 'r-')
xlabel(r'PE [m$^2$~s$^{-2}$]')
ylabel(r'z  [m]')

figure()
for n in range(start,end,step):
  plot(tke   [n,:], z, 'b-')
xlabel(r'TKE [m$^2$~s$^{-2}$]')
ylabel(r'z [m]')

figure()
plot(t, pe_total_sum - pe_total_sum[0], label='PE' )
plot(t, pe_botflux, label='PE bot')
plot(t, pe_topflux, label='PE top')
plot(t, tke_sum, label='TKE')
plot(t, pe_total_sum - pe_total_sum[0] + tke_sum, label='PE + TKE')
legend(loc=0, frameon=False)

