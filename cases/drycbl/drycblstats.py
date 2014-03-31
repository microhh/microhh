import numpy
import struct
import netCDF4

from pylab import *

start = 0
end   = 17
step  = 8

stats = netCDF4.Dataset("drycbl.default.0000000.nc","r")
t = stats.variables["t"][start:end]
z = stats.variables["z"][:]

b     = stats.variables["b"    ][start:end,:]
bsort = stats.variables["bsort"][start:end,:]

zsort = stats.variables["zsort"][start:end,:]

pe_total = stats.variables["pe_total"][start:end,:]
pe_bg    = stats.variables["pe_bg"   ][start:end,:]
pe_avail = stats.variables["pe_avail"][start:end,:]

tke = stats.variables["tke"][start:end,:]
#tkeref = 0.5*(stats.variables["u2"][start:end,:] + stats.variables["v2"][start:end,:] + 0.5*(stats.variables["w2"][start:end,0:-1]+stats.variables["w2"][start:end,1::]))
#ke  = stats.variables["ke" ][start:end,:]

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


