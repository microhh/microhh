import numpy
import struct
import netCDF4

from pylab import *

start = 0
end   = 15
step  = 5

stats = netCDF4.Dataset("drycbl.default.0000000.nc","r")
t = stats.variables["t"][start:end]
z = stats.variables["z"][:]

bt     = stats.variables["b"    ][start:end,:]
bsortt = stats.variables["bsort"][start:end,:]

# enable LaTeX plotting
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

close('all')
figure()
for n in range(start,end,step):
  plot(bt    [n,:], z, 'b-')
  plot(bsortt[n,:], z, 'r-')
xlabel(r'b [m~s$^{-2}$]')
ylabel(r'z [m]')

