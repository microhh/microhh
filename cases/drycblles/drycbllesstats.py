import numpy
import struct
import netCDF4

from pylab import *

start = 0
end   = 19
step  = 2
plotens = False
dx = 6400./256

stats = netCDF4.Dataset("drycblles.0000000.nc","r")
t  = stats.variables["t"][start:end]
z  = stats.variables["z"][:]
zh = stats.variables["zh"][:]
st  = stats.variables["s"][start:end,:]
evisct = stats.variables["evisc"][start:end,:]
u2t = stats.variables["u2"][start:end,:]
v2t = stats.variables["v2"][start:end,:]
w2t = stats.variables["w2"][start:end,:]
s2t = stats.variables["s2"][start:end,:]
sturbt = stats.variables["sw"][start:end,:]
sdifft = stats.variables["sdiff"][start:end,:]
sfluxt = stats.variables["sflux"][start:end,:]
sgradt = stats.variables["sgrad"][start:end,:]

s = numpy.mean(st,0)
evisc = numpy.mean(evisct,0)

u2 = numpy.mean(u2t,0)
v2 = numpy.mean(v2t,0)
w2 = numpy.mean(w2t,0)
s2 = numpy.mean(s2t,0)
 
sturb = numpy.mean(sturbt,0)
sdiff = numpy.mean(sdifft,0)
sflux = numpy.mean(sfluxt,0)

ht     = zeros(t.size)
wstart = zeros(t.size)
Ret    = zeros(t.size)
etat   = zeros(t.size)
for n in range(t.size):
  hindex    = find(sgradt[n,:] == max(sgradt[n,:]))
  ht[n]     = z[hindex[0]]
  wstart[n] = ((9.81/300.)*sfluxt[n,0]*ht[n])**(1./3.)
  Ret[n]    = wstart[n]*ht[n] / max(evisct[n,:])
  etat[n]   = (max(evisct[n,:])**3./((9.81/300.)*sfluxt[n,0]))**.25
  

# enable LaTeX plotting
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

close('all')
figure()
for n in range(start,end,step):
  plot(st[n,:], z)

figure()
for n in range(start,end,step):
  plot(evisct[n,:], z)

figure()
for n in range(start,end,step):
  plot(u2t[n,:], z)

figure()
for n in range(start,end,step):
  plot(w2t[n,:], zh)

figure()
for n in range(start,end,step):
  plot(sfluxt[n,:], zh)

figure()
plot(numpy.mean(sfluxt[-5:-1,:],0), zh, 'b-' )
plot(numpy.mean(sturbt[-5:-1,:],0), zh, 'b--')
plot(numpy.mean(sdifft[-5:-1,:],0), zh, 'b:' )
ylim(0., 1500.)

tlog   = linspace(1000.,8000.,100)
reflog = 8.*tlog**.5
figure()
loglog(t, ht, 'bo' )
loglog(tlog, reflog, 'k--')
xlabel(r'time')
ylabel(r'h')

figure()
loglog(t, Ret     , 'bo', label=r'$Re$')
loglog(t, Ret**.75, 'ro', label=r'$h / \eta$')
xlabel(r'time')
legend(loc=0, frameon=False)

figure()
plot(t, etat/dx, 'bo', label=r'$\eta / \Delta x$')
xlabel(r'time')
legend(loc=0, frameon=False)
