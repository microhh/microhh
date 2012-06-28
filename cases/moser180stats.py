import numpy
import struct
import netCDF4

from pylab import *

nx = 256
ny = 192
nz = 128

iter = 0
nt   = 1

# read Moser's data
Mosermean = numpy.loadtxt("chan180.means", skiprows=25)
Moserrey  = numpy.loadtxt("chan180.reystress", skiprows=25)

yplusMoser = Mosermean[:,1]
uavgMoser  = Mosermean[:,2]
uvarMoser  = Moserrey[:,2]
vvarMoser  = Moserrey[:,3]
wvarMoser  = Moserrey[:,4]

# read the grid data
n = nx*ny*nz

fin = open("grid.{:07d}".format(0),"rb")
raw = fin.read(nx*8)
x   = numpy.array(struct.unpack('{}d'.format(nx), raw))
raw = fin.read(nx*8)
xh  = numpy.array(struct.unpack('{}d'.format(nx), raw))
raw = fin.read(ny*8)
y   = numpy.array(struct.unpack('{}d'.format(ny), raw))
raw = fin.read(ny*8)
yh  = numpy.array(struct.unpack('{}d'.format(ny), raw))
raw = fin.read(nz*8)
z   = numpy.array(struct.unpack('{}d'.format(nz), raw))
raw = fin.read(nz*8)
zh  = numpy.array(struct.unpack('{}d'.format(nz), raw))
fin.close()

# read the 3d data and process it
uavgt = numpy.zeros((nt, nz))
vavgt = numpy.zeros((nt, nz))

uvart = numpy.zeros((nt, nz))
vvart = numpy.zeros((nt, nz))
wvart = numpy.zeros((nt, nz))

for t in range(nt):
  prociter = iter + 500*t
  print("Processing iter = {:07d}".format(prociter))

  fin = open("u.{:07d}".format(prociter),"rb")
  raw = fin.read(n*8)
  tmp = numpy.array(struct.unpack('{}d'.format(n), raw))
  u   = tmp.reshape((nz, ny, nx))
  fin.close()
  
  fin = open("v.{:07d}".format(prociter),"rb")
  raw = fin.read(n*8)
  tmp = numpy.array(struct.unpack('{}d'.format(n), raw))
  v   = tmp.reshape((nz, ny, nx))
  fin.close()
  
  fin = open("w.{:07d}".format(prociter),"rb")
  raw = fin.read(n*8)
  tmp = numpy.array(struct.unpack('{}d'.format(n), raw))
  w   = tmp.reshape((nz, ny, nx))
  fin.close()
  
  del(raw)
  del(tmp)

  uavgt[t,:] = numpy.mean(numpy.mean(u,2),1)
  vavgt[t,:] = numpy.mean(numpy.mean(v,2),1)

  for k in range(nz):
    uvart[t,k] = numpy.var(u[k,:,:] - uavgt[t,k])
    vvart[t,k] = numpy.var(v[k,:,:] - vavgt[t,k])
    wvart[t,k] = numpy.var(w[k,:,:])

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
for t in range(nt):
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
for t in range(nt):
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
 
