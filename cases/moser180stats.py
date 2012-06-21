import numpy
import struct
import netCDF4

from pylab import *

nx = 256
ny = 192
nz = 128
iter = 7000

# read the 3d data
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

fin = open("u.{:07d}".format(iter),"rb")
raw = fin.read(n*8)
tmp = numpy.array(struct.unpack('{}d'.format(n), raw))
u   = tmp.reshape((nz, ny, nx))
fin.close()

fin = open("v.{:07d}".format(iter),"rb")
raw = fin.read(n*8)
tmp = numpy.array(struct.unpack('{}d'.format(n), raw))
v   = tmp.reshape((nz, ny, nx))
fin.close()

fin = open("w.{:07d}".format(iter),"rb")
raw = fin.read(n*8)
tmp = numpy.array(struct.unpack('{}d'.format(n), raw))
w   = tmp.reshape((nz, ny, nx))
fin.close()

del(raw)
del(tmp)

# read Moser's data
Mosermean = numpy.loadtxt("chan180.means", skiprows=25)
Moserrey  = numpy.loadtxt("chan180.reystress", skiprows=25)

yplusMoser = Mosermean[:,1]
uavgMoser  = Mosermean[:,2]
uvarMoser  = Moserrey[:,2]
vvarMoser  = Moserrey[:,3]
wvarMoser  = Moserrey[:,4]

# process the 3d data
uavg  = numpy.mean(numpy.mean(u,2),1)
uvar  = numpy.zeros(nz)
vvar  = numpy.zeros(nz)
wvar  = numpy.zeros(nz)

for k in range(nz):
  uvar[k] = numpy.var(u[k,:,:] - uavg[k])
  vvar[k] = numpy.var(v[k,:,:])
  wvar[k] = numpy.var(w[k,:,:])

visc  = 1.0e-5
ustar = (visc * uavg[0] / z[0])**0.5

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
semilogx(yplus[starty:endy], uavg[starty:endy] / ustar, 'b-')
semilogx(yplusMoser, uavgMoser, 'k--', label="Moser")
semilogx(ypluslin, ulin, 'k:')
semilogx(ypluslog, ulog, 'k:')
xlabel('y+')
ylabel('u+')
legend(loc=2, frameon=False)
grid()
axis([0.3, 200, 0, 22])

figure()
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
 
