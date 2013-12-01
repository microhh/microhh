from pylab import *
import netCDF4

# load the data
stats = netCDF4.Dataset("ekman.0000000.nc","r")
t = stats.variables["t"][:]
z = stats.variables["z"][:]
u = stats.variables["u"][:,:]
v = stats.variables["v"][:,:]

# analytical solution
zsize = 5.
ug    = 1.
visc  = 0.1
fc    = 1.

gamma = (fc / (2.*visc))**.5
ustar = (ug*(visc*fc)**.5)**.5
h     = pi / gamma

uref = zeros(z.size)
vref = zeros(z.size)
for k in range(z.size):
  uref[k] = ug*(1.- exp(-gamma*z[k]) * cos(gamma*z[k]))
  vref[k] = ug*(    exp(-gamma*z[k]) * sin(gamma*z[k]))
dz = zsize / z.size 

uerror = sqrt(dz*sum((u[-1,:] - uref)**2.))
verror = sqrt(dz*sum((v[-1,:] - vref)**2.))

close('all')
# ekman spiral
figure()
plot(u[-1,:], v[-1,:])
plot(uref, vref, 'k--')
xlabel('u [m/s]')
ylabel('v [m/s]')

# vertical profiles
figure()
subplot(121)
plot(u[-1,:], z)
plot(uref[:], z, 'k--')
xlabel('u [m/s]')
ylabel('z [m]')
subplot(122)
plot(v[-1,:], z)
plot(vref[:], z, 'k--')
xlabel('v [m/s]')
ylabel('z [m]')

