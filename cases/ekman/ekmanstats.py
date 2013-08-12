from pylab import *
import netCDF4

# set the height
kmax  = 96
zsize = 4.

# analytical solution
ug = 1.
visc = 0.1
fc = 1.

gamma = (fc / (2.*visc))**.5
ustar = (ug*(visc*fc)**.5)**.5
h     = pi / gamma

uref = zeros(kmax)
vref = zeros(kmax)
Uref = zeros(kmax)

dz = zsize / kmax
z  = linspace(0.5*dz, zsize-0.5*dz, kmax)

for k in range(kmax):
  uref[k] = ug*(1.- exp(-gamma*z[k]) * cos(gamma*z[k]))
  vref[k] = ug*(    exp(-gamma*z[k]) * sin(gamma*z[k]))
  Uref[k] = (uref[k]**2. + vref[k]**2.)**.5

# load the simulation data
stats = netCDF4.Dataset("ekman.0000000.nc","r")
t = stats.variables["t"][:]
z = stats.variables["z"][:]
u = stats.variables["u"][:,:]
v = stats.variables["v"][:,:]

close('all')
figure()
plot(u[-1,:], v[-1,:])
plot(uref, vref, 'k--')

#figure()
#plot(u/ug,z)
#plot(v/ug,z)
#plot(U/ug,z)
