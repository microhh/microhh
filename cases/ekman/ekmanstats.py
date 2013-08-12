from pylab import *
import netCDF4

# set the height
kmax  = 64
zsize = 5.

# analytical solution
ug = 1.
visc = 0.1
fc = 1.

gamma = (fc / (2.*visc))**.5
ustar = (ug*(visc*fc)**.5)**.5
h     = pi / gamma

# load the simulation data
stats32 = netCDF4.Dataset("ekman32/ekman.0000000.nc","r")
t32 = stats32.variables["t"][:]
z32 = stats32.variables["z"][:]
u32 = stats32.variables["u"][:,:]
v32 = stats32.variables["v"][:,:]
uref32 = zeros(z32.size)
vref32 = zeros(z32.size)
for k in range(z32.size):
  uref32[k] = ug*(1.- exp(-gamma*z32[k]) * cos(gamma*z32[k]))
  vref32[k] = ug*(    exp(-gamma*z32[k]) * sin(gamma*z32[k]))
dz32 = zsize / 32

stats64 = netCDF4.Dataset("ekman64/ekman.0000000.nc","r")
t64 = stats64.variables["t"][:]
z64 = stats64.variables["z"][:]
u64 = stats64.variables["u"][:,:]
v64 = stats64.variables["v"][:,:]
uref64 = zeros(z64.size)
vref64 = zeros(z64.size)
for k in range(z64.size):
  uref64[k] = ug*(1.- exp(-gamma*z64[k]) * cos(gamma*z64[k]))
  vref64[k] = ug*(    exp(-gamma*z64[k]) * sin(gamma*z64[k]))
dz64 = zsize / 64

stats128 = netCDF4.Dataset("ekman128/ekman.0000000.nc","r")
t128 = stats128.variables["t"][:]
z128 = stats128.variables["z"][:]
u128 = stats128.variables["u"][:,:]
v128 = stats128.variables["v"][:,:]
uref128 = zeros(z128.size)
vref128 = zeros(z128.size)
for k in range(z128.size):
  uref128[k] = ug*(1.- exp(-gamma*z128[k]) * cos(gamma*z128[k]))
  vref128[k] = ug*(    exp(-gamma*z128[k]) * sin(gamma*z128[k]))
dz128 = zsize / 128

u32error  = sqrt(dz32 *sum((u32 [-1,:] - uref32 )**2.))
u64error  = sqrt(dz64 *sum((u64 [-1,:] - uref64 )**2.))
u128error = sqrt(dz128*sum((u128[-1,:] - uref128)**2.))
v32error  = sqrt(dz32 *sum((v32 [-1,:] - vref32 )**2.))
v64error  = sqrt(dz64 *sum((v64 [-1,:] - vref64 )**2.))
v128error = sqrt(dz128*sum((v128[-1,:] - vref128)**2.))
print("L2 errors u = ", u32error, u64error, u128error)
print("L2 errors v = ", v32error, v64error, v128error)
print("Convergence in u = ", (log(u128error) - log(u32error)) / (log(zsize/128.) - log(zsize/32.)))
print("Convergence in v = ", (log(v128error) - log(v32error)) / (log(zsize/128.) - log(zsize/32.)))

close('all')
# ekman spiral
figure()
plot(u32 [-1,:], v32 [-1,:])
plot(u64 [-1,:], v64 [-1,:])
plot(u128[-1,:], v128[-1,:])

# vertical profiles
figure()
subplot(121)
plot(u32 [-1,:], z32 , 'b-')
plot(u64 [-1,:], z64 , 'g-')
plot(u128[-1,:], z128, 'r-')
subplot(122)
plot(v32 [-1,:], z32 , 'b-')
plot(v64 [-1,:], z64 , 'g-')
plot(v128[-1,:], z128, 'r-')

# errors
figure()
subplot(121)
plot(u32 [-1,:] - uref32 , z32 , 'b-')
plot(u64 [-1,:] - uref64 , z64 , 'g-')
plot(u128[-1,:] - uref128, z128, 'r-')
subplot(122)
plot(v32 [-1,:] - vref32 , z32 , 'b-')
plot(v64 [-1,:] - vref64 , z64 , 'g-')
plot(v128[-1,:] - vref128, z128, 'r-')

