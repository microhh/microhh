from pylab import *
from microhh import *

t = 0

itot = 64
jtot = 64
ktot = 64

igc = 1
jgc = 1
kgc = 1

data = microhh(t, itot, jtot, ktot, igc, kgc, jgc)

figure(1)
pcolormesh(data.xh, data.z, data.u[:,1,:])
xlim(min(data.xh), max(data.xh))
ylim(min(data.z ), max(data.z ))
xlabel('x')
ylabel('z')
colorbar()
title('streamwise velocity u')

figure(2)
pcolormesh(data.x, data.z, data.v[:,1,:])
xlim(min(data.x), max(data.x))
ylim(min(data.z), max(data.z))
xlabel('x')
ylabel('z')
colorbar()
title('spanwise velocity v')

figure(3)
pcolormesh(data.x, data.zh, data.w[:,1,:])
xlim(min(data.x ), max(data.x ))
ylim(min(data.zh), max(data.zh))
xlabel('x')
ylabel('z')
colorbar()
title('vertical velocity w')

figure(4)
pcolormesh(data.x, data.z, data.p[:,1,:])
xlim(min(data.x), max(data.x))
ylim(min(data.z), max(data.z))
xlabel('x')
ylabel('z')
colorbar()
title('modified pressure p')

figure(5)
pcolormesh(data.x, data.z, data.s[:,1,:])
xlim(min(data.x), max(data.x))
ylim(min(data.z), max(data.z))
xlabel('x')
ylabel('z')
colorbar()
title('scalar s')

##### some means ######
umean = average(average(data.u,1),1)
vmean = average(average(data.v,1),1)
wmean = average(average(data.w,1),1)

uvar = zeros(ktot+2*kgc)
for k in range(ktot+2*kgc):
  uvar[k] = average((data.u[k,:,:] - umean[k])**2.)

vvar  = average(average(data.v**2.,1),1)
wvar  = average(average(data.w**2.,1),1)

figure(6)
plot(umean[1:ktot+1], data.z [1:ktot+1], label='u')
plot(vmean[1:ktot+1], data.z [1:ktot+1], label='v')
plot(wmean[2:ktot+1], data.zh[2:ktot+1], label='w')
xlabel('velocity')
ylabel('z')

uvar = zeros(ktot+2*kgc)
for k in range(ktot+2*kgc):
  uplane = average(data.u[k,:,:])
  uvar[k] = average((data.u[k,:,:] - uplane)**2.)

figure(7)
plot(uvar[1:ktot+1], data.z [1:ktot+1], label='uvar')
plot(vvar[1:ktot+1], data.z [1:ktot+1], label='vvar')
plot(wvar[2:ktot+1], data.zh[2:ktot+1], label='wvar')
xlabel('velocity variance')
ylabel('z')
legend(loc=0, frameon=False)

