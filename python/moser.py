from pylab import *
from microhh import *

t = 3000

itot = 96
jtot = 64
ktot = 48

data = microhh(t, itot, jtot, ktot)

figure(1)
pcolormesh(data.xh, data.z, data.u[:,1,:])
xlim(min(data.xh), max(data.xh))
ylim(min(data.z ), max(data.z ))
xlabel('x')
ylabel('z')
title('streamwise velocity')
colorbar()

figure(2)
pcolormesh(data.x, data.z, data.v[:,1,:])
xlim(min(data.x), max(data.x))
ylim(min(data.z), max(data.z))
xlabel('x')
ylabel('z')
title('spanwise velocity')

figure(3)
pcolormesh(data.x, data.zh, data.w[:,1,:])
xlim(min(data.x ), max(data.x ))
ylim(min(data.zh), max(data.zh))
xlabel('x')
ylabel('z')
title('vertical velocity')

figure(4)
pcolormesh(data.x, data.z, data.p[:,1,:])
xlim(min(data.x), max(data.x))
ylim(min(data.z), max(data.z))
xlabel('x')
ylabel('z')
title('pressure')

figure(5)
pcolormesh(data.x, data.z, data.s[:,1,:])
xlim(min(data.x), max(data.x))
ylim(min(data.z), max(data.z))
xlabel('x')
ylabel('z')
title('scalar')

##### some means ######
umean = average(average(data.u,1),1)
vmean = average(average(data.v,1),1)
wmean = average(average(data.w,1),1)

pmean = average(average(data.p,1),1)

uvar  = average(average(data.u**2.,1),1)
vvar  = average(average(data.v**2.,1),1)
wvar  = average(average(data.w**2.,1),1)

figure(6)
plot(umean, data.z , label='u')
plot(vmean, data.z , label='v')
plot(wmean, data.zh, label='w')
xlabel('velocity')
ylabel('height')
legend(loc=0, frameon=False)

figure(7)
plot(pmean, data.z)
xlabel('modified pressure')
ylabel('height')

uvar = zeros(ktot)
for k in range(ktot):
  uplane = average(data.u[k,:,:])
  uvar[k] = average((data.u[k,:,:] - uplane)**2.)

figure(8)
plot(uvar, data.z , label='u')
plot(vvar, data.z , label='v')
plot(wvar, data.zh, label='w')
xlabel('variance')
ylabel('height')
legend(loc=0, frameon=False)

