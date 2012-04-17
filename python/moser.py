from pylab import *
from microhh import *

t = 0

itot = 64
jtot = 64
ktot = 64

data = microhh(t, itot, jtot, ktot)

figure(1)
pcolormesh(data.xh, data.z, data.u[:,1,:])
xlim(min(data.xh), max(data.xh))
ylim(min(data.z ), max(data.z ))

figure(2)
pcolormesh(data.x, data.z, data.v[:,1,:])
xlim(min(data.x), max(data.x))
ylim(min(data.z), max(data.z))

figure(3)
pcolormesh(data.x, data.zh, data.w[:,1,:])
xlim(min(data.x ), max(data.x ))
ylim(min(data.zh), max(data.zh))

figure(4)
pcolormesh(data.x, data.z, data.p[:,1,:])
xlim(min(data.x), max(data.x))
ylim(min(data.z), max(data.z))

figure(5)
pcolormesh(data.x, data.z, data.s[:,1,:])
xlim(min(data.x), max(data.x))
ylim(min(data.z), max(data.z))

##### some means ######
umean = average(average(data.u,1),1)

uvar  = average(average(data.u**2.,1),1)
vvar  = average(average(data.v**2.,1),1)
wvar  = average(average(data.w**2.,1),1)

figure(6)
plot(umean, data.z)

uvar = zeros(ktot)
for k in range(ktot):
  uplane = average(data.u[k,:,:])
  uvar[k] = average((data.u[k,:,:] - uplane)**2.)

figure(7)
plot(uvar, data.z, label='uvar')
plot(vvar, data.z, label='vvar')
plot(wvar, data.zh, label='wvar')
legend(loc=0, frameon=False)

