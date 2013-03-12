from pylab import *
from microhh import *

t = 500

itot = 512
jtot = 1
ktot = 512

data = microhh(t, itot, jtot, ktot)

figure()
pcolormesh(data.xh, data.z, data.u[:,0,:])
xlim(min(data.xh), max(data.xh))
ylim(min(data.z ), max(data.z ))
xlabel('x')
ylabel('z')
title('streamwise velocity')
colorbar()

"""
figure()
pcolormesh(data.xh, data.y, data.u[0,:,:])
xlim(min(data.xh), max(data.xh))
ylim(min(data.y ), max(data.y ))
xlabel('x')
ylabel('y')
title('streamwise velocity')
colorbar()

figure(figsize=(16,6))
subplot(131)
pcolormesh(data.y, data.z, data.u[:,:,0])
xlim(min(data.y), max(data.y))
ylim(min(data.z), max(data.z))
xlabel('y')
ylabel('z')
title('streamwise velocity')
colorbar()
subplot(132)
pcolormesh(data.yh, data.z, data.v[:,:,0])
xlim(min(data.yh), max(data.yh))
ylim(min(data.z ), max(data.z ))
xlabel('y')
ylabel('z')
title('spanwise velocity')
colorbar()
subplot(133)
pcolormesh(data.y, data.zh, data.w[:,:,0])
xlim(min(data.y ), max(data.y ))
ylim(min(data.zh), max(data.zh))
xlabel('y ')
ylabel('zh')
title('vertical velocity')
colorbar()

figure()
pcolormesh(data.x, data.z, data.v[:,0,:])
xlim(min(data.x), max(data.x))
ylim(min(data.z), max(data.z))
xlabel('x')
ylabel('z')
title('spanwise velocity')
colorbar()

"""
figure()
pcolormesh(data.x, data.zh, data.w[:,0,:])
xlim(min(data.x ), max(data.x ))
ylim(min(data.zh), max(data.zh))
xlabel('x')
ylabel('z')
title('vertical velocity')
colorbar()

figure()
pcolormesh(data.x, data.z, data.p[:,0,:])
xlim(min(data.x), max(data.x))
ylim(min(data.z), max(data.z))
xlabel('x')
ylabel('z')
title('pressure')
colorbar()

figure()
pcolormesh(data.x, data.z, data.s[:,0,:])
xlim(min(data.x), max(data.x))
ylim(min(data.z), max(data.z))
xlabel('x')
ylabel('z')
title('scalar')
colorbar()

"""
umean = average(average(data.u,1),1)
vmean = average(average(data.v,1),1)
wmean = average(average(data.w,1),1)
pmean = average(average(data.p,1),1)

figure()
plot(umean, data.z , label='u')
plot(vmean, data.z , label='v')
plot(wmean, data.zh, label='w')
xlabel('velocity')
ylabel('height')
legend(loc=0, frameon=False)

figure()
plot(pmean, data.z)
xlabel('modified pressure')
ylabel('height')
"""

smean = average(average(data.s,1),1)

figure()
plot(smean, data.z, 'bo-')
xlabel('scalar')
ylabel('height')
