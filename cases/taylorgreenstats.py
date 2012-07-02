from pylab import *
from microhh import *

t=200

itot = 96
jtot = 1
ktot = 48

data = microhh(t, itot, jtot, ktot)

figure()
pcolormesh(data.xh, data.z, data.u[:,0,:])
xlim(min(data.xh), max(data.xh))
ylim(min(data.z ), max(data.z ))
xlabel('x')
ylabel('z')
title('streamwise velocity')
colorbar()

figure()
pcolormesh(data.x, data.z, data.v[:,0,:])
xlim(min(data.x), max(data.x))
ylim(min(data.z), max(data.z))
xlabel('x')
ylabel('z')
title('spanwise velocity')
colorbar()

figure()
pcolormesh(data.x, data.zh, data.w[:,0,:])
xlim(min(data.x ), max(data.x ))
ylim(min(data.zh), max(data.zh))
xlabel('x')
ylabel('z')
title('vertical velocity')
colorbar()

figure()
pcolormesh(data.x, data.zh, data.p[:,0,:])
xlim(min(data.x), max(data.x))
ylim(min(data.z), max(data.z))
xlabel('x')
ylabel('z')
title('modified pressure')
colorbar()

