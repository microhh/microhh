from pylab import *
from taylorgreenfunc import *

t    = 400
time = 1.
visc = (8.*pi**2. * 100.)**(-1.)

data32  = microhh(t,  32,  16, 'taylorgreen32' )
data64  = microhh(t,  64,  32, 'taylorgreen64' )
data128 = microhh(t, 128,  64, 'taylorgreen128')
data256 = microhh(t, 256, 128, 'taylorgreen256')

ref32  = getref(data32 .x, data32 .xh, data32 .z, data32 .zh, visc, time)
ref64  = getref(data64 .x, data64 .xh, data64 .z, data64 .zh, visc, time)
ref128 = getref(data128.x, data128.xh, data128.z, data128.zh, visc, time)
ref256 = getref(data256.x, data256.xh, data256.z, data256.zh, visc, time)

err32  = geterror(data32 , ref32 )
err64  = geterror(data64 , ref64 )
err128 = geterror(data128, ref128)
err256 = geterror(data256, ref256)

ns    = array([32, 64, 128, 256])
dxs   = 1./ns
errsu = array([err32.u, err64.u, err128.u, err256.u])
errsw = array([err32.w, err64.w, err128.w, err256.w])
errsp = array([err32.p, err64.p, err128.p, err256.p])

print('errors p', errsp)
if(t > 0):
  print('convergence u', (log(errsu[-1])-log(errsu[0])) / (log(dxs[-1])-log(dxs[0])) )
  print('convergence w', (log(errsw[-1])-log(errsw[0])) / (log(dxs[-1])-log(dxs[0])) )
print('convergence p', (log(errsp[-1])-log(errsp[0])) / (log(dxs[-1])-log(dxs[0])) )

off2 = 0.01
off4 = 0.002
slope2 = off2*(dxs[:] / dxs[0])**2.
slope4 = off4*(dxs[:] / dxs[0])**4.

close('all')

figure()
if(t > 0):
  loglog(dxs, errsu, 'bo-', label="u")
  loglog(dxs, errsw, 'go-', label="w")
loglog(dxs, errsp, 'ro-', label="p")
loglog(dxs, slope2, 'k--', label="2nd")
loglog(dxs, slope4, 'k:' , label="4th")
legend(loc=4, frameon=False)
#xlim(0.01, 0.2)

figure()
subplot(121)
pcolormesh(data32.xh, data32.z, data32.u[:,0,:])
xlim(min(data32.xh), max(data32.xh))
ylim(min(data32.z) , max(data32.z ))
xlabel('x')
ylabel('z')
title('u')
colorbar()
subplot(122)
pcolormesh(data32.xh, data32.z, ref32.u[:,0,:])
xlim(min(data32.xh), max(data32.xh))
ylim(min(data32.z ), max(data32.z ))
xlabel('x')
ylabel('z')
title('u ref')
colorbar()

figure()
subplot(121)
pcolormesh(data32.x, data32.zh, data32.w[:,0,:])
xlim(min(data32.x ), max(data32.x ))
ylim(min(data32.zh), max(data32.zh))
xlabel('x')
ylabel('z')
title('w')
colorbar()
subplot(122)
pcolormesh(data32.x, data32.zh, ref32.w[:,0,:])
xlim(min(data32.x ), max(data32.x ))
ylim(min(data32.zh), max(data32.zh))
xlabel('x')
ylabel('z')
title('w ref')
colorbar()

figure()
subplot(121)
pcolormesh(data32.x, data32.z, data32.p[:,0,:])
xlim(min(data32.x), max(data32.x))
ylim(min(data32.z), max(data32.z))
xlabel('x')
ylabel('z')
title('p')
colorbar()
subplot(122)
pcolormesh(data32.x, data32.z, ref32.p[:,0,:])
xlim(min(data32.x), max(data32.x))
ylim(min(data32.z), max(data32.z))
xlabel('x')
ylabel('z')
title('p ref')
colorbar()

figure()
pcolormesh(data32.x, data32.z, data32.p[:,0,:]-ref32.p[:,0,:])
xlim(min(data32.x), max(data32.x))
ylim(min(data32.z), max(data32.z))
xlabel('x')
ylabel('z')
title('u')
title('p err')
colorbar()

