from pylab import *
from taylorgreenfunc import *

t    = 400
time = 1.
visc = (8.*pi**2. * 100.)**(-1.)

ns    = array([32, 64, 128, 256])
dxs   = 1./ns

# 2nd order data
data32_2nd  = microhh(t,  32,  16, 'taylorgreen32_2nd' )
data64_2nd  = microhh(t,  64,  32, 'taylorgreen64_2nd' )
data128_2nd = microhh(t, 128,  64, 'taylorgreen128_2nd')
data256_2nd = microhh(t, 256, 128, 'taylorgreen256_2nd')

ref32_2nd  = getref(data32_2nd .x, data32_2nd .xh, data32_2nd .z, data32_2nd .zh, visc, time)
ref64_2nd  = getref(data64_2nd .x, data64_2nd .xh, data64_2nd .z, data64_2nd .zh, visc, time)
ref128_2nd = getref(data128_2nd.x, data128_2nd.xh, data128_2nd.z, data128_2nd.zh, visc, time)
ref256_2nd = getref(data256_2nd.x, data256_2nd.xh, data256_2nd.z, data256_2nd.zh, visc, time)

err32_2nd  = geterror(data32_2nd , ref32_2nd )
err64_2nd  = geterror(data64_2nd , ref64_2nd )
err128_2nd = geterror(data128_2nd, ref128_2nd)
err256_2nd = geterror(data256_2nd, ref256_2nd)

errsu_2nd = array([err32_2nd.u, err64_2nd.u, err128_2nd.u, err256_2nd.u])
errsw_2nd = array([err32_2nd.w, err64_2nd.w, err128_2nd.w, err256_2nd.w])
errsp_2nd = array([err32_2nd.p, err64_2nd.p, err128_2nd.p, err256_2nd.p])

print('errors p_2nd', errsp_2nd)
if(t > 0):
  print('convergence u_2nd', (log(errsu_2nd[-1])-log(errsu_2nd[0])) / (log(dxs[-1])-log(dxs[0])) )
  print('convergence w_2nd', (log(errsw_2nd[-1])-log(errsw_2nd[0])) / (log(dxs[-1])-log(dxs[0])) )
print('convergence p_2nd', (log(errsp_2nd[-1])-log(errsp_2nd[0])) / (log(dxs[-1])-log(dxs[0])) )

# 42 order data
data32_42  = microhh(t,  32,  16, 'taylorgreen32_42' )
data64_42  = microhh(t,  64,  32, 'taylorgreen64_42' )
data128_42 = microhh(t, 128,  64, 'taylorgreen128_42')
data256_42 = microhh(t, 256, 128, 'taylorgreen256_42')

ref32_42  = getref(data32_42 .x, data32_42 .xh, data32_42 .z, data32_42 .zh, visc, time)
ref64_42  = getref(data64_42 .x, data64_42 .xh, data64_42 .z, data64_42 .zh, visc, time)
ref128_42 = getref(data128_42.x, data128_42.xh, data128_42.z, data128_42.zh, visc, time)
ref256_42 = getref(data256_42.x, data256_42.xh, data256_42.z, data256_42.zh, visc, time)

err32_42  = geterror(data32_42 , ref32_42 )
err64_42  = geterror(data64_42 , ref64_42 )
err128_42 = geterror(data128_42, ref128_42)
err256_42 = geterror(data256_42, ref256_42)

errsu_42 = array([err32_42.u, err64_42.u, err128_42.u, err256_42.u])
errsw_42 = array([err32_42.w, err64_42.w, err128_42.w, err256_42.w])
errsp_42 = array([err32_42.p, err64_42.p, err128_42.p, err256_42.p])

print('errors p_42', errsp_42)
if(t > 0):
  print('convergence u_42', (log(errsu_42[-1])-log(errsu_42[0])) / (log(dxs[-1])-log(dxs[0])) )
  print('convergence w_42', (log(errsw_42[-1])-log(errsw_42[0])) / (log(dxs[-1])-log(dxs[0])) )
print('convergence p_42', (log(errsp_42[-1])-log(errsp_42[0])) / (log(dxs[-1])-log(dxs[0])) )

off2 = 0.01
off4 = 0.002
slope2 = off2*(dxs[:] / dxs[0])**2.
slope4 = off4*(dxs[:] / dxs[0])**4.

close('all')

figure()
if(t > 0):
  loglog(dxs, errsu_2nd, 'bo-', label="u_2nd")
  loglog(dxs, errsw_2nd, 'bv-', label="w_2nd")
  loglog(dxs, errsu_42 , 'go-', label="u_42" )
  loglog(dxs, errsw_42 , 'gv-', label="w_42" )
loglog(dxs, errsp_2nd, 'b^-', label="p_2nd")
loglog(dxs, errsp_42 , 'g^-', label="p_42" )
loglog(dxs, slope2, 'k--', label="2nd")
loglog(dxs, slope4, 'k:' , label="4th")
legend(loc=4, frameon=False)
#xlim(0.01, 0.2)

"""
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
"""
