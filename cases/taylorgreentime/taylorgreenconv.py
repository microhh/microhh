from pylab import *
from taylorgreenfunc import *

time = 1.
visc = (8.*pi**2. * 100.)**(-1.)

nts   = array([100, 200, 400])
dts   = 1./nts

# 2nd order data
data100  = microhh(100, 128, 64, 'taylorgreen100' )
data200  = microhh(200, 128, 64, 'taylorgreen200' )
data400  = microhh(400, 128, 64, 'taylorgreen400' )

ref100  = getref(data100.x, data100.xh, data100.z, data100.zh, visc, time)
ref200  = getref(data200.x, data200.xh, data200.z, data200.zh, visc, time)
ref400  = getref(data400.x, data400.xh, data400.z, data400.zh, visc, time)

err100  = geterror(data100, ref100)
err200  = geterror(data200, ref200)
err400  = geterror(data400, ref400)

errsu = array([err100.u, err200.u, err400.u])
errsw = array([err100.w, err200.w, err400.w])
errsp = array([err100.p, err200.p, err400.p])

print('errors u', errsu)
print('errors w', errsw)
print('errors p', errsp)
print('convergence u', (log(errsu[-1])-log(errsu[0])) / (log(dts[-1])-log(dts[0])) )
print('convergence w', (log(errsw[-1])-log(errsw[0])) / (log(dts[-1])-log(dts[0])) )
print('convergence p', (log(errsp[-1])-log(errsp[0])) / (log(dts[-1])-log(dts[0])) )
off2 = 0.01
off4 = 0.002
slope2 = off2*(dts[:] / dts[0])**2.
slope4 = off4*(dts[:] / dts[0])**4.

close('all')

figure()
loglog(dts, errsu, 'b-', label="u")
loglog(dts, errsw, 'g-', label="w")
loglog(dts, errsp, 'r-', label="p")
loglog(dts, slope2, 'k--', label="2nd")
loglog(dts, slope4, 'k:' , label="4th")
legend(loc=0, frameon=False)

