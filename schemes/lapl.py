from pylab import *

xsize = 1.
igc   = 3

def geterror(u, uref, n=2):
  error = sqrt(sum((uref - u)**2.) / u.size)
  return error

def refdata(n):
  x   = linspace(0., xsize, n+1)[0:n]
  u   = sin(2.*pi*x/xsize)
  return x, u

def laplx2nd(x, u):
  istart = igc
  iend   = u.size + igc
  dx     = xsize / u.size
  laplu    = zeros(u.size)
  ucalc  = zeros(u.size + 2*igc)
  ucalc[istart:iend] = u[:]

  # periodic bcs
  ucalc[0   :igc     ] = u[u.size-igc:u.size]
  ucalc[iend:iend+igc] = u[0:igc]

  for i in range(istart, iend):
    laplu[i-igc] = (ucalc[i-1] - 2.*ucalc[i] + ucalc[i+1]) / (dx**2.)

  lapluref = -4*pi**2./xsize**2. * sin(2.*pi*x/xsize)

  erru = geterror(laplu, lapluref)

  return laplu, lapluref, erru


def laplx4th(x, u):
  istart = igc
  iend   = u.size + igc
  dx     = xsize / u.size
  laplu    = zeros(u.size)
  ucalc  = zeros(u.size + 2*igc)
  ucalc[istart:iend] = u[:]

  # periodic bcs
  ucalc[0   :igc     ] = u[u.size-igc:u.size]
  ucalc[iend:iend+igc] = u[0:igc]

  for i in range(istart, iend):
    laplu[i-igc] = (-1.*(ucalc[i-2]+ucalc[i+2]) + 16.*(ucalc[i-1]+ucalc[i+1]) - 30.*ucalc[i]) / (12.*dx**2.)

  lapluref = -4*pi**2./xsize**2. * sin(2.*pi*x/xsize)

  erru = geterror(laplu, lapluref)

  return laplu, lapluref, erru

x8, u8 = refdata(8)
dx8    = xsize / 8.
lapl8_2nd, laplref8_2nd, err8_2nd = laplx2nd(x8, u8)
lapl8_4th, laplref8_4th, err8_4th = laplx4th(x8, u8)

x16, u16 = refdata(16)
dx16     = xsize / 16.
lapl16_2nd, laplref16_2nd, err16_2nd = laplx2nd(x16, u16)
lapl16_4th, laplref16_4th, err16_4th = laplx4th(x16, u16)

x32, u32 = refdata(32)
dx32     = xsize / 32.
lapl32_2nd, laplref32_2nd, err32_2nd = laplx2nd(x32, u32)
lapl32_4th, laplref32_4th, err32_4th = laplx4th(x32, u32)

x64, u64 = refdata(64)
dx64     = xsize / 64.
lapl64_2nd, laplref64_2nd, err64_2nd = laplx2nd(x64, u64)
lapl64_4th, laplref64_4th, err64_4th = laplx4th(x64, u64)

dxs  = array([dx8 , dx16 , dx32 , dx64 ])
errs_2nd = array([err8_2nd, err16_2nd, err32_2nd, err64_2nd])
errs_4th = array([err8_4th, err16_4th, err32_4th, err64_4th])

print('convergence 2nd', (log(errs_2nd[-1])-log(errs_2nd[0])) / (log(dxs[-1])-log(dxs[0])) )
print('convergence 4th', (log(errs_4th[-1])-log(errs_4th[0])) / (log(dxs[-1])-log(dxs[0])) )

off2 = 0.7
off4 = 0.04
slope2 = off2*(dxs[:] / dxs[0])**2.
slope4 = off4*(dxs[:] / dxs[0])**4.

close('all')
figure()
plot(x8 , lapl8_2nd , 'b-',  label="8_2nd" )
plot(x16, lapl16_2nd, 'g-',  label="16_2nd")
plot(x32, lapl32_2nd, 'r-',  label="32_2nd")
plot(x8 , lapl8_4th , 'b--', label="8_4th" )
plot(x16, lapl16_4th, 'g--', label="16_4th")
plot(x32, lapl32_4th, 'r--', label="32_4th")
legend(loc=0, frameon=False)

figure()
loglog(dxs, errs_2nd, 'bo-', label="lapl2nd")
loglog(dxs, errs_4th, 'go-', label="lapl4th")
loglog(dxs, slope2, 'k--', label="2nd")
loglog(dxs, slope4, 'k:' , label="4th")
legend(loc=4, frameon=False)
xlim(0.01, 0.2)

