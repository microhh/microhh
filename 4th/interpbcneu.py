from pylab import *

xsize = 1.
igc   = 3

def geterror(u, uref, n=2):
  error = sqrt(sum((uref - u)**2.) / u.size)
  return error

def refdata(n):
  x  = linspace(0.5*xsize/n, xsize-0.5*xsize/n, n)
  xi = linspace(0., xsize, n+1)[0:n]
  u = cos(2.*pi*x/xsize)
  return x, xi, u


def ix2nd(x, u):
  istart = igc
  iend   = u.size + igc
  dx     = xsize / u.size
  iu     = zeros(u.size)
  ucalc  = zeros(u.size + 2*igc)
  ucalc[istart:iend] = u[:]

  # periodic bcs
  # ucalc[0   :igc     ] = u[u.size-igc:u.size]
  # ucalc[iend:iend+igc] = u[0:igc]

  # non-periodic bc
  ucalc[istart-1] = ucalc[istart]
  ucalc[iend    ] = ucalc[iend-1]

  for i in range(istart, iend):
    iu[i-igc] = (ucalc[i] + ucalc[i-1]) / 2.

  iuref = cos(2.*pi*x/xsize)

  erru = geterror(iu, iuref)

  return iu, iuref, erru

def ix4th(x, u):
  istart = igc
  iend   = u.size + igc
  dx     = xsize / u.size
  iu     = zeros(u.size)
  ucalc  = zeros(u.size + 2*igc)
  ucalc[istart:iend] = u[:]

  # periodic bcs
  # ucalc[0   :igc     ] = u[u.size-igc:u.size]
  # ucalc[iend:iend+igc] = u[0:igc]

  # ghost cell
  ucalc[istart-1] = (21.*ucalc[istart] + 3.*ucalc[istart+1] - ucalc[istart+2]) / 23.
  ucalc[iend    ] = (21.*ucalc[iend-1] + 3.*ucalc[iend-2  ] - ucalc[iend-3  ]) / 23.

  i = istart
  iu[i-igc] = (5.*ucalc[i-1] + 15.*ucalc[i] - 5.*ucalc[i+1] + ucalc[i+2]) / 16.
  for i in range(istart+1, iend):
    iu[i-igc] = (- ucalc[i-2] + 9.*ucalc[i-1] + 9.*ucalc[i] - ucalc[i+1]) / 16.

  iuref = cos(2.*pi*x/xsize)

  erru = geterror(iu, iuref)

  return iu, iuref, erru


x8, xi8, u8 = refdata(8)
dx8         = xsize / 8.
i8_2nd, iref8_2nd, err8_2nd = ix2nd(xi8, u8)
i8_4th, iref8_4th, err8_4th = ix4th(xi8, u8)

x16, xi16, u16 = refdata(16)
dx16           = xsize / 16.
i16_2nd, iref16_2nd, err16_2nd = ix2nd(xi16, u16)
i16_4th, iref16_4th, err16_4th = ix4th(xi16, u16)

x32, xi32, u32 = refdata(32)
dx32           = xsize / 32.
i32_2nd, iref32_2nd, err32_2nd = ix2nd(xi32, u32)
i32_4th, iref32_4th, err32_4th = ix4th(xi32, u32)

x64, xi64, u64 = refdata(64)
dx64           = xsize / 64.
i64_2nd, iref64_2nd, err64_2nd = ix2nd(xi64, u64)
i64_4th, iref64_4th, err64_4th = ix4th(xi64, u64)

dxs  = array([dx8, dx16, dx32, dx64])
errs_2nd = array([err8_2nd, err16_2nd, err32_2nd, err64_2nd])
errs_4th = array([err8_4th, err16_4th, err32_4th, err64_4th])

print('convergence 2nd', (log(errs_2nd[-1])-log(errs_2nd[0])) / (log(dxs[-1])-log(dxs[0])) )
print('convergence 4th', (log(errs_4th[-1])-log(errs_4th[0])) / (log(dxs[-1])-log(dxs[0])) )

off2 = 0.2
off4 = 0.02
slope2 = off2*(dxs[:] / dxs[0])**2.
slope4 = off4*(dxs[:] / dxs[0])**4.

close('all')
figure()
plot(xi8 , i8_2nd , 'b-',  label="8_2nd" )
plot(xi16, i16_2nd, 'g-',  label="16_2nd")
plot(xi32, i32_2nd, 'r-',  label="32_2nd")
plot(xi8 , i8_4th , 'b--', label="8_4th" )
plot(xi16, i16_4th, 'g--', label="16_4th")
plot(xi32, i32_4th, 'r--', label="32_4th")
legend(loc=4, frameon=False)

figure()
loglog(dxs, errs_2nd, 'bo-', label="g2nd")
loglog(dxs, errs_4th, 'go-', label="g4th")
loglog(dxs, slope2, 'k--', label="2nd")
loglog(dxs, slope4, 'k:' , label="4th")
legend(loc=4, frameon=False)
xlim(0.01, 0.2)

