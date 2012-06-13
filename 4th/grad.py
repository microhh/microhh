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

def gx2nd(x, u):
  istart = igc
  iend   = u.size + igc
  dx     = xsize / u.size
  gu     = zeros(u.size)
  ucalc  = zeros(u.size + 2*igc)
  ucalc[istart:iend] = u[:]

  # periodic bcs
  ucalc[0   :igc     ] = u[u.size-igc:u.size]
  ucalc[iend:iend+igc] = u[0:igc]

  for i in range(istart, iend):
    gu[i-igc] = (ucalc[i+1] - ucalc[i-1]) / (2.*dx)

  guref = 2*pi/xsize * cos(2.*pi*x/xsize)

  erru = geterror(gu, guref)

  return gu, guref, erru

def gx4th(x, u):
  istart = igc
  iend   = u.size + igc
  dx     = xsize / u.size
  gu     = zeros(u.size)
  ucalc  = zeros(u.size + 2*igc)
  ucalc[istart:iend] = u[:]

  # periodic bcs
  ucalc[0   :igc     ] = u[u.size-igc:u.size]
  ucalc[iend:iend+igc] = u[0:igc]

  for i in range(istart, iend):
    gu[i-igc] = ((ucalc[i+3]-ucalc[i-3]) - 36.*(ucalc[i+2]-ucalc[i-2]) + 261.*(ucalc[i+1]-ucalc[i-1])) / (384.*dx)

  guref = 2*pi/xsize * cos(2.*pi*x/xsize)

  erru = geterror(gu, guref)

  return gu, guref, erru


x8, u8 = refdata(8)
dx8    = xsize / 8.
g8_2nd, gref8_2nd, err8_2nd = gx2nd(x8, u8)
g8_4th, gref8_4th, err8_4th = gx4th(x8, u8)

x16, u16 = refdata(16)
dx16     = xsize / 16.
g16_2nd, gref16_2nd, err16_2nd = gx2nd(x16, u16)
g16_4th, gref16_4th, err16_4th = gx4th(x16, u16)

x32, u32 = refdata(32)
dx32     = xsize / 32.
g32_2nd, gref32_2nd, err32_2nd = gx2nd(x32, u32)
g32_4th, gref32_4th, err32_4th = gx4th(x32, u32)

x64, u64 = refdata(64)
dx64     = xsize / 64.
g64_2nd, gref64_2nd, err64_2nd = gx2nd(x64, u64)
g64_4th, gref64_4th, err64_4th = gx4th(x64, u64)

dxs  = array([dx8, dx16, dx32, dx64])
errs_2nd = array([err8_2nd, err16_2nd, err32_2nd, err64_2nd])
errs_4th = array([err8_4th, err16_4th, err32_4th, err64_4th])

off2 = 0.2
off4 = 0.02
slope2 = off2*(dxs[:] / dxs[0])**2.
slope4 = off4*(dxs[:] / dxs[0])**4.

close('all')
figure()
plot(x8 , g8_2nd , 'b-',  label="8_2nd" )
plot(x16, g16_2nd, 'g-',  label="16_2nd")
plot(x32, g32_2nd, 'r-',  label="32_2nd")
plot(x8 , g8_4th , 'b--', label="8_4th" )
plot(x16, g16_4th, 'g--', label="16_4th")
plot(x32, g32_4th, 'r--', label="32_4th")
legend(loc=4, frameon=False)

figure()
loglog(dxs, errs_2nd, 'bo-', label="g2nd")
loglog(dxs, errs_4th, 'go-', label="g4th")
loglog(dxs, slope2, 'k--', label="2nd")
loglog(dxs, slope4, 'k:' , label="4th")
legend(loc=4, frameon=False)
xlim(0.01, 0.2)

