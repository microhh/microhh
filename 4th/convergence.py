from pylab import *

xsize = 1.
igc   = 3

def refdata(n):
  x = linspace(0., xsize, n+1)[0:n]
  u = sin(2.*pi*x/xsize)
  return x, u

def dgx4th(x, u):
  istart = igc
  iend   = u.size + igc
  dx     = xsize / u.size
  dgu    = zeros(u.size)
  ucalc  = zeros(u.size + 2*igc)
  ucalc[istart:iend] = u[:]

  # periodic bcs
  ucalc[0   :igc     ] = u[u.size-igc:u.size]
  ucalc[iend:iend+igc] = u[0:igc]

  for i in range(istart, iend):
    dgu[i-igc] = ((ucalc[i-3]+ucalc[i+3]) -54.*(ucalc[i-2]+ucalc[i+2]) + 783.*(ucalc[i-1]+ucalc[i+1]) - 1460.*ucalc[i]) / (576.*dx**2.)

  dguref = -4*pi**2./xsize**2. * sin(2.*pi*x/xsize)

  return dgu, dguref

x1024, u1024      = refdata(1024)
dg1024, dgref1024 = dgx4th(x1024, u1024)

close('all')
figure()
plot(x1024, dg1024   , label="1024")
plot(x1024, dgref1024, label="ref" )
legend(loc=0, frameon=False)

