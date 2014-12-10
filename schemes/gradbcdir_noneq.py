from pylab import *

xsize = 0.5
igc = 3
stretch = 2.75

def geterror(u, uref, n=2):
  error = sqrt(sum((uref - u)**2.) / u.size)
  #error = abs(max(uref - u))
  return error

def refdata(n):
  a = linspace(1./(n+1.), 1.-1./(n+1.), n)
  #x = 0.5 * (1. - stretch**a) / (1. - stretch)
  x = 0.25 + 0.25 * tanh(stretch * (2.*a - 1.) ) / tanh( stretch ) 
  u = sin(2.*pi*x/xsize)
  return x, u

def gx2nd(x, u):
  istart = igc
  iend   = u.size + igc
  dx     = xsize / u.size
  gu     = zeros(u.size)

  ucalc  = zeros(u.size + 2*igc)
  xcalc  = zeros(x.size + 2*igc)

  xcalc[istart:iend] = x[:]
  ucalc[istart:iend] = u[:]

  # # periodic bcs
  # ucalc[0   :igc     ] = u[u.size-igc:u.size]
  # ucalc[iend:iend+igc] = u[0:igc]

  # non-periodic bc
  xcalc[istart-1] = - xcalc[istart]
  xcalc[iend    ] = 2.*xsize - xcalc[iend-1]

  ucalc[istart-1] = - ucalc[istart]
  ucalc[iend    ] = - ucalc[iend-1]

  for i in range(istart, iend):
    gu[i-igc] = (ucalc[i+1] - ucalc[i-1]) / (xcalc[i+1] - xcalc[i-1])

  guref = 2.*pi/xsize * cos(2.*pi*x/xsize)

  erru = geterror(gu, guref)

  return gu, guref, erru

def gx4th(x, u):
  istart = igc
  iend   = u.size + igc
  dx     = xsize / u.size
  gu     = zeros(u.size)

  ucalc  = zeros(u.size + 2*igc)
  xcalc  = zeros(x.size + 2*igc)

  xcalc[istart:iend] = x[:]
  ucalc[istart:iend] = u[:]

  # periodic bcs
  #ucalc[0   :igc     ] = u[u.size-igc:u.size]
  #ucalc[iend:iend+igc] = u[0:igc]

  # ghost cell
  xcalc[istart-1] = -2.*xcalc[istart] + (1./3.)*xcalc[istart+1]
  xcalc[istart-2] = -9.*xcalc[istart] + 2.*xcalc[istart+1]
  #xcalc[istart-1] = -xcalc[istart  ]
  #xcalc[istart-2] = -xcalc[istart+1]

  xcalc[iend  ] = (8./3.)*xsize - 2.*xcalc[iend-1] + (1./3.)*xcalc[iend-2]
  xcalc[iend+1] = 8.*xsize - 9.*xcalc[iend-1] + 2.*xcalc[iend-2]
  #xcalc[iend  ] = 2.*xsize - xcalc[iend-1]
  #xcalc[iend+1] = 2.*xsize - xcalc[iend-2]

  ucalc[istart-1] = -2.*ucalc[istart] + (1./3.)*ucalc[istart+1]
  ucalc[istart-2] = -9.*ucalc[istart] + 2.*ucalc[istart+1]

  ucalc[iend  ] = -2.*ucalc[iend-1] + (1./3.)*ucalc[iend-2]
  ucalc[iend+1] = -9.*ucalc[iend-1] + 2.*ucalc[iend-2]

  ci0 = -1./16.; ci1 =   9./16.; ci2 =   9./16.; ci3 = -1./16.;
  bi0 =  5./16.; bi1 =  15./16.; bi2 =  -5./16.; bi3 =  1./16.;
  ti0 =  1./16.; ti1 =  -5./16.; ti2 =  15./16.; ti3 =  5./16.;
  cg0 =  1./24.; cg1 = -27./24.; cg2 =  27./24.; cg3 = -1./24.;

  i = istart
  gu[i-igc] = ( cg0*(bi0*ucalc[i-2] + bi1*ucalc[i-1] + bi2*ucalc[i  ] + bi3*ucalc[i+1])
              + cg1*(ci0*ucalc[i-2] + ci1*ucalc[i-1] + ci2*ucalc[i  ] + ci3*ucalc[i+1])
              + cg2*(ci0*ucalc[i-1] + ci1*ucalc[i  ] + ci2*ucalc[i+1] + ci3*ucalc[i+2])
              + cg3*(ci0*ucalc[i  ] + ci1*ucalc[i+1] + ci2*ucalc[i+2] + ci3*ucalc[i+3]) ) \
            / ( cg0*(bi0*xcalc[i-2] + bi1*xcalc[i-1] + bi2*xcalc[i  ] + bi3*xcalc[i+1])
              + cg1*(ci0*xcalc[i-2] + ci1*xcalc[i-1] + ci2*xcalc[i  ] + ci3*xcalc[i+1])
              + cg2*(ci0*xcalc[i-1] + ci1*xcalc[i  ] + ci2*xcalc[i+1] + ci3*xcalc[i+2])
              + cg3*(ci0*xcalc[i  ] + ci1*xcalc[i+1] + ci2*xcalc[i+2] + ci3*xcalc[i+3]) )
  for i in range(istart+1, iend-1):
    gu[i-igc] = ( cg0*(ci0*ucalc[i-3] + ci1*ucalc[i-2] + ci2*ucalc[i-1] + ci3*ucalc[i  ])
                + cg1*(ci0*ucalc[i-2] + ci1*ucalc[i-1] + ci2*ucalc[i  ] + ci3*ucalc[i+1])
                + cg2*(ci0*ucalc[i-1] + ci1*ucalc[i  ] + ci2*ucalc[i+1] + ci3*ucalc[i+2])
                + cg3*(ci0*ucalc[i  ] + ci1*ucalc[i+1] + ci2*ucalc[i+2] + ci3*ucalc[i+3]) ) \
              / ( cg0*(ci0*xcalc[i-3] + ci1*xcalc[i-2] + ci2*xcalc[i-1] + ci3*xcalc[i  ])
                + cg1*(ci0*xcalc[i-2] + ci1*xcalc[i-1] + ci2*xcalc[i  ] + ci3*xcalc[i+1])
                + cg2*(ci0*xcalc[i-1] + ci1*xcalc[i  ] + ci2*xcalc[i+1] + ci3*xcalc[i+2])
                + cg3*(ci0*xcalc[i  ] + ci1*xcalc[i+1] + ci2*xcalc[i+2] + ci3*xcalc[i+3]) )
  i = iend-1
  gu[i-igc] = ( cg0*(ci0*ucalc[i-3] + ci1*ucalc[i-2] + ci2*ucalc[i-1] + ci3*ucalc[i  ])
              + cg1*(ci0*ucalc[i-2] + ci1*ucalc[i-1] + ci2*ucalc[i  ] + ci3*ucalc[i+1])
              + cg2*(ci0*ucalc[i-1] + ci1*ucalc[i  ] + ci2*ucalc[i+1] + ci3*ucalc[i+2])
              + cg3*(ti0*ucalc[i-1] + ti1*ucalc[i  ] + ti2*ucalc[i+1] + ti3*ucalc[i+2]) ) \
            / ( cg0*(ci0*xcalc[i-3] + ci1*xcalc[i-2] + ci2*xcalc[i-1] + ci3*xcalc[i  ])
              + cg1*(ci0*xcalc[i-2] + ci1*xcalc[i-1] + ci2*xcalc[i  ] + ci3*xcalc[i+1])
              + cg2*(ci0*xcalc[i-1] + ci1*xcalc[i  ] + ci2*xcalc[i+1] + ci3*xcalc[i+2])
              + cg3*(ti0*xcalc[i-1] + ti1*xcalc[i  ] + ti2*xcalc[i+1] + ti3*xcalc[i+2]) )

  """
  cg0 = -1./12.; cg1 = 8./12.; cg2 = -8./12.; cg3 = 1./12.;
  for i in range(istart, iend):
    gu[i-igc] = (cg0*ucalc[i-2] + cg1*ucalc[i-1] + cg2*ucalc[i+1] + cg3*ucalc[i+2]) \
              / (cg0*xcalc[i-2] + cg1*xcalc[i-1] + cg2*xcalc[i+1] + cg3*xcalc[i+2])
              """

  guref = 2.*pi/xsize * cos(2.*pi*x/xsize)

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

x128, u128 = refdata(128)
dx128      = xsize / 128.
g128_2nd, gref128_2nd, err128_2nd = gx2nd(x128, u128)
g128_4th, gref128_4th, err128_4th = gx4th(x128, u128)

dxs  = array([dx8, dx16, dx32, dx64, dx128])
errs_2nd = array([err8_2nd, err16_2nd, err32_2nd, err64_2nd, err128_2nd])
errs_4th = array([err8_4th, err16_4th, err32_4th, err64_4th, err128_4th])

print('convergence 2nd', (log(errs_2nd[-1])-log(errs_2nd[0])) / (log(dxs[-1])-log(dxs[0])) )
print('convergence 4th', (log(errs_4th[-1])-log(errs_4th[0])) / (log(dxs[-1])-log(dxs[0])) )

off2 = 0.4
off3 = 0.3
off4 = 0.1
slope2 = off2*(dxs[:] / dxs[0])**2.
slope3 = off3*(dxs[:] / dxs[0])**3.
slope4 = off4*(dxs[:] / dxs[0])**4.

close('all')
figure()
plot(x8 , g8_2nd  - gref8_2nd ,'b-',  label="8_2nd" )
plot(x16, g16_2nd - gref16_2nd,'g-',  label="16_2nd")
plot(x32, g32_2nd - gref32_2nd,'r-',  label="32_2nd")
plot(x64, g64_2nd - gref64_2nd,'c-',  label="64_2nd")
plot(x8 , g8_4th  - gref8_4th ,'b--', label="8_4th" )
plot(x16, g16_4th - gref16_4th,'g--', label="16_4th")
plot(x32, g32_4th - gref32_4th,'r--', label="32_4th")
plot(x64, g64_4th - gref64_4th,'c--', label="64_4th")
#plot(x128, gref128_4th, 'k:', label="ref")
legend(loc=4, frameon=False)

figure()
loglog(dxs, errs_2nd, 'bo-', label="g2nd")
loglog(dxs, errs_4th, 'go-', label="g4th")
loglog(dxs, slope2, 'k--', label="2nd")
loglog(dxs, slope3, 'k-.', label="3rd")
loglog(dxs, slope4, 'k:' , label="4th")
legend(loc=4, frameon=False)
#xlim(0.01, 0.2)

