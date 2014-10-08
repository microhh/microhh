import numpy

# set the height
kmax  = 32
zsize = 4400.
dz = zsize / kmax

# set the height
z    = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl  = numpy.zeros(numpy.size(z))
qt   = numpy.zeros(numpy.size(z))
u    = numpy.zeros(numpy.size(z))
ug   = numpy.zeros(numpy.size(z))

for k in range(kmax):
  # temperature
  if(z[k] <= 50.):
    thl[k] = 299.0  + (z[k]     )*(301.5 -299.0 )/(50.)
    qt[k] = 15.20  + (z[k]     )*(15.17 -15.20 )/(50.)
  elif(z[k] <=  350.):
    thl[k] = 301.5  + (z[k]-  50.)*(302.5 -301.5 )/(350.-50.)
    qt[k] = 15.17  + (z[k]-  50.)*(14.98 -15.17 )/(350.-50.)
  elif(z[k] <=  650.):
    thl[k] = 302.5  + (z[k]- 350.)*(303.53-302.5 )/(650.-350.)
    qt[k] = 14.98  + (z[k]- 350.)*(14.80 -14.98 )/(650.-350.)
  elif(z[k] <=  700.):
    thl[k] = 303.53 + (z[k]- 650.)*(303.7 -303.53)/(700.-650.)
    qt[k] = 14.80  + (z[k]- 650.)*(14.70 -14.80 )/(700.-650.)
  elif(z[k] <= 1300.):
    thl[k] = 303.7  + (z[k]- 700.)*(307.13-303.7 )/(1300.-700.)
    qt[k] = 14.70  + (z[k]- 700.)*( 13.50-14.80 )/(1300.-700.)
  elif(z[k] <= 2500.):
    thl[k] = 307.13 + (z[k]-1300.)*(314.0 -307.13)/(2500.-1300.)
    qt[k] = 13.50  + (z[k]-1300.)*( 3.00 - 13.50)/(2500.-1300.)
  elif(z[k] <= 5500.):
    thl[k] = 314.0  + (z[k]-2500.)*(343.2 -314.0 )/(5500.-2500.)
    qt[k] =  3.00

  # u-wind component
  u[:] = 10.

  # ug-wind component
  ug[k] = 10.

# normalize profiles to SI
qt /= 1000.  # g to kg

# set the time series
t  = numpy.array([  0.,   4.,  6.5,  7.5,  10., 12.5, 14.5])

#H  = numpy.array([-30.,  90., 140., 140., 100., -10.,  -10]) # Original 
H  = numpy.array([  0.,  90., 140., 140., 100., -10.,    0]) # No negative fluxes
print('Negative sensible heat flux set to zero to prevent problems with Obuk. iteration!!')

LE = numpy.array([  5., 250., 450., 500., 420., 180.,    0])

advthl = numpy.array([ 0.   , 0.  ,  0.  , -0.08, -0.16, -0.16])
radthl = numpy.array([-0.125, 0.  ,  0.  ,  0.  ,  0.   , -0.1])
advqt  = numpy.array([ 0.08 , 0.02, -0.04, -0.10, -0.16, -0.30])

tls    = numpy.array([  0.,   3.,  6.,  9.,  12., 14.5])
thlls  = numpy.zeros((t.size, kmax))
qtls   = numpy.zeros((t.size, kmax))

# large scale forcings
for n in range(tls.size):
  tendthl = advthl[n] + radthl[n]
  tendqt  = advqt[n]
  for k in range(kmax):
    # temperature
    if(z[k] <= 1000.):
      thlls[n,k] = tendthl
      qtls[n,k] = tendqt
    else:
      thlls[n,k] = tendthl - (z[k]-1000.)*(tendthl)/(5500.-1000.)
      qtls[n,k] = tendqt - (z[k]-1000.)*(tendqt)/(5500.-1000.)
tls   *= 3600. # h to s
thlls /= 3600. # h to s
qtls  /= 3600. # h to s
qtls  /= 1000. # g to kg

# write the prof data to a file
proffile = open('arm.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s}\n'.format('z','thl','qt','u','ug'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E}\n'.format(z[k], thl[k], qt[k], u[k], ug[k]))
proffile.close()

# write the time data to a file
Rd  = 287.
cp  = 1005.
Lv  = 2.5e6
p0  = 97000.
rho = p0/(Rd*thl[0]*(1. + 0.61*qt[0]))
print('rho = ', rho)
t *= 3600. # h to s
sbotthl = H/(rho*cp)
sbotqt  = LE/(rho*Lv)
timefile = open('arm.time','w')
timefile.write('{0:^20s} {1:^20s} {2:^20s}\n'.format('t','sbot[thl]', 'sbot[qt]'))
for n in range(t.size):
  timefile.write('{0:1.14E} {1:1.14E} {2:1.14E}\n'.format(t[n], sbotthl[n], sbotqt[n]))
timefile.close()

# write the large scale forcing data to a file
timeproffile = open('thlls.timeprof','w')
timeproffile.write('{0:^20s} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E} {6:1.14E}\n'.format('z', tls[0], tls[1], tls[2], tls[3], tls[4], tls[5]))
for k in range(kmax):
  timeproffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E} {6:1.14E}\n'.format(z[k], thlls[0,k], thlls[1,k], thlls[2,k], thlls[3,k], thlls[4,k], thlls[5,k]))
timeproffile.close()

timeproffile = open('qtls.timeprof','w')
timeproffile.write('{0:^20s} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E} {6:1.14E}\n'.format('z', tls[0], tls[1], tls[2], tls[3], tls[4], tls[5]))
for k in range(kmax):
  timeproffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E} {6:1.14E}\n'.format(z[k], qtls[0,k], qtls[1,k], qtls[2,k], qtls[3,k], qtls[4,k], qtls[5,k]))
timeproffile.close()

