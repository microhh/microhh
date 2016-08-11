import numpy as np
import matplotlib.pylab as pl

pl.close('all')

# Get number of vertical levels and size from .ini file
with open('rain.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

# Well-mixed layer
zi1     = 1000     # m
thl0    = 293      # K
qt0     = 8e-3    # kg/kg

# Cloud layer:
dthl1   = 0        # K
dqt1    = 0    # kg/kg
dthldz1 = 3e-3  # K/m
dqtdz1  = -2.5e-6 # kg/kg/m
zi2     = 4000  # m

# Free troposphere
dthl2   = 4        # K
dqt2    = -2e-3    # kg/kg
dthldz2 = 10e-3  # K/m
dqtdz2  = 0     # kg/kg/m

# Profiles
z   = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl = np.zeros(z.size)
qt  = np.zeros(z.size)
u   = np.zeros(z.size)
v   = np.zeros(z.size)
ug  = np.zeros(z.size)
vg  = np.zeros(z.size)

for k in range(z.size):
    if (z[k] < zi1):
        thl[k] = thl0
        qt [k] = qt0
    elif (z[k] < zi2):
        thl[k] = thl0 + dthl1 + (z[k] - zi1) * dthldz1 
        qt [k] = qt0  + dqt1  + (z[k] - zi1) * dqtdz1
    else:
        thl[k] = thl0 + dthl1 + (zi2 - zi1) * dthldz1 + dthl2 + (z[k] - zi2) * dthldz2 
        qt [k] = qt0  + dqt1  + (zi2 - zi1) * dqtdz1  + dqt2  + (z[k] - zi2) * dqtdz2 

qt[qt<0] = 0 # :)

pl.figure()
pl.subplot(121)
pl.plot(thl,z)
pl.subplot(122)
pl.plot(qt,z)

# write the data to a file
proffile = open('rain.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s} {5:^20s} {6:^20s}\n'.format('z','thl','qt','u','ug','v','vg'))
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E} {6:1.14E}\n'.format(z[k], thl[k], qt[k], u[k], ug[k], v[k], vg[k]))
proffile.close()

ep = 287.04 / 461.5 

# Surface settings
def esat(T):
    c0 = 0.6105851e+03; c1 = 0.4440316e+02; c2 = 0.1430341e+01; c3 = 0.2641412e-01 
    c4 = 0.2995057e-03; c5 = 0.2031998e-05; c6 = 0.6936113e-08; c7 = 0.2564861e-11 
    c8 = -.3704404e-13 

    if (np.size(T)>1):
        x = T-273.15
        x[x<-80]=-80
    else:
        x  = max(-80.,T-273.15)
    return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

def qsat(p, T):
    return ep*esat(T)/(p-(1.-ep)*esat(T))

# Very rough estimation p
p = 1e5 - 9.81 * z
T = thl * (p/1e5)**(287/1004)
qs = qsat(p,T)

pl.figure()
pl.subplot(121)
pl.plot(qt, z, label='qt')
pl.plot(qs, z, label='qs')
pl.legend()

pl.subplot(122)
pl.plot((qt/qs)*100, z, label='qt')


#ps  = 101540.
#SST = 299.8 
#ths = SST / (ps/1.e5)**(287.04/1005.)
#qs  = qsat(ps, SST) 
#print('sbot[thl]=%f, sbot[qt]=%f'%(ths, qs))
#
#if(True):
#    # TMP: sounding UCLA-LES
#    ucla_ps  = np.array([     0,  740.,  3260.,    4000.  ])
#    ucla_ts  = np.array([ 297.9,  297.9,  312.6644, 317.0 ])
#    ucla_rts = np.array([  16.0,   13.8,    2.4,      1.8 ])
#    ucla_us  = np.array([  -9.9,   -8.42,  -3.38,    -1.9 ])
#    ucla_vs  = np.array([  -3.8,   -3.8,   -3.8,     -3.8 ])
#
#    qs_tmp = ucla_rts * 1e-3 # g/kg -> kg/kg
#    es_tmp = ps * (qs_tmp) / (-ep*qs_tmp + qs_tmp + ep)
#    ucla_rts2 = ep * es_tmp / (ps - es_tmp)
#
#    print(ucla_rts2)
#
#    import pylab as pl
#    pl.close('all')
#
#    pl.figure()
#    pl.subplot(331)
#    pl.plot(thl, z)
#    pl.plot(ucla_ts, ucla_ps, 'o')
#    
#    pl.subplot(332)
#    pl.plot(qt*1000., z)
#    pl.plot(ucla_rts, ucla_ps, 'o')
#    
#    pl.subplot(333)
#    pl.plot(u, z)
#    pl.plot(ucla_us, ucla_ps, 'o')
#    
#    pl.subplot(334)
#    pl.plot(v, z)
#    pl.plot(ucla_vs, ucla_ps, 'o')
#    
#    pl.subplot(335)
#    pl.plot(wls, z)
#    
#    pl.subplot(336)
#    pl.plot(thlls, z)
#    
#    pl.subplot(337)
#    pl.plot(qtls, z)
