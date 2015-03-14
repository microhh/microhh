import numpy as np

# Get number of vertical levels and size from .ini file
with open('rico.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

# set the height
z     = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl   = np.zeros(z.size)
qt    = np.zeros(z.size)
u     = np.zeros(z.size)
ug    = np.zeros(z.size)
v     = np.zeros(z.size)
vg    = np.zeros(z.size)
wls   = np.zeros(z.size)
thlls = np.zeros(z.size)
qtls  = np.zeros(z.size)

for k in range(kmax):

    # Liquid water potential temperature
    if(z[k] < 740.):
        thl[k] = 297.9
    else:
        thl[k] = 297.9 + (317.0 - 297.9)/(4000. - 740.) * (z[k] - 740.) 

    # Specific humidity
    if(z[k] < 740.):
        qt[k] = 16.0 + (13.8 - 16.0) / 740. * z[k]
    elif(z[k] < 3260.):
        qt[k] = 13.8 + (2.4 - 13.8) / (3260. - 740.) * (z[k] - 740.) 
    else:
        qt[k] = 2.4 + (1.8 - 2.4)/(4000. - 3260.) * (z[k] - 3260.) 

    # Subsidence
    if(z[k] < 2260):
        wls[k] = -0.005 * (z[k] / 2260.)
    else:
        wls[k] = -0.005

    # U and V component wind
    u[k]  = -9.9 + 2.0e-3 * z[k]
    ug[k] = u[k]
    v[k]  = -3.8
    vg[k] = v[k]

    # Advective and radiative tendency thl
    thlls[k] = -2.5 / 86400.

    # Advective tendency qt
    if(z[k] < 2980):
        qtls[k] = -1.0 / 86400. + (1.3456/ 86400.) * z[k] / 2980.
    else:
        qtls[k] = 4e-6

# normalize profiles to SI
qt  /= 1000.
qtls/= 1000.

# write the data to a file
proffile = open('rico.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s} {5:^20s} {6:^20s} {7:^20s} {8:^20s} {9:^20s}\n'.format('z','thl','qt','u','ug','v','vg','wls','thlls','qtls'))
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E} {6:1.14E} {7:1.14E} {8:1.14E} {9:1.14E}\n'.format(z[k], thl[k], qt[k], u[k], ug[k], v[k], vg[k], wls[k], thlls[k], qtls[k]))
proffile.close()

# Surface settings
def esat(T):
    c0 = 0.6105851e+03; c1 = 0.4440316e+02; c2 = 0.1430341e+01; c3 = 0.2641412e-01 
    c4 = 0.2995057e-03; c5 = 0.2031998e-05; c6 = 0.6936113e-08; c7 = 0.2564861e-11 
    c8 = -.3704404e-13 
    x  = max(-80.,T-273.15)
    return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

def qsat(p, T):
    ep = 287.04 / 461.5 
    return ep*esat(T)/(p-(1-ep)*esat(T))

ps  = 101540.
SST = 299.8 
ths = SST / (ps/1.e5)**(287.04/1005.)
qs  = qsat(ps, SST) 
print('sbot[thl]=%f, sbot[qt]=%f'%(ths, qs))

if(False):
    # TMP: sounding UCLA-LES
    ucla_ps  = ([    0,  740., 3260., 4000.,  ])
    ucla_ts  = ([297.9, 297.9, 312.6644, 317.0,])
    ucla_rts = ([ 16.0,  13.8,   2.4,   1.8,   ])
    ucla_us  = ([ -9.9, -8.42,  -3.38, -1.9,   ])
    ucla_vs  = ([ -3.8,  -3.8,  -3.8,  -3.8,   ])

    import pylab as pl
    pl.close('all')

    pl.figure()
    pl.subplot(331)
    pl.plot(thl, z)
    pl.plot(ucla_ts, ucla_ps, 'o')
    
    pl.subplot(332)
    pl.plot(qt, z)
    pl.plot(ucla_rts, ucla_ps, 'o')
    
    pl.subplot(333)
    pl.plot(u, z)
    pl.plot(ucla_us, ucla_ps, 'o')
    
    pl.subplot(334)
    pl.plot(v, z)
    pl.plot(ucla_vs, ucla_ps, 'o')
    
    pl.subplot(335)
    pl.plot(wls, z)
    
    pl.subplot(336)
    pl.plot(thlls, z)
    
    pl.subplot(337)
    pl.plot(qtls, z)
