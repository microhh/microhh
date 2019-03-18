import numpy

# Get number of vertical levels and size from .ini file
with open('dycoms.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

# set the height
z     = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl   = numpy.zeros(numpy.size(z))
qt    = numpy.zeros(numpy.size(z))
u     = numpy.zeros(numpy.size(z))
ug    = numpy.zeros(numpy.size(z))
v     = numpy.zeros(numpy.size(z))
vg    = numpy.zeros(numpy.size(z))
wls   = numpy.zeros(numpy.size(z))
thlls = numpy.zeros(numpy.size(z))
qtls  = numpy.zeros(numpy.size(z))

for k in range(kmax):
    # temperature
    if(z[k] <= 840.):
        thl[k] = 289.0
    else:
        thl[k] = 297.5 + (z[k]-840.)** (1./3.)

    # specific humidity
    if(z[k] <= 840.):
        qt[k] = 1e-3*9.0
    else:
        qt[k] = 1.e-3*1.5

    wls[k] = -3.75E-6 * z[k]

    # u-wind component
    u[k] = 6

    # ug-wind component
    ug[k] = 7

    # u-wind component
    v[k] = -4.25

    # ug-wind component
    vg[k] = -5.5



# write the data to a file
proffile = open('dycoms.prof', 'w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s} {5:^20s} {6:^20s} {7:^20s}\n'.format('z', 'thl', 'qt', 'u', 'ug', 'v', 'vg', 'wls'))
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E} {6:1.14E} {7:1.14E}\n'.format(z[k], thl[k], qt[k], u[k], ug[k], v[k], vg[k], wls[k]))
proffile.close()
