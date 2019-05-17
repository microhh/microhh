import numpy

# Get number of vertical levels and size from .ini file
with open('drycblles.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

dthetadz = 0.003

# set the height
z    = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
th   = numpy.zeros(numpy.size(z))
thls = numpy.zeros(numpy.size(z))
wls  = numpy.zeros(numpy.size(z))

# linearly stratified profile
for k in range(kmax):
    th  [k] = 300. + dthetadz*z[k]
    thls[k] = 2.*(z[k]/zsize - 0.5) / 3600.
    wls [k] = -0.01*(z[k]/zsize)

"""
# well mixed profile with jump
h    = 1000.
dth  = 10.
dthz = 100.

for k in range(kmax):
    if(z[k] <= h - 0.5*dthz):
        th[k] = 300.
    elif(z[k] <= h + 0.5*dthz):
        th[k] = 300. + dth/dthz * (z[k]-(h-0.5*dthz))
    else:
        th[k] = 300. + dth + dthetadz*(z[k]-(h+0.5*dthz))
    thls[k] = 2.*(z[k]/zsize - 0.5) / 3600.
    wls [k] = -0.01*(z[k]/zsize)
"""

# write the data to a file
proffile = open('drycblles.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s}\n'.format('z','th','thls', 'wls'))
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E}\n'.format(z[k], th[k], thls[k], wls[k]))
proffile.close()
