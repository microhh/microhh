import numpy

# Get number of vertical levels and size from .ini file
with open('eady.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

dthetadz = 0.004

# set the height
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)

fc = 1.e-3
dudz = 1e-2

# linearly stratified profile
th = 300. + dthetadz*z
u = dudz*z
ug = u.copy()
print("dthetady_ls = {0}".format(-dudz*fc))

# write the data to a file
proffile = open('eady.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s}\n'.format('z','th','u','ug'))
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E}\n'.format(z[k], th[k], u[k], ug[k]))
proffile.close()
