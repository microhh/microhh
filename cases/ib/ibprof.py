import numpy as np

# Get number of vertical levels and size from .ini file
with open('ib.ini') as f:
    for line in f:
        if (line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if (line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

# Set the height.
z = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
u = 0.1*np.ones(np.size(z))

# Write the data to a file.
proffile = open('ib.prof','w')
proffile.write('{0:^20s} {1:^20s}\n'.format('z','u'))
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E}\n'.format(z[k], u[k]))
proffile.close()
