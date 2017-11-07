import numpy as np
import matplotlib.pylab as pl

# Help function to write the MicroHH input files
def write_output(file_name, data):

    f    = open(file_name, 'w')
    size = data[list(data.keys())[0]].size

    # Write header
    for var in data.keys():
        f.write('{0:^21s} '.format(var))
    f.write('\n')

    # Write data
    for k in range(size):
        for var in data.keys():
            f.write('{0:+1.14E} '.format(data[var][k]))
        f.write('\n')

    f.close()

# Help function to write the MicroHH input files
def write_time_profs(file_name, z, times, data):
    """ Write time varying input profiles for MicroHH, e.g. large scale forcings """
    f = open(file_name, 'w')

    f.write('{0:^21s} '.format('z'))
    for time in times:
        f.write('{0:^21.0f} '.format(time))
    f.write('\n')

    for k in range(z.size):
        f.write('{0:+1.14E} '.format(z[k]))
        for i in range(times.size):
            f.write('{0:+1.14E} '.format(data[i,k]))
        f.write('\n')

    f.close()

# Get number of vertical levels and size from .ini file
with open('ls_force.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

# Vertical profiles
dz  = zsize / kmax
z   = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl = 300 + 0.003*z

# Time varying forcings
# ---------------------
# Subsidence
time_p    = np.array([0,5400,10800])
wls       = np.zeros((time_p.size, z.size))
wls[-1,:] = -1e-5 * z

# Large scale advection moisture
qtls       = np.zeros((time_p.size, z.size))
qtls[-1,:] = 0.0000005

# Nudging
nudgefac = np.zeros(z.size)
nudgefac[z>1000] = 1. 
snudge = np.zeros((time_p.size, z.size))
snudge[-1,:] = 1.

# Geostrophic wind
ug = np.zeros((time_p.size, z.size))
vg = np.zeros((time_p.size, z.size))
ug[-1,:] = 2.
vg[-1,:] = -2

# Time dependent surface properties
time_s  = np.linspace(0, 10800, 8)
wthlbot = 0.1 * np.sin(np.pi*time_s / 10800)
pbot    = np.linspace(1e5, 1.02e5, 8)

# Write input MicroHH
# ---------------------
# Vertical profiles
data = {'z':z, 'thl':thl, 'nudgefac':nudgefac}
write_output('ls_force.prof', data)

# Time dependent surface properties
data = {'t':time_s, 'sbot[thl]':wthlbot, 'pbot':pbot}
write_output('ls_force.time', data)

# Time varying subsidence
write_time_profs('wls.timeprof', z, time_p, wls)

# Time varying large scale advection
write_time_profs('snudge.timeprof', z, time_p, snudge)

# Time varying geostrophic wind
write_time_profs('ug.timeprof', z, time_p, ug)
write_time_profs('vg.timeprof', z, time_p, vg)
