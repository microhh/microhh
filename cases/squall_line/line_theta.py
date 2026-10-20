import numpy as np
import netCDF4 as nc
import math
import matplotlib.pyplot as plt
import os
import microhh_tools as mht

# Get grid information from .ini file
ini = mht.Read_namelist('squall_line.ini')
ktot = ini['grid']['ktot']
zsize = ini['grid']['zsize']
jtot = ini['grid']['jtot']
ysize = ini['grid']['ysize']
itot = ini['grid']['itot']
xsize = ini['grid']['xsize']
rndseed = ini['fields']['rndseed']

dx = xsize / itot
dy = ysize / jtot

ds_input = nc.Dataset("squall_line_input.nc")
z = ds_input.variables['z'][:]

# line position (m)
xbub = (1./2.) * xsize
zbub = 2500.

# Bubble size (m)
# the bubble size is equal to the buble position as we start from (0, 0)
lxbub = xbub
lzbub  = zbub

# Bubble Amplitude (K)
bubamp = -5
bub_rnd = 0.1

iminbub = int(0)
ibub_1 = int((xbub/dx)/3)
ibub_2 = int((xbub/dx)/3*2)
imaxbub = int(xbub/dx + 1)
kminbub = int(0)
kmaxbub = np.where(z > zbub)[0][0]
print(kmaxbub)

# Open the file
filename = "thl.0000000"
float_type = np.float64 if os.path.getsize(filename) // (2*itot + 2*jtot + 2*ktot) == 8 else np.float32

thl = np.fromfile(filename, dtype=float_type).reshape(ktot, jtot, itot)

np.random.seed(rndseed)
rnd_array = np.random.rand(ktot, jtot, itot)

for k in range(kminbub, kmaxbub):
    for i in range(ibub_1, ibub_2):
        dist_vert = z[k]/lzbub
        dist_hor  = ((xbub/3)*2 - i*dx)/(lxbub/3)
        dist = math.sqrt(((xbub - i*dx)/lxbub)**2 + ((zbub - z[k])/lzbub)**2)
        if (dist_vert < 1.0):
            if (dist_hor < 1.0):
                for j in range(0, jtot):
                    thl[k, j, i] = thl[k, j, i] \
                                    + bubamp * (1 - dist_vert) * (1 - dist_hor) \
                                    + rnd_array[k, j, i] * bub_rnd
    for i in range(ibub_2, imaxbub):
        dist_vert = z[k]/lzbub
        dist = math.sqrt(((xbub - i*dx)/lxbub)**2 + ((zbub - z[k])/lzbub)**2)
        if (dist_vert < 1.0):
            for j in range(0, jtot):
                thl[k, j, i] = thl[k, j, i] \
                               + bubamp * (1 - dist_vert) \
                               + rnd_array[k, j, i] * bub_rnd

thl.tofile("thl.0000000")

# reduce moisture in the cold pool to avoid any condensation
filename = "qt.0000000"
float_type = np.float64 if os.path.getsize(filename) // (2*itot + 2*jtot + 2*ktot) == 8 else np.float32

qt = np.fromfile(filename, dtype=float_type).reshape(ktot, jtot, itot)

for k in range(kminbub, kmaxbub):
    for i in range(ibub_1, imaxbub):
        dist_vert = z[k]/lzbub
        if (dist_vert < 1.0):
            for j in range(0, jtot):
                qt[k, j, i] = qt[k, j, i] - 5e-3 * (1 - dist_vert)
                qt[k, j, i] = max(qt[k, j, i], 0)

    for i in range(iminbub, ibub_1):
        dist_vert = z[k]/lzbub
        dist_hor  = (xbub/3 - i*dx)/(lxbub/3)
        if (dist_vert < 1.0):
            if (dist_hor < 1.0):
                for j in range(0, jtot):
                    qt[k, j, i] = qt[k, j, i] - 5e-3 * (1 - dist_vert) * (1 - dist_hor)
                    qt[k, j, i] = max(qt[k, j, i], 0)

qt.tofile("qt.0000000")

# plot
xh = np.arange(0.0, xsize, dx)
yh = np.arange(0.0, ysize, dy)

fig, ax=plt.subplots(1, 1)
cp = ax.contourf(xh, yh, thl[0, :, :])
fig.colorbar(cp)    # Add a colorbar to a plot
ax.set_title('thl at z = 0')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')

fig2, ax=plt.subplots(1, 1)
cp = ax.contourf(xh, yh, qt[0, :, :])
fig2.colorbar(cp)    # Add a colorbar to a plot
ax.set_title('qt at z = 0')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
plt.show()

