import numpy as np
import netCDF4 as nc
import math
import matplotlib.pyplot as plt
import os

# Get grid information from .ini file
with open('weisman_klemp.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            ktot = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])
        if(line.split('=')[0]=='jtot'):
            jtot = int(line.split('=')[1])
        if(line.split('=')[0]=='ysize'):
            ysize = float(line.split('=')[1])
        if(line.split('=')[0]=='itot'):
            itot = int(line.split('=')[1])
        if(line.split('=')[0]=='xsize'):
            xsize = float(line.split('=')[1])

dx = xsize / itot
dy = ysize / jtot
dz = zsize / ktot

# Bubble position (m)
xbub = (1./3.) * xsize
ybub = (1./2.) * ysize
zbub = 1400.

# Bubble size (m)
lxybub = 10000.
lzbub  = 1400.

# Bubble Amplitude (K)
bubamp = 2

iminbub = int(round((xbub - lxybub)/dx))
imaxbub = int(round((xbub + lxybub)/dx)) + 1
jminbub = int(round((ybub - lxybub)/dy))
jmaxbub = int(round((ybub + lxybub)/dy)) + 1
kminbub = int(round((zbub - lzbub)/dz))
kmaxbub = int(round((zbub + lzbub)/dz))  + 1

# Open the file
filename = "thl.0000000"
float_type = np.float64 if os.path.getsize(filename) // (2*itot + 2*jtot + 2*ktot) == 8 else np.float32

thl = np.fromfile(filename, dtype=float_type).reshape(ktot, jtot, itot)

for k in range(kminbub, kmaxbub):
    for j in range(jminbub, jmaxbub):
        for i in range(iminbub, imaxbub):
            dist = math.sqrt( ((xbub - i*dx)/lxybub)**2 + ((ybub - j*dy)/lxybub)**2 + ((zbub - k*dz)/lzbub)**2)
            if (dist < 1.0):
                thl[k, j, i] = thl[k, j, i] + bubamp * np.cos(dist*np.pi/2)**2 + bubamp * 0.01*(np.random.rand() - 0.5)

thl.tofile("thl.0000000")

"""
xh = np.arange(0.0, xsize, dx)
yh = np.arange(0.0, ysize, dy)

fig,ax=plt.subplots(1, 1)
cp = ax.contourf(xh, yh, thl[int(round(zbub/dz)), :, :])
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Filled Contours Plot')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
plt.show()
"""
