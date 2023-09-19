import numpy as np
import netCDF4 as nc
import math
import matplotlib.pyplot as plt

# Get grid information from .ini file
with open('weisman_klemp.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])
        if(line.split('=')[0]=='jtot'):
            jmax = int(line.split('=')[1])
        if(line.split('=')[0]=='ysize'):
            ysize = float(line.split('=')[1])
        if(line.split('=')[0]=='itot'):
            imax = int(line.split('=')[1])
        if(line.split('=')[0]=='xsize'):
            xsize = float(line.split('=')[1])

dx = xsize / imax
dy = ysize / jmax
dz = zsize / kmax

# Bubble position (m)
xbub = 38400.
ybub = 0.5*ysize
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
float_type = np.float32

thl = np.fromfile("thl.0000000", dtype=float_type).reshape(kmax, jmax, imax)

for k in range(kminbub,kmaxbub):
    for j in range(jminbub,jmaxbub):
        for i in range(iminbub,imaxbub):
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
