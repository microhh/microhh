import numpy as np
import math
import matplotlib.pylab as pl

pl.ion()
pl.close('all')

from microhh_tools import *

# Class to write field to restart file
def Write_field(data, path, per_slice=True):
    if(per_slice): # Write level by level (less memory hungry.....)
        nz = data.shape[0]
        ny = data.shape[1]
        nx = data.shape[2]
        fout  = open(path, "wb")
        for k in range(nz):
            tmp  = data[k,:,:].reshape(nx*ny)
            tmp2 = st.pack('{0}{1}d'.format('<', tmp.size), *tmp) 
            fout.write(tmp2)
        fout.close()
    else: # Write entire field at once (memory hungry....)
        tmp  = data.reshape(data.size)
        tmp2 = st.pack('{0}{1}d'.format('<', tmp.size), *tmp) 
        fout = open(path, "wb")
        fout.write(tmp2)
        fout.close()  

def Patch(itot, jtot, x, y, patch_xoffs, patch_yoffs, patch_xh, patch_yh, patch_xr, patch_yr, patch_xi, patch_yi, patch_dim):
    pattern = np.zeros((jtot, itot))

    for i in range(itot):
        for j in range(jtot):
            xmod = math.fmod(x[i]-patch_xoffs, patch_xh) 
            ymod = math.fmod(y[j]-patch_yoffs, patch_yh)
            
            errvalx = 0.5 - 0.5*math.erf(2.*(abs(2.*xmod - patch_xh) - patch_xr) / patch_xi)
            errvaly = 0.5 - 0.5*math.erf(2.*(abs(2.*ymod - patch_yh) - patch_yr) / patch_yi) if patch_dim == 2 else 1

            pattern[j,i] = errvalx*errvaly;

    return pattern

# Read namelist and grid
nl = Read_namelist()
gr = Read_grid(nl.grid.itot, nl.grid.jtot, nl.grid.ktot, nl.grid.zsize)

# Calculate patch
patch = Patch(nl.grid.itot, nl.grid.jtot, gr.x, gr.y, 0, 0, nl.boundary.patch_xh, nl.boundary.patch_yh,
              nl.boundary.patch_xr, nl.boundary.patch_yr, nl.boundary.patch_xi, nl.boundary.patch_yi, nl.boundary.patch_dim)

pl.figure()
pl.plot(patch[:,64])

# Case settings
zi0  = 1000     # initial BL height
th0  = 290      # initial BL potential temperature
dth0 = 2        # initial BL  "    "      "    "   jump
gth  = 6.e-3     # potential temperature lapse rate

qt0  = 7.e-3     # specific humidity
dqt0 = -2e-3
gqt  = -2e-6

# Define base state temperature and moisture profile
th = np.zeros(gr.ktot)
qt = np.zeros(gr.ktot)
for k in range(gr.ktot):
    if(gr.z[k] < zi0):
        th[k] = th0
        qt[k] = qt0
    else:
        th[k] = th0 + dth0 + (gr.z[k] - zi0) * gth
        qt[k] = qt0 + dqt0 + (gr.z[k] - zi0) * gqt

# Very rough estimate pressure
p  = nl.thermo.pbot - 9.81 * gr.z
T  = th * exner(p)
qs = qsat(p,T)
RH = qt/qs

# Rural temperatue profile
th2 = th.copy()
th2[gr.z<150] -= 10.

T2 = th2 * exner(p)
qs2 = qsat(p,T2)
RH2 = qt/qs2

pl.figure()
pl.subplot(131)
pl.plot(th, gr.z)
pl.plot(th2, gr.z)
pl.subplot(132)
pl.plot(qt*1000, gr.z)
pl.subplot(133)
pl.plot(RH*100, gr.z)
pl.plot(RH2*100, gr.z)

# Write profiles
# write the data to a file
proffile = open('foggy_city.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s}\n'.format('z','thl','qt','u','v','ug','vg'))
for k in range(gr.z.size):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E}\n'.format(gr.z[k], th[k], qt[k]))
proffile.close()

# Write restart files
th3d = np.zeros((nl.grid.ktot, nl.grid.jtot, nl.grid.itot))
qt3d = np.zeros((nl.grid.ktot, nl.grid.jtot, nl.grid.itot))

for i in range(gr.itot):
    for j in range(gr.jtot):
        th3d[:,j,i] = patch[j,i] * th + (1-patch[j,i]) * th2
        qt3d[:,j,i] = qt

Write_field(th3d, 'thl.0000000')
Write_field(qt3d, 'qt.0000000')
