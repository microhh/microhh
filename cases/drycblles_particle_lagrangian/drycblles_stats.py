import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

plt.close('all')

# Read raw particle dump.
ds = xr.open_dataset('particle_dump.0000000.h5')

# Don't use stretched model grid, but equidistant grid for binning.
zsize = 3200
ktot = 128
dz = zsize / ktot

zh = np.arange(0, zsize+0.1, dz)
z = np.arange(dz/2, zsize, dz)

# Count particles per height bin.
c1, _ = np.histogram(ds.z[0 ], bins=zh)
c2, _ = np.histogram(ds.z[-1], bins=zh)

# Plot!
plt.figure()
plt.plot(c1, z, label='t=0 s')
plt.plot(c2, z, label='t=10800 s')
plt.legend()
plt.xlabel('# particles (-)')
plt.ylabel('z (m)')
