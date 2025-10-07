import matplotlib.pyplot as plt
import colormaps as cmaps
import xarray as xr
import numpy as np
import glob

plt.close('all')
#plt.ioff()
plt.ion()

ds = xr.open_dataset('particle_dump.0000000.h5', engine='h5netcdf')

n_plot = 500
cmap = cmaps.WhiteBlueGreenYellowRed
cc = cmap(np.linspace(0, 1, n_plot))

"""
max_length = 50
for t in range(0, 551, 1):
    print(t)
    fig=plt.figure(figsize=(19.2, 10.8))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax=plt.gca()
    ax.set_axis_off()

    t0 = max(0, t-max_length)
    t1 = t+2
    
    for i in range(n_plot):
        plt.plot(x[i,t0:t1], z[i,t0:t1], color=cc[i], alpha=0.5)

    plt.xlim(0, 3200)
    plt.ylim(0, 900)

    plt.savefig(f'figs/fig{t:05d}.png', dpi=100)

    plt.close('all')
"""

fig=plt.figure()
for i in range(n_plot):
    plt.plot(ds.x[:,i], ds.z[:,i], color=cc[i], alpha=1)
