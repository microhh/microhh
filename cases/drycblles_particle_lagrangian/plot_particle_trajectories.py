import matplotlib.pyplot as plt
import colormaps as cmaps
import numpy as np
import glob

plt.close('all')
#plt.ioff()
plt.ion()

float_type = np.float64

files = [f for f in glob.glob('particles.0*') if f != 'particles.0000000.h5']
files.sort()

n_time = len(files)
n_part = 1

x = np.zeros((n_part, n_time), dtype=float_type)
y = np.zeros((n_part, n_time), dtype=float_type)
z = np.zeros((n_part, n_time), dtype=float_type)

for t,f in enumerate(files):
    raw = np.fromfile(f)

    x[:,t] = raw[0:n_part]
    y[:,t] = raw[n_part:2*n_part]
    z[:,t] = raw[2*n_part:]

n_plot = 1
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
#plt.style.use('dark_background')
#fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
for i in range(n_plot):
    plt.plot(x[i,:], z[i,:], color='k', alpha=1)
