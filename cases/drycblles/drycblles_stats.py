import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

start = 0
end = 19
step = 2

stats = nc.Dataset('drycblles.default.0000000.nc', 'r')
default = stats.groups['default']
thermo = stats.groups['thermo']

t = stats.variables['time'][start:end]
z = stats.variables['z'][:]
zh = stats.variables['zh'][:]

st = thermo.variables['th'][start:end, :]
evisct = default.variables['evisc'][start:end, :]
u2t = default.variables['u_2'][start:end, :]
v2t = default.variables['v_2'][start:end, :]
w2t = default.variables['w_2'][start:end, :]
s2t = thermo.variables['th_2'][start:end, :]
sturbt = thermo.variables['th_w'][start:end, :]
sdifft = thermo.variables['th_diff'][start:end, :]
sfluxt = thermo.variables['th_flux'][start:end, :]
sgradt = thermo.variables['th_grad'][start:end, :]

s = np.mean(st, axis=0)
evisc = np.mean(evisct, axis=0)

u2 = np.mean(u2t, axis=0)
v2 = np.mean(v2t, axis=0)
w2 = np.mean(w2t, axis=0)
s2 = np.mean(s2t, axis=0)

sturb = np.mean(sturbt, axis=0)
sdiff = np.mean(sdifft, axis=0)
sflux = np.mean(sfluxt, axis=0)

ht = np.zeros(t.size)
wstart = np.zeros(t.size)
for n in range(t.size):
    hindex = np.argmin(abs(sgradt[n, :] - max(sgradt[n, :])))
    ht[n] = z[hindex]
    wstart[n] = ((9.81/300.)*sfluxt[n, 0]*ht[n])**(1./3.)

plt.close('all')
plt.figure()
for n in range(start, end, step):
    plt.plot(st[n, :], z)
plt.xlabel(r'$\theta [K]$')
plt.ylabel(r'$z [m]$')

plt.figure()
for n in range(start, end, step):
    plt.plot(evisct[n, :], z)
plt.xlabel(r'$K_m [m^2 s^{-1}]$')
plt.ylabel(r'$z [m]$')

plt.figure()
for n in range(start, end, step):
    plt.plot(u2t[n, :], z)
plt.xlabel(r'$u^2 [m^2 s^{-2}]$')
plt.ylabel(r'$z [m]$')

plt.figure()
for n in range(start, end, step):
    plt.plot(w2t[n, :], zh)
plt.xlabel(r'$w^2 [m^2 s^{-2}]$')
plt.ylabel(r'$z [m]$')

plt.figure()
for n in range(start, end, step):
    plt.plot(sfluxt[n, :], zh)
plt.xlabel(r'$\overline{w\theta} [K m s^{-1}]$')
plt.ylabel(r'$z [m]$')

plt.figure()
plt.plot(np.mean(sfluxt[-5:-1, :], 0), zh, 'b-')
plt.plot(np.mean(sturbt[-5:-1, :], 0), zh, 'b--')
plt.plot(np.mean(sdifft[-5:-1, :], 0), zh, 'b:')
plt.ylim(0., 1500.)
plt.xlabel(r'$\overline{w\theta} [K m s^{-1}]$')
plt.ylabel(r'$z [m]$')

plt.figure()
plt.plot(t, ht)
plt.xlabel(r'$time [s]$')
plt.ylabel(r'$h [m]$')

plt.show()
