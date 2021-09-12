import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

start = 0
end = 100
step = 10

stats = nc.Dataset("orlando.default.0010800.nc", "r")
t = stats.variables["time"][start:end]
z = stats.variables["z"][:]
zh = stats.variables["zh"][:]
st = stats.groups["thermo"].variables["thl"][start:end, :]
evisct = stats.groups["default"].variables["evisc"][start:end, :]
ut = stats.groups["default"].variables["u"][start:end, :]
vt = stats.groups["default"].variables["v"][start:end, :]
u2t = stats.groups["default"].variables["u_2"][start:end, :]
v2t = stats.groups["default"].variables["v_2"][start:end, :]
w2t = stats.groups["default"].variables["w_2"][start:end, :]
s2t = stats.groups["thermo"].variables["thl_2"][start:end, :]
sturbt = stats.groups["thermo"].variables["thl_w"][start:end, :]
sdifft = stats.groups["thermo"].variables["thl_diff"][start:end, :]
sfluxt = stats.groups["thermo"].variables["thl_flux"][start:end, :]
sgradt = stats.groups["thermo"].variables["thl_grad"][start:end, :]

s = np.mean(st, 0)
#evisc = np.mean(evisct, 0)

u2 = np.mean(u2t, 0)
u = np.mean(ut, 0)
v2 = np.mean(v2t, 0)
w2 = np.mean(w2t, 0)
s2 = np.mean(s2t, 0)

sturb = np.mean(sturbt, 0)
sdiff = np.mean(sdifft, 0)
sflux = np.mean(sfluxt, 0)

ht = np.zeros(t.size)
wstart = np.zeros(t.size)
for n in range(t.size):
    hindex = np.argmin(abs(sgradt[n, :] - max(sgradt[n, :])))
    ht[n] = z[hindex]
    wstart[n] = ((9.81/300.)*sfluxt[n, 0]*ht[n])**(1./3.)

"""
# enable LaTeX plotting
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
"""

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
    plt.plot(ut[n, :], z)
plt.xlabel(r'$u [m s^{-1}]$')
plt.ylabel(r'$z [m]$')

plt.figure()
for n in range(start, end, step):
    plt.plot(vt[n, :], z)
plt.xlabel(r'$v [m s^{-1}]$')
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

thermo = stats.groups['thermo']
thl = thermo.variables['thl'][:]
qr = thermo.variables['qr'][:]
qt = thermo.variables['qt'][:]
ql = thermo.variables['ql'][:]
ql_frac = thermo.variables['ql_frac'][:]
qflux = thermo.variables['qt_flux'][:]
sfluxt = thermo.variables["thl_flux"][:]
bx = thermo.variables['b'][:]


thl_2 = thermo.variables['thl_2'][:]
f,ax = plt.subplots(1)
for th in thl:
    ax.plot(th[:], z)
ax.set_ylabel('z (m)')
ax.set_xlabel('theta (K)')
f,ax = plt.subplots(1)
for q in qr:
    ax.plot(q[:], z)
ax.set_ylabel('z (m)')
ax.set_xlabel('qr (kg/kg)')
f,ax = plt.subplots(1)
for q in qt:
    ax.plot(1e3*q[:], z)
ax.set_ylabel('z (m)')
ax.set_xlabel('qt (g/kg)')
f,ax = plt.subplots(1)
for q in ql:
    ax.plot(1e3*q[:], z)
ax.set_ylabel('z (m)')
ax.set_xlabel('ql (g/kg)')
f,ax = plt.subplots(1)
for q in ql_frac:
    ax.plot(q[:], z)
ax.set_ylabel('z (m)')
ax.set_xlabel('ql_frac (-)')

f,ax = plt.subplots(1)
for q in bx:
    ax.plot(q[:], z)
ax.set_ylabel('z (m)')
ax.set_xlabel('Buoyancy (m/s2)')

f,ax = plt.subplots(1)
for q in sfluxt:
    ax.plot(q[:], zh)
ax.set_ylabel('z (m)')
ax.set_xlabel('Heat Flux (Km/s)')

f,ax = plt.subplots(1)
for q in qflux:
    ax.plot(1e3*q[:], zh)
ax.set_ylabel('z (m)')
ax.set_xlabel('Moisture Flux ((g/kg).(m/s)')
plt.show()

