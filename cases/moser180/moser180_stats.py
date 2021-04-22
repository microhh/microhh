import numpy as np
import netCDF4 as nc
from matplotlib.pyplot import *

sample_size = 60
plotens = True

stats = nc.Dataset("moser180.default.0000000.nc","r")
end = len(stats.variables["time"][:])
stats.close()
start = max(0, end - sample_size)

# Read the viscosity.
with open('moser180.ini') as f:
    for line in f:
        if line.split('=')[0]=='visc':
            visc = float(line.split('=')[1])

# Read Moser's data.
Mosermean = np.loadtxt("chan180.means", skiprows=25)
Moserrey  = np.loadtxt("chan180.reystress", skiprows=25)
Moseru2   = np.loadtxt("chan180.uubal", skiprows=25)
Moserv2   = np.loadtxt("chan180.wwbal", skiprows=25)
Moserw2   = np.loadtxt("chan180.vvbal", skiprows=25)
Mosertke  = np.loadtxt("chan180.kbal" , skiprows=25)
Moseruw   = np.loadtxt("chan180.uvbal", skiprows=25)

yplusMoser = Mosermean[:,1]
uavgMoser  = Mosermean[:,2]
uvarMoser  = Moserrey[:,2]
vvarMoser  = Moserrey[:,3]
wvarMoser  = Moserrey[:,4]

u2_shearMoser = Moseru2[:,3]
u2_turbMoser  = Moseru2[:,6]
u2_viscMoser  = Moseru2[:,7]
u2_dissMoser  = Moseru2[:,2]
u2_rdstrMoser = Moseru2[:,4]

v2_turbMoser  = Moserv2[:,6]
v2_viscMoser  = Moserv2[:,7]
v2_dissMoser  = Moserv2[:,2]
v2_rdstrMoser = Moserv2[:,4]

w2_turbMoser  = Moserw2[:,6]
w2_viscMoser  = Moserw2[:,7]
w2_dissMoser  = Moserw2[:,2]
w2_presMoser  = Moserw2[:,5]
w2_rdstrMoser = Moserw2[:,4]

tke_shearMoser = Mosertke[:,3]
tke_turbMoser  = Mosertke[:,6]
tke_viscMoser  = Mosertke[:,7]
tke_dissMoser  = Mosertke[:,2]
tke_presMoser  = Mosertke[:,5]

uw_shearMoser = Moseruw[:,3]
uw_presMoser  = Moseruw[:,5]
uw_turbMoser  = Moseruw[:,6]
uw_viscMoser  = Moseruw[:,7]
uw_dissMoser  = Moseruw[:,2]
uw_rdstrMoser = Moseruw[:,4]

stats = nc.Dataset("moser180.default.0000000.nc","r")
t  = stats.variables["time"] [start:end]
z  = stats.variables["z"] [:]
zh = stats.variables["zh"][:]
uavgt = stats.groups["default"].variables["u"] [start:end,:]
vavgt = stats.groups["default"].variables["v"] [start:end,:]
uvart = stats.groups["default"].variables["u_2"][start:end,:]
vvart = stats.groups["default"].variables["v_2"][start:end,:]
wvart = stats.groups["default"].variables["w_2"][start:end,:]

uwt    = stats.groups["default"].variables["u_w"]   [start:end,:]
udifft = stats.groups["default"].variables["u_diff"][start:end,:]
ufluxt = stats.groups["default"].variables["u_flux"][start:end,:]

# variance budgets
u2_sheart = stats.groups["budget"].variables["u2_shear"][start:end,:]
u2_turbt  = stats.groups["budget"].variables["u2_turb"] [start:end,:]
u2_visct  = stats.groups["budget"].variables["u2_visc"] [start:end,:]
u2_disst  = stats.groups["budget"].variables["u2_diss"] [start:end,:]
u2_rdstrt = stats.groups["budget"].variables["u2_rdstr"][start:end,:]

v2_turbt  = stats.groups["budget"].variables["v2_turb"] [start:end,:]
v2_visct  = stats.groups["budget"].variables["v2_visc"] [start:end,:]
v2_disst  = stats.groups["budget"].variables["v2_diss"] [start:end,:]
v2_rdstrt = stats.groups["budget"].variables["v2_rdstr"][start:end,:]

w2_turbt  = stats.groups["budget"].variables["w2_turb"] [start:end,:]
w2_visct  = stats.groups["budget"].variables["w2_visc"] [start:end,:]
w2_disst  = stats.groups["budget"].variables["w2_diss"] [start:end,:]
w2_prest  = stats.groups["budget"].variables["w2_pres"] [start:end,:]
w2_rdstrt = stats.groups["budget"].variables["w2_rdstr"][start:end,:]

# tke budget
tke_sheart = stats.groups["budget"].variables["tke_shear"][start:end,:]
tke_turbt  = stats.groups["budget"].variables["tke_turb"] [start:end,:]
tke_visct  = stats.groups["budget"].variables["tke_visc"] [start:end,:]
tke_disst  = stats.groups["budget"].variables["tke_diss"] [start:end,:]
tke_prest  = stats.groups["budget"].variables["tke_pres"] [start:end,:]

# stress budgets
uw_sheart = stats.groups["budget"].variables["uw_shear"][start:end,:]
uw_visct  = stats.groups["budget"].variables["uw_visc"] [start:end,:]
uw_turbt  = stats.groups["budget"].variables["uw_turb"] [start:end,:]
uw_disst  = stats.groups["budget"].variables["uw_diss"] [start:end,:]
uw_prest  = stats.groups["budget"].variables["uw_pres"] [start:end,:]
uw_rdstrt = stats.groups["budget"].variables["uw_rdstr"][start:end,:]

utotavgt = (uavgt**2. + vavgt**2.)**.5
ustart = (ufluxt[:,0]**2.)**.25

uavg = np.mean(uavgt,0)
vavg = np.mean(vavgt,0)
uvar = np.mean(uvart,0)
vvar = np.mean(vvart,0)
wvar = np.mean(wvart,0)

uw    = np.mean(uwt,0)
udiff = np.mean(udifft,0)
uflux = np.mean(ufluxt,0)

u2_shear = np.mean(u2_sheart,0)
u2_turb  = np.mean(u2_turbt ,0)
u2_visc  = np.mean(u2_visct ,0)
u2_diss  = np.mean(u2_disst ,0)
u2_rdstr = np.mean(u2_rdstrt,0)
u2_resid = u2_shear + u2_turb + u2_visc + u2_diss + u2_rdstr

v2_turb  = np.mean(v2_turbt ,0)
v2_visc  = np.mean(v2_visct ,0)
v2_diss  = np.mean(v2_disst ,0)
v2_rdstr = np.mean(v2_rdstrt,0)
v2_resid = v2_turb + v2_visc + v2_diss + v2_rdstr

w2_turb  = np.mean(w2_turbt ,0)
w2_visc  = np.mean(w2_visct ,0)
w2_diss  = np.mean(w2_disst ,0)
w2_pres  = np.mean(w2_prest ,0)
w2_rdstr = np.mean(w2_rdstrt,0)
w2_resid = w2_turb + w2_visc + w2_diss + w2_pres + w2_rdstr

tke_shear = np.mean(tke_sheart,0)
tke_turb  = np.mean(tke_turbt ,0)
tke_visc  = np.mean(tke_visct ,0)
tke_diss  = np.mean(tke_disst ,0)
tke_pres  = np.mean(tke_prest ,0)
tke_resid = tke_shear + tke_turb + tke_visc + tke_diss + tke_pres

uw_shear = np.mean(uw_sheart,0)
uw_turb  = np.mean(uw_turbt ,0)
uw_visc  = np.mean(uw_visct ,0)
uw_diss  = np.mean(uw_disst ,0)
uw_pres  = np.mean(uw_prest ,0)
uw_rdstr = np.mean(uw_rdstrt,0)
uw_resid = uw_shear + uw_turb + uw_visc + uw_diss + uw_pres + uw_rdstr

utotavg = np.mean(utotavgt,0)
ustar   = np.mean(ustart)

print('Re_tau = %.2f' % (ustar / visc))

# create the theoretical lines
ypluslin = np.arange(0.5,15., 0.1)
ypluslog = np.arange(5.,800., 1.)
ulin = ypluslin
ulog = 2.5 * np.log( ypluslog ) + 5.

yplus  = z  * ustar / visc
yplush = zh * ustar / visc

starty = 0
endy   = z.size // 2

close('all')
figure()
if plotens:
    for n in range(end-start):
        semilogx(yplus[starty:endy], utotavgt[n,starty:endy] / ustar, color='#cccccc')
semilogx(yplus[starty:endy], utotavg[starty:endy] / ustar, 'bo-', label='u')
semilogx(yplusMoser, uavgMoser, 'k--', label="Moser")
semilogx(ypluslin, ulin, 'k:')
semilogx(ypluslog, ulog, 'k:')
xlabel('y+')
ylabel('u+')
legend(loc=2, frameon=False)
grid()
axis([0.1, 800, 0, 25])

figure()
if plotens:
    for n in range(end-start):
        plot(yplus [starty:endy], (uvart[n,starty:endy] / ustar**2.)**0.5, color='#cccccc')
        plot(yplus [starty:endy], (vvart[n,starty:endy] / ustar**2.)**0.5, color='#cccccc')
        plot(yplush[starty:endy], (wvart[n,starty:endy] / ustar**2.)**0.5, color='#cccccc')
plot(yplus [starty:endy], (uvar[starty:endy] / ustar**2.)**0.5, 'b-', label='u')
plot(yplus [starty:endy], (vvar[starty:endy] / ustar**2.)**0.5, 'g-', label='v')
plot(yplush[starty:endy], (wvar[starty:endy] / ustar**2.)**0.5, 'r-', label='w')
plot(yplusMoser, np.sqrt(uvarMoser), 'k--', label="Moser")
plot(yplusMoser, np.sqrt(vvarMoser), 'k--')
plot(yplusMoser, np.sqrt(wvarMoser), 'k--')
xlabel('y+')
ylabel('rms')
legend(loc=0, frameon=False)
grid()
axis([0, 120, 0, 4.0])

figure()
if plotens:
    for n in range(end-start):
        plot(yplus[starty:endy], u2_sheart[n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], u2_turbt [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], u2_visct [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], u2_disst [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], u2_rdstrt[n,starty:endy] * visc / ustar**4., color='#cccccc')
plot(yplus[starty:endy], u2_shear[starty:endy] * visc / ustar**4., 'b-', label='S')
plot(yplus[starty:endy], u2_turb [starty:endy] * visc / ustar**4., 'g-', label='Tt')
plot(yplus[starty:endy], u2_visc [starty:endy] * visc / ustar**4., 'c-', label='Tv')
plot(yplus[starty:endy], u2_diss [starty:endy] * visc / ustar**4., 'r-', label='D')
plot(yplus[starty:endy], u2_rdstr[starty:endy] * visc / ustar**4., 'm-', label='P')
plot(yplus[starty:endy], u2_resid[starty:endy] * visc / ustar**4., 'k-', label='resid')
plot(yplusMoser, u2_shearMoser, 'k--', label="Moser")
plot(yplusMoser, u2_turbMoser , 'k--')
plot(yplusMoser, u2_viscMoser , 'k--')
plot(yplusMoser, u2_dissMoser , 'k--')
plot(yplusMoser, u2_rdstrMoser, 'k--')
xlabel('y+')
ylabel('Rxx')
legend(loc=0, frameon=False)
grid()
axis([0, 200, -0.6, 0.6])

figure()
if plotens:
    for n in range(end-start):
        plot(yplus[starty:endy], v2_turbt [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], v2_visct [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], v2_disst [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], v2_rdstrt[n,starty:endy] * visc / ustar**4., color='#cccccc')
plot(yplus[starty:endy], v2_turb [starty:endy] * visc / ustar**4., 'g-', label='Tt')
plot(yplus[starty:endy], v2_visc [starty:endy] * visc / ustar**4., 'c-', label='Tv')
plot(yplus[starty:endy], v2_diss [starty:endy] * visc / ustar**4., 'r-', label='D')
plot(yplus[starty:endy], v2_rdstr[starty:endy] * visc / ustar**4., 'm-', label='P')
plot(yplus[starty:endy], v2_resid[starty:endy] * visc / ustar**4., 'k-', label='resid')
plot(yplusMoser, v2_turbMoser , 'k--', label="Moser")
plot(yplusMoser, v2_viscMoser , 'k--')
plot(yplusMoser, v2_dissMoser , 'k--')
plot(yplusMoser, v2_rdstrMoser, 'k--')
xlabel('y+')
ylabel('Ryy')
legend(loc=0, frameon=False)
grid()
axis([0, 200, -0.06, 0.06])

figure()
if plotens:
    for n in range(end-start):
        plot(yplush[starty:endy], w2_turbt [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplush[starty:endy], w2_visct [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplush[starty:endy], w2_disst [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplush[starty:endy], w2_prest [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplush[starty:endy], w2_rdstrt[n,starty:endy] * visc / ustar**4., color='#cccccc')
plot(yplush[starty:endy], w2_turb [starty:endy] * visc / ustar**4., 'g-', label='Tt')
plot(yplush[starty:endy], w2_visc [starty:endy] * visc / ustar**4., 'c-', label='Tv')
plot(yplush[starty:endy], w2_diss [starty:endy] * visc / ustar**4., 'r-', label='D')
plot(yplush[starty:endy], w2_pres [starty:endy] * visc / ustar**4., 'y-', label='Tp')
plot(yplush[starty:endy], w2_rdstr[starty:endy] * visc / ustar**4., 'm-', label='P')
plot(yplush[starty:endy], w2_resid[starty:endy] * visc / ustar**4., 'k-', label='resid')
plot(yplusMoser, w2_turbMoser , 'k--', label="Moser")
plot(yplusMoser, w2_viscMoser , 'k--')
plot(yplusMoser, w2_dissMoser , 'k--')
plot(yplusMoser, w2_presMoser , 'k--')
plot(yplusMoser, w2_rdstrMoser, 'k--')
xlabel('y+')
ylabel('Rzz')
legend(loc=0, frameon=False)
grid()
axis([0, 200, -0.06, 0.06])

figure()
if plotens:
    for n in range(end-start):
        plot(yplus[starty:endy], tke_sheart[n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], tke_turbt [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], tke_visct [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], tke_disst [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], tke_prest [n,starty:endy] * visc / ustar**4., color='#cccccc')
plot(yplus[starty:endy], tke_shear[starty:endy] * visc / ustar**4., 'b-', label='S')
plot(yplus[starty:endy], tke_turb [starty:endy] * visc / ustar**4., 'g-', label='Tt')
plot(yplus[starty:endy], tke_visc [starty:endy] * visc / ustar**4., 'c-', label='Tv')
plot(yplus[starty:endy], tke_diss [starty:endy] * visc / ustar**4., 'r-', label='D')
plot(yplus[starty:endy], tke_pres [starty:endy] * visc / ustar**4., 'y-', label='Tp')
plot(yplus[starty:endy], tke_resid[starty:endy] * visc / ustar**4., 'k-', label='resid')
plot(yplusMoser, tke_shearMoser, 'k--', label="Moser")
plot(yplusMoser, tke_turbMoser , 'k--')
plot(yplusMoser, tke_viscMoser , 'k--')
plot(yplusMoser, tke_dissMoser , 'k--')
plot(yplusMoser, tke_presMoser , 'k--')
xlabel('y+')
ylabel('tke')
legend(loc=0, frameon=False)
grid()
axis([0, 200, -0.3, 0.3])

figure()
if plotens:
    for n in range(end-start):
        plot(yplus[starty:endy], uw_sheart[n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], uw_turbt [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], uw_visct [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], uw_disst [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], uw_prest [n,starty:endy] * visc / ustar**4., color='#cccccc')
        plot(yplus[starty:endy], uw_rdstrt[n,starty:endy] * visc / ustar**4., color='#cccccc')
plot(yplus[starty:endy], uw_shear[starty:endy] * visc / ustar**4., 'b-', label='S')
plot(yplus[starty:endy], uw_turb [starty:endy] * visc / ustar**4., 'g-', label='Tt')
plot(yplus[starty:endy], uw_visc [starty:endy] * visc / ustar**4., 'c-', label='Tv')
plot(yplus[starty:endy], uw_diss [starty:endy] * visc / ustar**4., 'r-', label='D')
plot(yplus[starty:endy], uw_pres [starty:endy] * visc / ustar**4., 'y-', label='Tp')
plot(yplus[starty:endy], uw_rdstr[starty:endy] * visc / ustar**4., 'm-', label='P')
plot(yplus[starty:endy], uw_resid[starty:endy] * visc / ustar**4., 'k-', label='resid')
plot(yplusMoser, uw_shearMoser, 'k--', label="Moser")
plot(yplusMoser, uw_turbMoser , 'k--')
plot(yplusMoser, uw_viscMoser , 'k--')
plot(yplusMoser, uw_dissMoser , 'k--')
plot(yplusMoser, uw_presMoser , 'k--')
plot(yplusMoser, uw_rdstrMoser, 'k--')
xlabel('y+')
ylabel('Rxz')
legend(loc=0, frameon=False)
grid()
axis([0, 200, -0.1, 0.1])

figure()
if plotens:
    for n in range(end-start):
        plot(zh, ufluxt[n,:] / ustar**2., color='#cccccc')
        plot(zh, uwt   [n,:] / ustar**2., color='#cccccc')
        plot(zh, udifft[n,:] / ustar**2., color='#cccccc')
plot(zh, uflux / ustar**2., 'b-' , label='total flux')
plot(zh, uw    / ustar**2., 'b--', label='turbulent flux')
plot(zh, udiff / ustar**2., 'b:' , label='diffusive flux')
xlabel('y+')
ylabel('uflux')
legend(loc=0, frameon=False)
grid()
axis([0., 2., -1.1, 1.1])
show()
