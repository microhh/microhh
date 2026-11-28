import matplotlib.pyplot as pl
import netCDF4 as nc4
import xarray as xr
import numpy as np

pl.close('all')

def xr_read_all(f, groups=['default'], decode_times=True):
    # Read all NetCDF groups into a single Dataset.
    dss = [xr.open_dataset(f, decode_times=decode_times)]
    for group in groups:
        dss.append(xr.open_dataset(f, group=group, decode_times=decode_times))
    return xr.merge(dss)

# Reference data windtunnel
u_wt,  z_wt1 = np.loadtxt('reference_data/u.txt', delimiter=',', unpack=True)
su_wt, z_wt2 = np.loadtxt('reference_data/sigma_u.txt', delimiter=',', unpack=True)
sw_wt, z_wt3 = np.loadtxt('reference_data/sigma_w.txt', delimiter=',', unpack=True)
ke_wt, z_wt4 = np.loadtxt('reference_data/ke.txt', delimiter=',', unpack=True)

# LES results
s1 = xr_read_all('canopy.default.0000000.nc', groups=['default', 'budget'], decode_times=False)
s1m = s1.sel(time=slice(1500, 2400)).mean(dim='time')

# Friction velocity at canopy top:
ustar_zc = ((s1m.u_flux[10]**2 + s1m.v_flux[10]**2))**(1/4)

# Canopy height:
zc = 10.

# Mean scaled velocity profile
pl.figure()
pl.plot(s1m.u/ustar_zc, s1m.z/zc, label='LES')
pl.scatter(u_wt, z_wt1, color='k', facecolor='none', label='wind tunnel')
pl.plot([0,8], [1,1], 'k:')
pl.xlabel(r'$u/u_*$ (-)')
pl.ylabel(r'$z/z_c$ (-)')
pl.xlim(0, 8)
pl.ylim(0, 3)
pl.legend()


# Variance profiles
pl.figure()
pl.subplot(221)
pl.plot(np.sqrt(s1m.u_2) / ustar_zc, s1m.z/zc)
pl.scatter(su_wt, z_wt2, color='k', facecolor='none')
pl.plot([0,8], [1,1], 'k:')
pl.xlabel(r'$\sigma_u/u_*$ (-)')
pl.ylabel(r'$z/z_c$ (-)')
pl.ylim(0, 6)
pl.xlim(0, 2.5)

pl.subplot(222)
pl.plot(np.sqrt(s1m.v_2) / ustar_zc, s1m.z/zc)
pl.xlabel(r'$\sigma_v/u_*$ (-)')
pl.ylabel(r'$z/z_c$ (-)')
pl.ylim(0, 6)
pl.xlim(0, 2)

pl.subplot(223)
pl.plot(np.sqrt(s1m.w_2) / ustar_zc, s1m.zh/zc)
pl.xlabel(r'$\sigma_w/u_*$ (-)')
pl.scatter(sw_wt, z_wt3, color='k', facecolor='none')
pl.plot([0,8], [1,1], 'k:')
pl.ylabel(r'$z/z_c$ (-)')
pl.ylim(0, 6)
pl.xlim(0, 1.5)

pl.subplot(224)
pl.plot(s1m.ke / ustar_zc**2, s1m.z/zc)
pl.scatter(ke_wt, z_wt4, color='k', facecolor='none')
pl.plot([0,8], [1,1], 'k:')
pl.xlabel(r'$KE/u_*^2$ (-)')
pl.ylabel(r'$z/z_c$ (-)')
pl.ylim(0, 6)
pl.xlim(0, 50)


# TKE budget
pl.figure()
pl.plot(s1m.tke_shear, s1m.z/zc, label='shear')
pl.plot(s1m.tke_turb, s1m.z/zc, label='turb')
pl.plot(s1m.tke_diff, s1m.z/zc, label='diff')
pl.plot(s1m.tke_pres, s1m.z/zc, label='pres')
pl.legend()
pl.ylim(0,4)
pl.ylabel(r'$z/z_c$ (-)')
pl.xlabel(r'dTKE/dt (m$^2$ s$^{-3}$)')
