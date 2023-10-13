import matplotlib.pyplot as pl
import netCDF4 as nc4
import numpy as np
import xarray as xr

pl.close('all')

def xr_read_all(f, groups=None):

    # Get list of NetCDF groups:
    nc = nc4.Dataset(f)
    if groups is None:
        groups = list(nc.groups)

    # Check if NetCDF file has meaningful time units.
    if nc.variables['time'].units == 'seconds since start':
        decode_times = False
    else:
        decode_times = True

    nc.close()

    # Read all groups into a single Dataset.
    dss = [xr.open_dataset(f, decode_times=decode_times)]
    for group in groups:
        dss.append(xr.open_dataset(f, group=group, decode_times=decode_times))
    return xr.merge(dss)


if __name__ == '__main__':

    ds = xr_read_all('cabauw.default.0000000.nc')

    pl.figure(figsize=(8,7))

    nrows = 3
    ncols = 3
    sp = 1

    pl.subplot(nrows, ncols, sp); sp+=1
    pl.plot(ds.time, ds.rain_rate*3600)
    pl.ylabel(r'rain rate (mm h$^{-1}$)')

    step = 6*3600//300
    ysize = 7000
    nt = ds.time.size//step+1
    cc = pl.cm.turbo(np.linspace(0, 1, nt))

    pl.subplot(nrows, ncols, sp); sp+=1
    for i,t in enumerate(range(0, ds.time.size, step)):
       pl.plot(ds.ql[t,:]*1000, ds.z, color=cc[i], label=ds.time.values[t])
    pl.legend()
    pl.ylim(0, ysize)
    pl.ylabel('z (m)')
    pl.xlabel(r'$q_l$ (g kg$^{-1}$)')

    pl.subplot(nrows, ncols, sp); sp+=1
    for i,t in enumerate(range(0, ds.time.size, step)):
       pl.plot(ds.qi[t,:]*1000, ds.z, color=cc[i])
    pl.ylim(0, ysize)
    pl.ylabel('z (m)')
    pl.xlabel(r'$q_i$ (g kg$^{-1}$)')

    ax=pl.subplot(nrows, ncols, sp); sp+=1
    for i,t in enumerate(range(0, ds.time.size, step)):
       pl.plot(ds.qr[t,:]*1000, ds.z, color=cc[i])
    pl.ylim(0, ysize)
    pl.ylabel('z (m)')
    pl.xlabel(r'$q_r$ (g kg$^{-1}$)')

    pl.subplot(nrows, ncols, sp); sp+=1
    for i,t in enumerate(range(0, ds.time.size, step)):
       pl.plot(ds.qs[t,:]*1000, ds.z, color=cc[i])
    pl.ylim(0, ysize)
    pl.ylabel('z (m)')
    pl.xlabel(r'$q_s$ (g kg$^{-1}$)')

    pl.subplot(nrows, ncols, sp); sp+=1
    for i,t in enumerate(range(0, ds.time.size, step)):
       pl.plot(ds.qg[t,:]*1000, ds.z, color=cc[i])
    pl.ylim(0, ysize)
    pl.ylabel('z (m)')
    pl.xlabel(r'$q_g$ (g kg$^{-1}$)')

    pl.subplot(nrows, ncols, sp); sp+=1
    for i,t in enumerate(range(0, ds.time.size, step)):
       pl.plot(ds.vqr[t,:], ds.z, color=cc[i])
    pl.ylim(0, ysize)
    pl.ylabel('z (m)')
    pl.xlabel(r'$v_{qr}$ (m s$^{-1}$)')

    pl.subplot(nrows, ncols, sp); sp+=1
    pl.plot(ds.time, ds.qr[:,0]*1000, label='z=12.5 m')
    pl.plot(ds.time, ds.qr[:,10]*1000, label='z=300 m')
    pl.plot(ds.time, ds.qr[:,20]*1000, label='z=600 m')
    pl.xlabel('time')
    pl.ylabel(r'$q_r$ (g kg-1)')
    pl.legend()

    pl.tight_layout()
