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
    pl.subplot(321)
    pl.plot(ds.time, ds.rr*3600)
    pl.ylabel(r'rain rate (mm h$^{-1}$)')

    step = 36
    ysize = 7000

    pl.subplot(322)
    for i,t in enumerate(range(0, ds.time.size, step)):
       pl.plot(ds.ql[t,:]*1000, ds.z, label=ds.time.values[t])
    pl.legend()
    pl.ylim(0, ysize)
    pl.ylabel('z (m)')
    pl.xlabel(r'$q_l$ (kg kg$^{-1}$)')

    pl.subplot(323)
    for i,t in enumerate(range(0, ds.time.size, step)):
       pl.plot(ds.qi[t,:]*1000, ds.z)
    pl.ylim(0, ysize)
    pl.ylabel('z (m)')
    pl.xlabel(r'$q_i$ (kg kg$^{-1}$)')

    pl.subplot(324)
    for i,t in enumerate(range(0, ds.time.size, step)):
       pl.plot(ds.qr[t,:]*1000, ds.z)
    pl.ylim(0, ysize)
    pl.ylabel('z (m)')
    pl.xlabel(r'$q_r$ (kg kg$^{-1}$)')

    pl.subplot(325)
    for i,t in enumerate(range(0, ds.time.size, step)):
       pl.plot(ds.qs[t,:]*1000, ds.z)
    pl.ylim(0, ysize)
    pl.ylabel('z (m)')
    pl.xlabel(r'$q_s$ (kg kg$^{-1}$)')

    pl.subplot(326)
    for i,t in enumerate(range(0, ds.time.size, step)):
       pl.plot(ds.qg[t,:]*1000, ds.z)
    pl.ylim(0, ysize)
    pl.ylabel('z (m)')
    pl.xlabel(r'$q_g$ (kg kg$^{-1}$)')

    pl.tight_layout()
