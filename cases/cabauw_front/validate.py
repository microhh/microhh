import netCDF4 as nc4
import xarray as xr
import numpy as np
import datetime

import matplotlib.pyplot as pl
import matplotlib.colors as mcolors
import matplotlib.dates as mdates

pl.close('all')

def xr_read_all(f, groups=None):
    """
    Read all (or selection of) NetCDF groups from NetCDF
    file using xarray. If `groups=None`, all NetCDF groups
    are loaded.
    """
    nc = nc4.Dataset(f)
    if groups is None:
        groups = list(nc.groups)

    # Check of NetCDF file has meaningful time units.
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

def format_ax(ax=None, dt=3):
    if ax is None:
        ax = pl.gca()

    minor_loc = mdates.HourLocator(byhour=np.arange(0,25))
    major_loc = mdates.HourLocator(byhour=np.arange(0,25,dt))

    ax.xaxis.set_minor_locator(minor_loc)
    ax.xaxis.set_major_locator(major_loc)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

def color_y_ax(which, color, ax=None):
    if ax is None:
        ax = pl.gca()

    ax.yaxis.label.set_color(color)
    ax.spines[which].set_color(color)
    ax.tick_params(axis='y', which='both', colors=color)

def classify_uhh(ql, qi, qr):
    out = np.zeros_like(ql)

    out[ql>1e-9] = 3
    out[qi>1e-9] = 4
    out[qr>1e-9] = 2

    return out

colors = {
        0: '#e0e0e0',   # clear sky
        1: '#d6604d',   # cloud droplets only
        2: '#a6d96a',   # drizzle/rain
        3: '#1a9850',   # drizzle/rain/cloud
        4: '#4575b4',   # ice
        5: '#fdae61',   # ice/supercooled droplets
        6: '#abd9e9',   # melting ice
        7: '#ffffff',   # melting ice/cloud droplets
        8: '#f1b6da',   # aerosol
        9: '#de77ae',   # insects
        10:'#c51b7d'}   # aerosol&insects

names = [
        'Clear sky','Cloud droplets','Drizzle/rain','Drizzle/rain & cloud',
        'Ice','Ice/supercooled droplets','Melting ice',
        'Melting ice/cloud droplets','Aerosols','Insects','Aerosols/insects']

cmap = mcolors.ListedColormap(colors.values())

date = datetime.datetime(2016, 8, 11)
#date = datetime.datetime(2017, 12, 11)

cb_path = '/home/scratch1/meteo_data/observations/cabauw/old/'
#cb_path = '/Users/bart/meteo/observations/cabauw/'

# Read the cloud net classification
cn = xr.open_dataset('{path}/cloudnet/{year:04d}{month:02d}{day:02d}_cabauw_classification.nc'.format(
        path=cb_path, year=date.year, month=date.month, day=date.day))

clw = xr.open_dataset('{path}/cloudnet/{year:04d}{month:02d}{day:02d}_cabauw_lwc-scaled-adiabatic.nc'.format(
        path=cb_path, year=date.year, month=date.month, day=date.day))

sfc = xr.open_dataset('{path}/surface_meteo_lc1/cesar_surface_meteo_lc1_t10_v1.0_{year:04d}{month:02d}.nc'.format(
        path=cb_path, year=date.year, month=date.month))

sfc = sfc.sel(time=slice(date, date+datetime.timedelta(hours=24)))

#
# MicroHH statistics
#
uhh = xr_read_all('cabauw.default.0000000.nc')


if True:
    """
    Pcolormesh of ql/qi/qr/etc.
    """

    pl.figure(figsize=(12,7))

    ylim = (0, 10000)
    fields = ['l', 'i', 'r', 's', 'g', 'h']

    ii = 1
    for t in ['q', 'n']:
        for f in fields:
            fld = '{}{}'.format(t,f)

            if t == 'q':
                thres = 1e-9
            else:
                thres = 1

            if fld not in ['nl', 'ni']:
                ax = pl.subplot(2, len(fields), ii); ii+=1
                pl.title(fld, loc='left')
                pl.pcolormesh(
                        uhh.time, uhh.z, uhh[fld].T,
                        norm=mcolors.PowerNorm(gamma=0.2),
                        cmap=pl.cm.gist_earth_r)
                pl.ylim(ylim)
                pl.colorbar(location='top', shrink=0.7)
                pl.grid()
                format_ax(ax, dt=6)
            else:
                ii+=1

    pl.tight_layout()


if True:
    """
    Precipitation rate vs obs
    """
    pl.figure()
    pl.subplot(121)
    pl.plot(uhh.time, uhh.rain_rate*3600, label='rain')
    pl.plot(uhh.time, uhh.snow_rate*3600, label='snow')
    pl.plot(uhh.time, uhh.graupel_rate*3600, label='graupel')
    pl.plot(uhh.time, uhh.hail_rate*3600, label='hail')
    pl.scatter(sfc.time, sfc.RAIN*6, label='Cabauw obs')
    pl.legend()
    pl.ylabel('rr (mm h-1)')

    pl.subplot(122)
    pl.plot(uhh.time, np.cumsum(uhh.rain_rate)*300, label='rain')
    pl.plot(uhh.time, np.cumsum(uhh.snow_rate)*300, label='snow')
    pl.plot(uhh.time, np.cumsum(uhh.graupel_rate)*300, label='graupel')
    pl.plot(uhh.time, np.cumsum(uhh.hail_rate)*300, label='hail')
    pl.scatter(sfc.time, np.cumsum(sfc.RAIN), label='Cabauw obs')
    pl.legend()
    pl.ylabel('sum(rr) (mm)')

    pl.legend()


if False:
    """
    rain/snow rate and cumulative values.
    """

    pl.figure()
    pl.subplot(121)
    pl.plot(uhh.time, uhh.snow_rate*3600*10)
    pl.ylabel('snow rate (mm h-1)')

    pl.subplot(122)
    pl.plot(uhh.time, np.cumsum(uhh.snow_rate)*300*10)
    pl.ylabel('accumulated snow (mm)')

    pl.tight_layout()




if False:
    """
    Cloudnet/MicroHH data/figure
    """

    fig=pl.figure(figsize=(8,5))

    ax = fig.add_axes([0.08, 0.4, 0.7, 0.55])
    pl.title('{}: observations'.format(date.isoformat()), loc='left', fontsize=8)
    pl.pcolormesh(cn.time, cn.height/1000., cn.target_classification.T, cmap=cmap, vmin=-0.5, vmax=10.5)
    pl.ylabel(r'$z$ (km)')
    pl.ylim(0,11)
    pl.xlim(date, date+datetime.timedelta(days=1))
    format_ax()
    ax.set_xticklabels([])

    cax = fig.add_axes([0.8, 0.45, 0.015, 0.45])
    cb=pl.colorbar(cax=cax, ticks=np.arange(0,11))
    cb.ax.set_yticklabels(names, fontsize=8)
    cb.ax.tick_params(which='both', size=0)

    ax1 = fig.add_axes([0.08, 0.08, 0.7, 0.24])

    pl.plot(clw.time, clw.lwp, color='tab:red')
    pl.ylabel(r'LWP (kg m$^{-2}$)')
    pl.xlim(date, date+datetime.timedelta(days=1))
    pl.ylim(0,3)

    ax2=ax1.twinx()
    pl.plot(sfc.time, sfc.RAIN/600*3600, color='tab:blue')
    pl.ylabel(r'rr (mm h$^{-1}$)')
    pl.xlim(date, date+datetime.timedelta(days=1))
    pl.ylim(0,5)
    ax2.spines['right'].set_visible(True)

    format_ax(ax1)
    color_y_ax('left', 'tab:red',  ax1)
    color_y_ax('right', 'tab:blue', ax2)
    ax2.spines['left'].set_color('tab:red')

    # MicroHH
    fig=pl.figure(figsize=(8,5))

    ax = fig.add_axes([0.08, 0.4, 0.7, 0.55])
    pl.title('{}: MicroHH + SB06'.format(date.isoformat()), loc='left', fontsize=8)
    pl.pcolormesh(uhh.time, uhh.z/1000., classify_uhh(
        uhh.ql, uhh.qi, uhh.qr).T, cmap=cmap, vmin=-0.5, vmax=10.5)
    pl.ylabel(r'$z$ (km)')
    pl.ylim(0,11)
    pl.xlim(date, date+datetime.timedelta(days=1))
    format_ax()
    ax.set_xticklabels([])

    cax = fig.add_axes([0.8, 0.45, 0.015, 0.45])
    cb=pl.colorbar(cax=cax, ticks=np.arange(0,11))
    cb.ax.set_yticklabels(names, fontsize=8)
    cb.ax.tick_params(which='both', size=0)

    ax1 = fig.add_axes([0.08, 0.08, 0.7, 0.24])

    pl.plot(uhh.time, uhh.ql_path, color='tab:red')
    pl.ylabel(r'LWP (kg m$^{-2}$)')
    pl.xlim(date, date+datetime.timedelta(days=1))
    pl.ylim(0,3)

    ax2=ax1.twinx()
    pl.plot(uhh.time, uhh.rain_rate*3600, color='tab:blue')
    pl.ylabel(r'rr (mm h$^{-1}$)')
    pl.xlim(date, date+datetime.timedelta(days=1))
    pl.ylim(0,5)
    ax2.spines['right'].set_visible(True)

    format_ax(ax1)
    color_y_ax('left', 'tab:red',  ax1)
    color_y_ax('right', 'tab:blue', ax2)
    ax2.spines['left'].set_color('tab:red')


if False:
    """
    Outline of ql/qi/qr/etc.
    """

    pl.figure(figsize=(10,5))

    ylim = (0, 10000)
    fields = ['l', 'i', 'r', 's', 'g', 'h']

    ii = 1
    for t in ['q', 'n']:
        for f in fields:
            fld = '{}{}'.format(t,f)

            if t == 'q':
                thres = 1e-9
            else:
                thres = 1

            if fld not in ['nl', 'ni']:
                ax = pl.subplot(2, len(fields), ii); ii+=1
                pl.title(fld, loc='left')
                pl.contour(uhh.time, uhh.z, uhh[fld].T>thres, [0.9, 1.1], colors='b')
                pl.ylim(ylim)
                format_ax(ax, dt=6)
            else:
                ii+=1

    pl.tight_layout()



