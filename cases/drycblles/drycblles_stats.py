import xarray as xr
import matplotlib.pyplot as plt
import netCDF4 as nc
from matplotlib.dates import DateFormatter

# dir = '/data/lafe/20170808_col/64_cols/scalar/'
fname = 'drycblles.default.0000000.nc'
groups = [''] + list(nc.Dataset(fname).groups.keys())
prof_dict = {}
for group in groups:
    prof_dict[group] = xr.open_dataset(fname, group=group)
prof = xr.merge(prof_dict.values()).sel(z=slice(0,1500),zh=slice(0,1500),time=slice("2000-01-01T00:30:00",None))

prof.assign_coords(time = prof.time.dt.hour + prof.time.dt.minute/60)
profmn = prof.mean(dim='time', keep_attrs=True)

fig,axes = plt.subplots(2,2,figsize=(6,6),sharey=True)
i=-1
for var in ['u_2','v_2','w_2','th_2']:
    i+=1
    ax = axes.flatten()[i]
    try:
        profmn[var].plot(y='z', ax=ax)
    except:
        profmn[var].plot(y='zh',ax=ax)
    ax.set_ylabel('')
for i in range(2):
    axes[i,0].set_ylabel('z [m]')
fig.tight_layout()

for var in ['th_w','th_diff','th_flux']:
    profmn[var].plot(y='zh', label=var)
plt.legend()    

fig, axes = plt.subplots(2,3,figsize=(10,8),sharex=True)
dfmt = DateFormatter("%H:%M")
i=-1
for var in ['th','evisc','u_2','v_2','w_2','th_flux']:
    i+=1
    ax = axes.flatten()[i]
    prof[var].plot(x='time', ax=ax)
    prof['zi'].plot(x='time',ax=ax,ls='--')
    ax.set_ylabel('')
    ax.xaxis.set_major_formatter(dfmt)
for i in range(2):
    axes[i,0].set_ylabel('z [m]')
fig.tight_layout()
 