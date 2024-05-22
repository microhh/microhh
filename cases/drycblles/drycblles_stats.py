import xarray as xr
import matplotlib.pyplot as plt
import netCDF4 as nc

fname = 'drycblles.default.0000000.nc'
groups = [''] + list(nc.Dataset(fname).groups.keys())
prof_dict = {}
for group in groups:
    prof_dict[group] = xr.open_dataset(fname, group=group, decode_times=False)
prof = xr.merge(prof_dict.values()).sel(z=slice(0,1500),zh=slice(0,1500),time=slice(1800,None))

profmn = prof.mean(dim='time', keep_attrs=True)

fig,axes = plt.subplots(2,2,figsize=(6,6),sharey=True)
vars = ['u_2','v_2','w_2','th_2']
for i, var in enumerate(vars):
    ax = axes.flatten()[i]
    profmn[var].plot(y=profmn[var].dims[0], ax=ax)
    ax.set_ylabel('')

for i in range(2):
    axes[i,0].set_ylabel('z [m]')
fig.tight_layout()

plt.figure()
for var in ['th_w','th_diff','th_flux']:
    profmn[var].plot(y=profmn[var].dims[0], label=var)
plt.legend()
plt.tight_layout()    

fig, axes = plt.subplots(2,3,figsize=(10,8),sharex=True)
vars = ['th','evisc','u_2','v_2','w_2','th_flux']
for i, var in enumerate(vars):
    ax = axes.flatten()[i]
    prof[var].plot(x='time', ax=ax)
    prof['zi'].plot(x='time', ax=ax,ls='--')
    ax.set_ylabel('')
for i in range(2):
    axes[i,0].set_ylabel('z [m]')
fig.tight_layout()

plt.show()
