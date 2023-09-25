import xarray as xr
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

prof = {}
profmn = {}

cond_samps = ['default']
plotens = True

for cond_samp in cond_samps:
    fname = 'arm.{}.0000000.nc'.format(cond_samp)

    groups = [''] + list(nc.Dataset(fname).groups.keys())
    prof_dict = {}
    for group in groups:
        prof_dict[group] = xr.open_dataset(fname, group=group, decode_times=False)
    prof[cond_samp] = xr.merge(prof_dict.values()).sel(time=slice(10800,None))
    prof[cond_samp]['TKE'] = 0.5*(prof[cond_samp]['u_2'] + prof[cond_samp]['v_2'] + prof[cond_samp]['w_2'])
    prof[cond_samp] = prof[cond_samp].where(prof[cond_samp].apply(np.fabs)<1e10)
    for var in ['ql', 'ql_frac']:
       prof[cond_samp][var] = prof[cond_samp][var].where(prof[cond_samp][var] > 0)
    profmn[cond_samp] = prof[cond_samp].mean(dim='time', keep_attrs=True)

def plotprof(var):
    col = ['k', 'r', 'b', 'g', 'm', 'c']
    fig = plt.figure()
    for i, cs in enumerate(cond_samps):
        if(plotens):
            prof[cs][var].plot.line(y=prof[cs][var].dims[1], color=col[i], alpha=0.1)
        profmn[cs][var].plot(y=profmn[cs][var].dims[0], color = col[i], label=cs)
    plt.legend()


vars = ['thl', 'qt', 'thv', 'u', 'v', 'ql', 'ql_frac', 'thl_flux', 'qt_flux', 'thv_flux', 'u_2', 'v_2', 'w_2', 'thl_2', 'qt_2', 'thv_2']
cond_samp = 'default'
for var in vars:
    fig = plt.figure()
    prof[cond_samp][var].plot(y=prof[cond_samp][var].dims[1])
