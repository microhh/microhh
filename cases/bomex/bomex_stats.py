import xarray as xr
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

prof = {}
profmn = {}
plotens = False
cond_samps = ['default', 'ql', 'qlcore']
for cond_samp in cond_samps:
  fname = 'bomex.{}.0000000.nc'.format(cond_samp)

  groups = [''] + list(nc.Dataset(fname).groups.keys())
  prof_dict = {}
  for group in groups:
      prof_dict[group] = xr.open_dataset(fname, group=group, decode_times=False)
  prof[cond_samp] = xr.merge(prof_dict.values()).sel(time=slice(10800,None))
  prof[cond_samp]['TKE'] = 0.5*(prof[cond_samp]['u_2'] + prof[cond_samp]['v_2'] + prof[cond_samp]['w_2'])
  prof[cond_samp] = prof[cond_samp].where(prof[cond_samp].apply(np.fabs)<1e10)

  profmn[cond_samp] = prof[cond_samp].mean(dim='time', keep_attrs=True)

def plotfig(var):
  fig = plt.figure()
  for cs in cond_samps:
    if(plotens):
      prof[cs][var].plot.line(y=prof[cs][var].dims[1], color='#eeeeee')
    profmn[cs][var].plot(y=profmn[cs][var].dims[0], label=cs)
  plt.legend()

vars = ['area','thl', 'qt', 'thv', 'u', 'v', 'ql', 'ql_frac', 'thl_flux', 'qt_flux', 'thv_flux', 'u_2', 'v_2', 'w_2', 'thl_2', 'qt_2', 'thv_2']
for var in vars:
  plotfig(var)
