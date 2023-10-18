import matplotlib.pyplot as pl
import xarray as xr
import numpy as np

pl.close('all')

ds_in = xr.open_dataset('dmin_wetgrowth_lookup.nc')

dT = 0.5
temp_out = np.arange(ds_in.T[0], ds_in.T[-1]+1e-3, dT)

# Set coordinates for interpolation.
ds_in = ds_in.assign_coords({
        'nqw': ds_in.qw,
        'nqi': ds_in.qi,
        'ntemp': ds_in.T,
        'npres': ds_in.p})

# Interpolate!
ds_out = ds_in.interp({'ntemp': temp_out})

"""
pl.figure()
pl.subplot(121)
pl.imshow(ds_in.Dmin_wetgrowth_table[0,0,:,:])

pl.subplot(122)
pl.imshow(ds_out.Dmin_wetgrowth_table[0,0,:,:])
"""

comp = dict(zlib=True, complevel=5)
encoding = {var: comp for var in ds_out.data_vars}

ds_out.to_netcdf('dmin_wetgrowth_lookup_61.nc', encoding=encoding)
