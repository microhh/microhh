import matplotlib.pyplot as pl
import xarray as xr
import numpy as np

pl.close('all')

# Original LUT from ICON repository.
include_file = 'phillips_nucleation_2010.incf'

ni = 101
nj = 101

arrays = {
    'afrac_dust': np.zeros((ni, nj), dtype=np.float64),
    'afrac_soot': np.zeros((ni, nj), dtype=np.float64),
    'afrac_orga': np.zeros((ni, nj), dtype=np.float64)}

with open(include_file) as f:
    for l in f.readlines():

        # (i,j) indices.
        # NOTE: In Fortran, the arrays are allocated as (0:100),
        #       so no conversion to C indexing needed...
        indices = l.split('(')[1].split(')')[0]
        i0 = int(indices.split(',')[0].split(':')[0])
        i1 = int(indices.split(',')[0].split(':')[1])+1
        j  = int(indices.split(',')[1])

        # Name.
        name = l.split(' ')[1].split('(')[0]

        # Data.
        data = l.split('/')[1][:-1].split(',')
        data = np.array([float(x) for x in data])

        # Write into arrays.
        arrays[name][i0:i1, j] = data

# Save in NetCDF format.
data_vars = {
        'afrac_dust': (['i', 'j'], arrays['afrac_dust']),
        'afrac_soot': (['i', 'j'], arrays['afrac_soot']),
        'afrac_orga': (['i', 'j'], arrays['afrac_orga'])}

# Define coordinates.
coords = {
        'i': (['i'], np.arange(ni)),
        'j': (['j'], np.arange(nj))}

# Define global attributes.
attrs = {'History': 'Converted from `phillips_nucleation_2010.incf` using `phillips_nucleation_to_nc.py` in microhh_root/misc'}

# Create dataset.
ds = xr.Dataset(
        data_vars = data_vars, 
        coords = coords, 
        attrs = attrs)

# Save with compression.
comp = dict(zlib=True, complevel=5)
encoding = {var: comp for var in ds.data_vars}

ds.to_netcdf('phillips_nucleation_2010.nc', encoding=encoding)
