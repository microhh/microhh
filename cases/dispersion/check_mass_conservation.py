import numpy as np
import xarray as xr

import microhh_tools as mht

def xr_read_all(f, groups=['default'], decode_times=True):
    # Read all NetCDF groups into a single Dataset.
    dss = [xr.open_dataset(f, decode_times=decode_times)]
    for group in groups:
        dss.append(xr.open_dataset(f, group=group, decode_times=decode_times))
    return xr.merge(dss)

dtype = np.float32

ini = mht.Read_namelist('dispersion.ini')

stats = xr_read_all('dispersion.default.0000000.nc', groups=['default', 'thermo'], decode_times=False)
rhoref = stats.rhoref.values
zh = stats.zh.values
dz = zh[1:] - zh[:-1]

# Short-cuts.
itot = ini['grid']['itot']
jtot = ini['grid']['jtot']
ktot = ini['grid']['ktot']

xsize = ini['grid']['xsize']
ysize = ini['grid']['ysize']

dx = xsize / itot
dy = ysize / jtot

# Use `s2`, which used period BCs, so no mass is lost.
scalars = ini['source']['sourcelist']
strength = np.array(ini['source']['strength'])

for time in range(0, ini['time']['endtime']+1, ini['time']['savetime']):
    correct_mass = strength * time

    for i,scalar in enumerate(scalars):
        fld = np.fromfile('{0}.{1:07d}'.format(scalar, time), dtype)
        fld = fld.reshape((ktot, jtot, itot))
        mass = np.sum(rhoref[:,None,None] * fld * dx * dy * dz[:,None,None])

        print('Time = {0:5d}, scalar \"{1}\": correct mass = {2:8.2f} kg, integral field = {3:8.2f} kg'.format(time, scalar, correct_mass[i], mass))
