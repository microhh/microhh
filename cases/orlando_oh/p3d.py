from netCDF4 import *
import pylab as pl
with Dataset("oh.nc") as f:
    x = f.variables['x'][:]
    y = f.variables['y'][:]
    z = f.variables['z'][:]
    oh = f.variables['oh'][:]




