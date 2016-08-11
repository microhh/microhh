import numpy as np
import matplotlib.pylab as pl
import netCDF4 as nc4

pl.close('all')

class Read_statistics:
    def __init__(self, stat_file):
        f = nc4.Dataset(stat_file, 'r')
       
        # Dictionaries holding the units and variable names
        self.units = {}
        self.name  = {}

        for var in f.variables:
            setattr(self, var, f.variables[var].__array__())
            self.units[var] = f.variables[var].getncattr('units')
            self.name[var]  = f.variables[var].getncattr('long_name')
        f.close()

sc = Read_statistics('rain.patch_high.0000000.nc')
sr = Read_statistics('rain.patch_low.0000000.nc')



pl.figure()
pl.subplot(221)
pl.plot(sc.thl[-1,:], sc.z, label='city')
pl.plot(sr.thl[-1,:], sr.z, label='rural')
pl.legend()

pl.subplot(222)
pl.plot(sc.qt[-1,:], sc.z, label='city')
pl.plot(sr.qt[-1,:], sr.z, label='rural')

pl.subplot(223)
pl.plot(sc.ql[-1,:], sc.z, label='city')
pl.plot(sr.ql[-1,:], sr.z, label='rural')

pl.subplot(224)
pl.plot(sc.qr[-1,:], sc.z, label='city')
pl.plot(sr.qr[-1,:], sr.z, label='rural')
