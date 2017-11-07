import numpy as np
import matplotlib.pylab as pl
import netCDF4 as nc4

pl.close('all')

class Read_stat:
    def __init__(self, file_name):
        f = nc4.Dataset(file_name)
        for v in f.variables:
            setattr(self, v, f.variables[v][:])

s = Read_stat('ls_force.default.0000000.nc')

pl.figure()
pl.subplot(231)
pl.plot(s.t, s.phh[:,0])
pl.ylabel('Surface pressure (Pa)')
pl.xlabel('time (s)')

pl.subplot(232)
pl.plot(s.t, s.thlflux[:,0])
pl.ylabel('wthl_s (K m s-1)')
pl.xlabel('time (s)')

pl.subplot(233)
pl.plot(s.thl[0,:], s.z, label='t={}'.format(s.t[0]))
pl.plot(s.thl[18,:], s.z, label='t={}'.format(s.t[18]), dashes=[4,2])
pl.plot(s.thl[-1,:], s.z, label='t={}'.format(s.t[-1]))
pl.legend()
pl.xlabel('thl (K)')
pl.ylabel('z (m)')

pl.subplot(234)
pl.plot(s.qt[0,:], s.z)
pl.plot(s.qt[18,:], s.z, dashes=[4,2])
pl.plot(s.qt[-1,:], s.z)
pl.xlabel('qt (kg kg-1)')
pl.ylabel('z (m)')

pl.subplot(235)
pl.plot(s.s[0,:], s.z)
pl.plot(s.s[18,:], s.z, dashes=[4,2])
pl.plot(s.s[-1,:], s.z)
pl.xlabel('scalar (-)')
pl.ylabel('z (m)')

pl.subplot(236)
pl.plot(s.u[0,:], s.z)
pl.plot(s.u[18,:], s.z, dashes=[4,2])
pl.plot(s.u[-1,:], s.z)
pl.xlabel('u (m s-1)')
pl.ylabel('z (m)')

pl.tight_layout()
