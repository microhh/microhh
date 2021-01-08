import matplotlib.pyplot as pl
import numpy as np

class Read_afgl:
    def __init__(self, afgl_sounding):
        Mw = 18.016
        Md = 28.966

        snd = np.loadtxt(afgl_sounding)

        self.z   = snd[:,0] * 1000.
        self.p   = snd[:,1] * 100
        self.T   = snd[:,2]
        self.dn  = snd[:,3]
        self.h2o = snd[:,4]
        self.o3  = snd[:,5]
        self.n2o = snd[:,6]
        self.co  = snd[:,7]
        self.ch4 = snd[:,8]

        self.thl = self.T / (self.p / 1e5)**(287.04/1005.)
        self.qt  = self.h2o * Mw/Md / 1.e6

        es = 0.611e3 * np.exp(17.2694 * (self.T - 273.16) / (self.T - 35.86))
        qs = 0.622 * es / (self.p * 100.)
        self.RH = self.qt/qs


class Grid:
    def __init__(self, ktot, dz0):
        self.ktot = ktot
        self.dz0  = dz0

        self.z = np.zeros(ktot)
        self.dz = np.zeros(ktot)
        self.zsize = None

    def plot(self):
        pl.figure()
        pl.title(r'$z_\mathrm{{size}}$ = {0:.1f} m'.format(self.zsize), loc='left')
        pl.plot(self.dz, self.z, '-x')
        pl.xlabel(r'$\Delta z$ (m)')
        pl.ylabel(r'$z$ (m)')


class Grid_linear_stretched(Grid):
    def __init__(self, ktot, dz0, alpha):
        Grid.__init__(self, ktot, dz0)

        self.dz[:]  = dz0 * (1 + alpha)**np.arange(ktot)
        self.zh     = np.zeros(ktot+1)
        self.zh[1:] = np.cumsum(self.dz)
        self.z[:]   = 0.5 * (self.zh[1:] + self.zh[:-1])
        self.zsize  = self.zh[-1]
