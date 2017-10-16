import numpy as np
import matplotlib.pylab as pl
from scipy.special import erf


def esat(T):
    return 0.611e3 * np.exp(17.2694 * (T - 273.16) / (T - 35.86))

def qsat(T,p):
    return 0.622 * esat(T) / p

# Create stretched grid
class Stretched_grid:
    def __init__(self, kmax, nloc1, nbuf1, dz1, dz2):
        dn         = 1./kmax
        n          = np.linspace(dn, 1.-dn, kmax)
        nloc1     *= dn
        nbuf1     *= dn
        dzdn1      = dz1/dn
        dzdn2      = dz2/dn

        dzdn       = dzdn1 + 0.5*(dzdn2-dzdn1)*(1. + np.tanh((n-nloc1)/nbuf1))
        self.dz    = dzdn*dn

        self.kmax  = kmax
        self.z     = np.zeros(self.dz.size)
        stretch    = np.zeros(self.dz.size)

        self.z[0]  = 0.5*self.dz[0]
        stretch[0] = 1.

        for k in range(1, self.kmax):
              self.z [k] = self.z[k-1] + 0.5*(self.dz[k-1]+self.dz[k])
              stretch[k] = self.dz[k]/self.dz[k-1]

        self.zsize = self.z[kmax-1] + 0.5*self.dz[kmax-1]
        print('Grid: kmax=%i, zsize=%f'%(kmax,self.zsize))

    def plot(self):
        pl.figure()
        pl.plot(self.dz, self.z, '-x')
        pl.xlabel('dz (m)')
        pl.ylabel('z (m)')


if __name__ == '__main__':

    import microhh_tools as mht

    # Create stretched grid
    #gr = Stretched_grid(96, 40, 10, 2.0, 4.0)
    #gr = Stretched_grid(64, 40, 10, 5.0, 10.0)
    #gr = Stretched_grid(72,  50, 10, 3.0, 20.0)
    #gr = Stretched_grid(144, 120, 20, 1.0, 15.0)
    gr = Stretched_grid(192, 180, 10, 0.5, 15.0)
    gr.plot()

    # Write back vertical extent domain
    mht.replace_namelist_value('zsize', gr.zsize)
    mht.replace_namelist_value('ktot', gr.z.size)

    # Create initial vertical profiles
    thl = np.ones(gr.z.size)*318
    qt  = 6.4 - gr.z * 0.0030
    qt /= 1000

    u   = np.ones(gr.z.size)*6.25
    ug  = np.ones(gr.z.size)*3

    v   = np.ones(gr.z.size)*3.6
    vg  = np.ones(gr.z.size)*7

    # Large scale tendencies thl and qt
    thlls = -1.395 + 0.00089 * gr.z
    qtls  = 0.0015 + 0.0002  * gr.z

    # To MicroHH units (hours -> sec, g kg-1 -> kg kg-1
    thlls /= 3600.
    qtls  /= (3600. * 1000.)

    # Remove tendency in/near SBL
    zt   = 100       # Height of transition
    dz   = 100       # Depth of transition layer
    fac  = 0.5 + 0.5*erf((gr.z-zt)/(0.25*dz))

    thlls2 = thlls #* fac
    qtls2  = qtls  #* fac

    # Write to .prof file as input for MicroHH
    variables = {'z':gr.z, 'thl':thl, 'qt':qt, 'u':u, 'ug':ug, 'v':v, 'vg':vg, 'thlls':thlls2, 'qtls':qtls2}
    mht.write_output('must.prof', variables, gr.z.size)

    # Create surface time dependent input
    time = np.arange(0.5,5.7501,0.25)       # h UTC
    # Surface potential temperature fitted to observations
    ths_func = np.array([ -2.89610375e-01,   3.46630020e+00,  -1.46448338e+01, 3.25322260e+02])
    th1cm_fit = ths_func[3] + ths_func[2]*time + ths_func[1]*time**2 + ths_func[0]*time**3
    time -= 0.5                             # h since start of experiment
    time *= 3600                            # s   "    "           "

    pbot = 87000
    exnbot = (pbot / 1e5)**(287.04/1005.)
    Tbot = th1cm_fit * exnbot
    qsbot = qsat(Tbot, pbot)

    RHbot = (0.18 + 0.07 * (time/time.max()))
    qbot = RHbot * qsbot

    variables = {'t':time, 'sbot[thl]':th1cm_fit, 'sbot[qt]':qbot}
    mht.write_output('must.time', variables, time.size)
