import numpy as np
import matplotlib.pylab as pl
from scipy.special import erf

from microhh_tools import Stretched_grid, write_output, replace_namelist_value

def esat(T):
    return 0.611e3 * np.exp(17.2694 * (T - 273.16) / (T - 35.86))

def qsat(T,p):
    return 0.622 * esat(T) / p

# Create stretched grid
#gr = Stretched_grid(96, 40, 10, 2.0, 4.0)
gr = Stretched_grid(64, 40, 10, 5.0, 10.0)

#gr.plot()

# Write back vertical extent domain
replace_namelist_value('zsize', gr.zsize)
replace_namelist_value('ktot', gr.z.size)

# Create initial vertical profiles
thl = np.ones(gr.z.size)*318
qt  = 6.4 - gr.z * 0.0045
qt /= 1000.

u   = np.ones(gr.z.size)*5
ug  = np.zeros(gr.z.size)

v   = np.ones(gr.z.size)*3
vg  = np.zeros(gr.z.size)

# Large scale tendencies thl and qt
#thlls = -1.395 + 0.00069 * gr.z
#qtls  = 0.0015 + 0.0004  * gr.z
thlls = -1.395 + 0.00089 * gr.z
qtls  = 0.0015 + 0.0005  * gr.z

# To MicroHH units (hours -> sec, g kg-1 -> kg kg-1
thlls /= 3600.
qtls  /= (3600. * 1000.)

# Remove tendency in/near SBL
zt   = 100       # Height of transition
dz   = 100       # Depth of transition layer
fac  = 0.5 + 0.5*erf((gr.z-zt)/(0.25*dz))

thlls2 = thlls * fac
qtls2  = qtls  * fac

# Write to .prof file as input for MicroHH
variables = {'z':gr.z, 'thl':thl, 'qt':qt, 'u':u, 'ug':ug, 'v':v, 'vg':vg, 'thlls':thlls2, 'qtls':qtls2}
write_output('must.prof', variables, gr.z.size)

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
qbot = 1.0 * qsbot

variables = {'t':time, 'sbot[thl]':th1cm_fit, 'sbot[qt]':qbot}
write_output('must.time', variables, time.size)






