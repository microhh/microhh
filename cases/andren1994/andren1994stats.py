import numpy as np
from pylab import *
from readmicrohh import *

close('all')

fc    = 0.0001
kappa = 0.4
ug    = 10.
vg    = 0.

# averaging period profiles:
t0 = 7 * fc**-1
t1 = 10* fc**-1

# Read non-averaged statistics
ct = read_microhh('andren1994.default.0000000.nc')
# Read statistics averaged over 7-10 * fc
c  = read_microhh('andren1994.default.0000000.nc', t0=t0, t1=t1, average=True)

# Figure 3a & 3b
# ---------------------------------
if(True):
    tf = ct.t * fc;
    Cu = np.zeros(tf.size)
    Cv = np.zeros(tf.size)
    for t in range(tf.size):
        for k in range(c.z.size):
            Cu[t] += (ct.v[t,k] - vg) * (c.zh[k+1] - c.zh[k])
            Cv[t] += (ct.u[t,k] - ug) * (c.zh[k+1] - c.zh[k])
        Cu[t] *= -fc / ct.uflux[t,0]
        Cv[t] *=  fc / ct.vflux[t,0]
    
    figure()
    subplot(121)
    title('Fig 3a')
    plot(tf,Cu)
    xlabel('t*fc')
    ylabel('Cu')

    subplot(122)
    title('Fig 3b')
    plot(tf,Cv)
    xlabel('t*fc')
    ylabel('Cv')

# Figure 6
# ---------------------------------
if(True):
    zp  = c.z * fc / c.ustar
    zhp = c.zh * fc / c.ustar

    figure()
    subplot(121)
    title('Fig 6a')
    plot(c.uflux        /c.ustar**2.,zhp, label='total')
    plot(c.uw           /c.ustar**2.,zhp, label='resolved')
    plot((c.uflux-c.uw) /c.ustar**2.,zhp, label='sgs')
    legend(frameon=False)
    xlabel('uw/u*^2')
    ylabel('z*fc/u*')

    subplot(122)
    title('Fig 6b')
    plot(c.vflux        /c.ustar**2.,zhp, label='total')
    plot(c.vw           /c.ustar**2.,zhp, label='resolved')
    plot((c.vflux-c.vw) /c.ustar**2.,zhp, label='sgs')
    legend(frameon=False)
    xlabel('vw/u*^2')
    ylabel('z*fc/u*')

if(True):
    kappa = 0.4
    phi = kappa * c.zh / c.ustar * c.ugrad

    figure()
    plot(phi, zhp)
    legend(frameon=False)
    xlabel('phi_m')
    ylabel('z*fc/u*')
