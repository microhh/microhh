import matplotlib.pyplot as pl
import numpy as np
from readmicrohh import *

pl.close('all')

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

if(True):
    # Figure 3a & 3b

    tf = ct.t * fc;
    Cu = np.zeros(tf.size)
    Cv = np.zeros(tf.size)
    for t in range(tf.size):
        for k in range(c.z.size):
            Cu[t] += (ct.v[t,k] - vg) * (c.zh[k+1] - c.zh[k])
            Cv[t] += (ct.u[t,k] - ug) * (c.zh[k+1] - c.zh[k])
        Cu[t] *= -fc / ct.u_flux[t,0]
        Cv[t] *=  fc / ct.v_flux[t,0]
    
    pl.figure()
    pl.subplot(121)
    pl.title('Fig 3a')
    pl.plot(tf,Cu)
    pl.xlabel('t*fc')
    pl.ylabel('Cu')

    pl.subplot(122)
    pl.title('Fig 3b')
    pl.plot(tf,Cv)
    pl.xlabel('t*fc')
    pl.ylabel('Cv')


if(True):
    # Figure 6
    zp  = c.z  * fc / c.ustar
    zhp = c.zh * fc / c.ustar

    pl.figure()
    pl.subplot(121)
    pl.title('Fig 6a')
    pl.plot(c.u_flux         / c.ustar**2., zhp, label='total')
    pl.plot(c.u_w            / c.ustar**2., zhp, label='resolved')
    pl.plot((c.u_flux-c.u_w) / c.ustar**2., zhp, label='sgs')
    pl.legend(frameon=False)
    pl.xlabel('uw/u*^2')
    pl.ylabel('z*fc/u*')

    pl.subplot(122)
    pl.title('Fig 6b')
    pl.plot(c.v_flux         / c.ustar**2., zhp, label='total')
    pl.plot(c.v_w            / c.ustar**2., zhp, label='resolved')
    pl.plot((c.v_flux-c.v_w) / c.ustar**2., zhp, label='sgs')
    pl.legend(frameon=False)
    pl.xlabel('vw/u*^2')
    pl.ylabel('z*fc/u*')


if(True):
    kappa = 0.4
    phi = kappa * c.zh / c.ustar * ((c.u_grad)**2 + (c.v_grad)**2)**.5

    pl.figure()
    pl.plot(phi, zhp)
    pl.legend(frameon=False)
    pl.xlabel('phi_m')
    pl.ylabel('z*fc/u*')
