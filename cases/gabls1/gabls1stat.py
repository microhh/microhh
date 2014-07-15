## Plots all the statistics/figures from Beare et al. 2006 (BLM)

import numpy as np
from pylab import *
from read_microhh import * 

close('all')

# Profiles potential temp, wind speed (u)
# heat and momentum flux
# Read averaged MicroHH profiles:
g128  = read_microhh('gabls1.default.nc',  t0=28800, t1=32400, average=True)

figure()

ax=subplot(221)
plot(g128.th,g128.z)
ylim(0,250)
xlabel('theta [K]')
ylabel('z [m]')

ax=subplot(222)
plot(g128.u,g128.z)
ylim(0,250)
xlabel('u [m/s]')
ylabel('z [m]')

ax=subplot(223)
plot(g128.thflux,g128.zh)
ylim(0,250)
xlabel('wtheta [Km/s]')
ylabel('z [m]')
    
ax=subplot(224)
plot(g128.uflux,g128.zh)
ylim(0,250)
xlabel('uw [m2/s2]')
ylabel('z [m]')

# Time series
# Read all MicroHH data:
g128  = read_microhh('gabls1.default.nc')

figure()
    
ax=subplot(131)
plot(g128.t/3600.,g128.thflux[:,0])
xlabel('time [h]')
ylabel('wthetas [Kms-1]')

ax=subplot(132)
plot(g128.t/3600.,g128.ustar[:])
ylim(0,0.5)
xlabel('time [h]')
ylabel('u* [ms-1]')

ax=subplot(133)
plot(g128.t/3600.,g128.obuk[:])
ylim(0,300) 
xlabel('time [h]')
ylabel('Obukhov length [m]')


