#######################################
# script to compare microhh solutions #
# for Prandtl slope case with         #
# analytical solutions                #
#######################################
import numpy as np
import netCDF4 as nc
from pylab import *

rcParams['text.usetex'] = True

# read statistics file
stats = nc.Dataset('prandtlslope.default.0000000.nc')

# grab variables from statistics file
ucom = stats.variables['u'][:]
buoy = stats.variables['b'][:]
outt = stats.variables['t'][:]
zlev = stats.variables['z'][:]

# find indices for 40<time<80
avgind = np.where(outt>=40.0)[0]

# calculate means over 40<time<80
ucom_mean = np.mean(ucom[avgind,:],axis=0)
buoy_mean = np.mean(buoy[avgind,:],axis=0)

# constants
N     = 1.
Bs    = 0.005
nu    = 5.E-4
alpha = np.pi/6.

# scaling values
L  = nu**(0.5) * N**(-0.5) * np.sin(alpha)**(-0.5)
B  = Bs * L / nu
V  = B / N
zh = zlev / L
zh_ana = np.arange(0,15.1,0.1)
uh = ucom_mean / V
bh = buoy_mean / B

# analytical values
uh_ana = np.sqrt(2) * np.sin(zh_ana/np.sqrt(2)) * np.exp(-zh_ana/np.sqrt(2))
bh_ana = np.sqrt(2) * np.cos(zh_ana/np.sqrt(2)) * np.exp(-zh_ana/np.sqrt(2))

# plot analytical & numerical
figure(figsize=(12,6))
ax1 = subplot(121)
ax2 = subplot(122)

ax1.set_xlabel(r'$\mathit{u_n}$',fontsize=30)
ax1.set_ylabel(r'$\mathit{z_n}$',fontsize=30)
ax1.axvline(x=0, color='k', ls=':')
ax1.plot(uh_ana, zh_ana, color='k',ls='-',lw=3,label='Analytical')
ax1.plot(uh, zh, color='gray', ls='--',lw=3,label='Numerical')
ax1.set_ylim(0,15)
ax1.set_yticks(np.arange(0,16,3))
for tick in ax1.xaxis.get_major_ticks():
	tick.label.set_fontsize(21) 
for tick in ax1.yaxis.get_major_ticks():
	tick.label.set_fontsize(21) 
ax1.legend()

ax2.set_xlabel(r'$\mathit{b_n}$',fontsize=30)
ax2.set_ylabel(r'$\mathit{z_n}$',fontsize=30)
ax2.axvline(x=0, color='k', ls=':')
ax2.plot(bh_ana, zh_ana, color='k',ls='-',lw=3,label='Analytical')
ax2.plot(bh, zh, color='gray', ls='--',lw=3,label='Numerical')
ax2.set_ylim(0,15)
ax2.set_yticks(np.arange(0,16,3))
ax2.set_xticks(np.arange(-0.4,1.61,0.4))
for tick in ax2.xaxis.get_major_ticks():
	tick.label.set_fontsize(21)
for tick in ax2.yaxis.get_major_ticks():
	tick.label.set_fontsize(21) 
ax2.legend()

subplots_adjust(left=0.07, right=0.98, top=0.97, bottom = 0.13,wspace=0.25)
savefig('prandtlslope.png')
