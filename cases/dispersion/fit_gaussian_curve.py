import matplotlib.pyplot as pl
import matplotlib.colors as colors
from scipy.optimize import curve_fit
import xarray as xr
import numpy as np

pl.close('all')

def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
     
def gauss_fit(x, y):
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma])
    return popt,pcov

ds = xr.open_dataset('s1_path.xy.nc')
dsm = ds.sel(time=slice(7200, 43200)).mean(dim='time')

x = dsm.x.values
y = dsm.y.values

# Fit Gaussian curve at some distances from the source:
x_trans = np.arange(2500, 12800, 2000)
N = x_trans.size
colors = pl.cm.plasma(np.linspace(0, 1, N))

pl.figure(figsize=(6,4))
for i in range(N):
    ii = np.argmin(np.abs(x - x_trans[i]))
    s = dsm.s1_path[:,ii].values

    popt, pcov = gauss_fit(y, s)
    pl.plot(y/1000, s, color=colors[i], label='LES, x={} km'.format(x_trans[i]/1000))
    pl.plot(y/1000, gauss(y, *popt), ':', color=colors[i], label='Fit, sigma={0:.3f} km'.format(popt[-1]/1000))

pl.legend(ncol=1, fontsize=9)
pl.xlabel('y (km)')
pl.ylabel(r'vertical integral s (kg m$^{-2}$)')
pl.tight_layout()
pl.savefig('gaussian_fit.png')
