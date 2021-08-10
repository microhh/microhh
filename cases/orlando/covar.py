from netCDF4 import *
import numpy as np
import pylab as pl
for i in range(4,11):
    with Dataset("oh.nc") as f:
        oh = f.variables['oh'][i]
    with Dataset("isop.nc") as f:
        isop = f.variables['isop'][i]
    k = 5
    f,ax = pl.subplots(2,sharex=True)
    ax[0].imshow(isop[k,:,:],origin='lower')
    ax[1].imshow(oh[k,:,:],origin='lower')
    cv = []
    for k in range(24):
      a = isop[k,:,:].reshape(24*24)
      b = oh[k,:,:].reshape(24*24)
      cv.append(np.cov(a,b)[0,1]*100/(a.mean()*b.mean()))
    f,ax = pl.subplots()
    ax.plot(cv,range(24))






