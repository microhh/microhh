# Class to read in the entire profile file from MicroHH
# Examples:
# read_microhh('case.nc') #= entire file
# read_microhh('case.nc',t0=7200,t1=10800) #=read only hours 2-3
# read_microhh('case.nc',average=True,t0=7200,t1=10800) #= average profs over hours 2-3

import numpy as np
from netCDF4 import Dataset

class read_microhh:
    def __init__(self, ncfile, t0=-1, t1=-1, average=False):
        nc = Dataset(ncfile)

        # find times closest to t0,t1:
        self.t = nc.variables['t'][:]  
        t0 = (np.abs(self.t-t0)).argmin()   if t0 != -1 else 0
        t1 = (np.abs(self.t-t1)).argmin()+1 if t1 != -1 else self.t.size
        print('%s, av=%i, t=%.1f-%.1f s'%(ncfile,average,self.t[t0],self.t[t1-1]))

        # Read in all variables from NetCDF:
        for var in nc.variables:
            tmp = nc.variables[var]
            shp = tmp.shape

            if(np.size(shp) == 1):
                setattr(self,var,tmp[:])
            elif(np.size(shp) == 2):
                if(average):
                    setattr(self,var,np.mean(tmp[t0:t1,:],axis=0))
                else:
                    setattr(self,var,tmp[t0:t1,:])

        self.t = nc.variables['t'][t0:t1]   

