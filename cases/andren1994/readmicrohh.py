"""
Class to read in the entire profile file from MicroHH
Examples:
>> r=read_microhh('case.nc') #= entire file
>> r=read_microhh('case.nc', t0=7200, t1=10800) #=read only hours 2-3
>> r=read_microhh('case.nc', average=True, t0=7200, t1=10800) #= average profs over hours 2-3

with option to write averaged (1D) data back to NetCDF:
>> r=read_microhh('case.nc', average=True, t0=7200, t1=10800)
>> r.writeback()
This will create a file "case.nc.averaged.nc"
"""

import numpy as np
from netCDF4 import Dataset,default_fillvals
from time import gmtime, strftime

class read_microhh:
    def __init__(self, ncfile, t0=-1, t1=-1, average=False):
        self.nc = Dataset(ncfile)

        #if(t0 != -1 and t1 != -1):
        # find times closest to t0,t1:
        self.t = self.nc.variables['time'][:]

        t0 = (np.abs(self.t-t0)).argmin()   if t0 != -1 else 0
        t1 = (np.abs(self.t-t1)).argmin()+1 if t1 != -1 else self.t.size
        print('%s, av=%i, t=%.1f-%.1f s'%(ncfile,average,self.t[t0],self.t[t1-1]))

        self.z  = self.nc.variables['z'][:]
        self.zh = self.nc.variables['zh'][:]

        # Read in all variables from NetCDF:
        for var in self.nc['default'].variables:
            tmp = self.nc['default'].variables[var]
            shp = tmp.shape

            if(np.size(shp) == 1):
                if(average and tmp.dimensions[0] == 'time'):
                    setattr(self,var,np.mean(tmp[t0:t1]))
                else:
                    setattr(self,var,tmp[:])
            elif(np.size(shp) == 2):
                if(average):
                    setattr(self,var,np.mean(tmp[t0:t1,:],axis=0))
                else:
                    setattr(self,var,tmp[t0:t1,:])

        # Calculate net terms budget (if available):
        try:
            self.u2_net  = self.u2_shear  + self.u2_turb  + self.u2_visc  + self.u2_diss  + self.u2_rdstr 
            self.v2_net  = self.v2_shear  + self.v2_turb  + self.v2_visc  + self.v2_diss  + self.v2_rdstr 
            self.w2_net  =                  self.w2_turb  + self.w2_visc  + self.w2_diss  + self.w2_rdstr + self.w2_pres  
            self.tke_net = self.tke_shear + self.tke_turb + self.tke_visc + self.tke_diss                 + self.tke_pres  
        except:
            self.u2_net  = -1
            self.v2_net  = -1
            self.w2_net  = -1
            self.tke_net = -1

        # Store some input parameters:
        self.ncfile = ncfile
        self.average = average

    # Option to write back the averaged fields to NetCDF
    def writeback(self):
        if(self.average):
            nc2 = Dataset(self.ncfile+'.averaged.nc','w')

            # Copy back dimensions
            for dim in self.nc.dimensions:
                if(dim != 't'):
                    nc2.createDimension(dim,len(self.nc.dimensions[dim]))

            # Copy back averaged data 
            for var in self.nc.variables:
                dim = self.nc.variables[var].dimensions[-1]
                if(dim != 't'):
                    nc2.createVariable(var,'f8',(dim,),fill_value=default_fillvals['f8'])
                    nc2.variables[var].units = self.nc.variables[var].units
                    # TO-DO MicroHH: fix long_name / longname problem 
                    #nc2.variables[var].long_name = self.nc.variables[var].long_name
                    nc2.variables[var][:] = getattr(self,var)
            nc2.close()


#c  = read_microhh('andren1994.default.0000000.nc', t0=70000, t1=100000, average=True)

