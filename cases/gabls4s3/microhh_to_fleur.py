# 
# Script to convert the output from MicroHH to the GABLS4 LES specifications
#

import numpy as np
from netCDF4 import Dataset
import struct  as st
from pylab import * ## TMP

# load the init script to get variables like ug, vg, Ts
from gabls4s3init import *

# Settings NetCDF
default_fill_value = -9999

# Settings case
institute  = 'mp'
model      = 'MicroHH'
version    = 'v01'
experiment = 'exp1'

# Constants
rho = 1.     # Boussinesq
cp  = 1005. 
Lv  = 2.5e6

# Case specific settings
ps  = 65100
z0m = 0.01
z0h = 0.001

def key_nearest(array,value):
    return (np.abs(array-value)).argmin()

def add_global_attr(nc):
    nc.setncattr('Model', 'MicroHH (github.com/microhh)')
    nc.setncattr('Reference', 'van Heerwaarden et al, in preparation for Geosci. Model Dev.')
    nc.setncattr('Contact', 'Bart van Stratum (bart.vanstratum@mpimet.mpg.de), Chiel van Heerwaarden (chiel.vanheerwaarden@mpimet.mpg.de)')
    nc.setncattr('Diffusion', 'Smagorinsky-Lilly, Cs=0.1, tPr=1/3, L=(dxdydz)^1/3, with wall-damping (Mason and Thomson, 1992) with n=2')
    nc.setncattr('Advection', 'momentum and scalar: 2nd order centered with 4th order interpolations')
    nc.setncattr('Time step', 'variable, to maintain CFL=1.2')

def interpol(field, z, height):
    k = key_nearest(z, height)
    if(z[k] > height): k -= 1

    frac = (height - z[k]) / (z[k+1] - z[k])

    return (1.-frac) * field[:,k] + frac * field[:,k+1]

def interpol_tke(u2, v2, w2, zf, zh, height):
    kf  = key_nearest(zf, height)
    if(zf[kf] > height): kf -= 1
    kh = key_nearest(zh, height)
    if(zh[kh] > height): kh -= 1

    fracf = (height - zf[kf]) / (zf[kf+1] - zf[kf])
    frach = (height - zh[kh]) / (zh[kh+1] - zh[kh])

    u2i = (1.-fracf) * u2[:,kf] + fracf * u2[:,kf+1]
    v2i = (1.-fracf) * v2[:,kf] + fracf * v2[:,kf+1]
    w2i = (1.-frach) * w2[:,kh] + frach * w2[:,kh+1]

    return 0.5 * (u2i + v2i + w2i)

class read_stat:
    def __init__(self, ncfile):
        self.nc      = Dataset(ncfile)

        # Read in all variables from NetCDF:
        for var in self.nc.variables:
            tmp = self.nc.variables[var]
            shp = tmp.shape

            if(np.size(shp) == 1):
                setattr(self,var,tmp[:])
            elif(np.size(shp) == 2):
                setattr(self,var,tmp[:,:])

    def write_time_series(self):
        # Read list with required variables:
        var_file  = np.genfromtxt('variables_ts.txt', delimiter=',', dtype='str', autostrip=True)
        var_name  = var_file[:,0]
        var_lname = var_file[:,1]
        var_unit  = var_file[:,2]

        # Create new NetCDF file
        nc_ts = Dataset('gabls4_time_les_%s_%s_%s_stage3_%s.nc'%(institute, model, experiment, version), 'w')
      
        # Add global attributes:
        add_global_attr(nc_ts)
       
        # Create dimension:
        dim_t = nc_ts.createDimension('time', self.t.size)
  
        # Create variables:
        self.time_ts = nc_ts.createVariable('time', 'f4', ('time'), fill_value=default_fill_value)
        self.time_ts.setncattr('units', 's')
        self.time_ts.setncattr('long_name', 'seconds since start of experiment')

        for i in range(np.size(var_name)):
            tmp = nc_ts.createVariable(var_name[i], 'f4', 'time', fill_value=default_fill_value)
            tmp[:] = default_fill_value
            tmp.setncattr('units', var_unit[i])
            tmp.setncattr('long_name', var_lname[i])
            setattr(self, var_name[i]+'_ts', tmp)

        # Fill them (only the ones which are available, other get fill_value automatically):
        self.time_ts[:]    = self.t
        self.shf_ts[:]     = self.thflux[:,0] * rho * cp
        self.ustar_ts[:]   = self.ustar
        self.psurf_ts[:]   = ps
        self.thsurf_ts[:]  = np.interp(self.t, s3.t, s3.ths) # interpolated from model input 

        # Calculate ABL depth based on bulk Richardson number (no moisture):
        figure()
        for t in range(self.t.size):
            Rib = 9.81 * self.z / self.thsurf_ts[t] * (self.th[t,:] - self.thsurf_ts[t]) / (self.u[t,:]**2. + self.v[t,:]**2.)
            kzi = np.where(Rib > 0.25)[0][0]
            self.hpbl_ts[t] = self.z[kzi]

        self.z0m_ts[:]     = z0m
        self.z0h_ts[:]     = z0h

        # Quantities interpolated to SCM grid
        self.th2m_ts[:]    = interpol(self.th,     self.z, 2.0 )
        self.th3m_ts[:]    = interpol(self.th,     self.z, 3.3 )
        self.th9m_ts[:]    = interpol(self.th,     self.z, 8.8 )
        self.th18m_ts[:]   = interpol(self.th,     self.z, 17.9)
        self.th25m_ts[:]   = interpol(self.th,     self.z, 25.3)
        self.th33m_ts[:]   = interpol(self.th,     self.z, 32.7)
        self.th42m_ts[:]   = interpol(self.th,     self.z, 41.9)

        self.u10m_ts[:]    = interpol(self.u,      self.z, 10. )
        self.u3m_ts[:]     = interpol(self.u,      self.z, 3.3 )
        self.u9m_ts[:]     = interpol(self.u,      self.z, 8.8 )
        self.u18m_ts[:]    = interpol(self.u,      self.z, 17.9)
        self.u25m_ts[:]    = interpol(self.u,      self.z, 25.3)
        self.u33m_ts[:]    = interpol(self.u,      self.z, 32.7)
        self.u42m_ts[:]    = interpol(self.u,      self.z, 41.9)

        self.v10m_ts[:]    = interpol(self.v,      self.z, 10. )
        self.v3m_ts[:]     = interpol(self.v,      self.z, 3.3 )
        self.v9m_ts[:]     = interpol(self.v,      self.z, 8.8 )
        self.v18m_ts[:]    = interpol(self.v,      self.z, 17.9)
        self.v25m_ts[:]    = interpol(self.v,      self.z, 25.3)
        self.v33m_ts[:]    = interpol(self.v,      self.z, 32.7)
        self.v42m_ts[:]    = interpol(self.v,      self.z, 41.9)

        self.uw_3m_ts[:]   = interpol(self.uflux,  self.zh, 3.3)
        self.vw_3m_ts[:]   = interpol(self.vflux,  self.zh, 3.3)
        self.wth_3m_ts[:]  = interpol(self.thflux, self.zh, 3.3)
        self.TKE_3m_ts[:]  = interpol_tke(self.u2, self.v2, self.w2, self.z, self.zh, 3.3)

        self.uw_7m_ts[:]   = interpol(self.uflux,  self.zh, 7.03)
        self.vw_7m_ts[:]   = interpol(self.vflux,  self.zh, 7.03)
        self.wth_7m_ts[:]  = interpol(self.thflux, self.zh, 7.03)
        self.TKE_7m_ts[:]  = interpol_tke(self.u2, self.v2, self.w2, self.z, self.zh, 7.03)

        self.uw_15m_ts[:]  = interpol(self.uflux,  self.zh, 15.43)
        self.vw_15m_ts[:]  = interpol(self.vflux,  self.zh, 15.43)
        self.wth_15m_ts[:] = interpol(self.thflux, self.zh, 15.43)
        self.TKE_15m_ts[:] = interpol_tke(self.u2, self.v2, self.w2, self.z, self.zh, 15.43)

        self.uw_23m_ts[:]  = interpol(self.uflux,  self.zh, 22.79)
        self.vw_23m_ts[:]  = interpol(self.vflux,  self.zh, 22.79)
        self.wth_23m_ts[:] = interpol(self.thflux, self.zh, 22.79)
        self.TKE_23m_ts[:] = interpol_tke(self.u2, self.v2, self.w2, self.z, self.zh, 22.79)

        self.uw_30m_ts[:]  = interpol(self.uflux,  self.zh, 30.15)
        self.vw_30m_ts[:]  = interpol(self.vflux,  self.zh, 30.15)
        self.wth_30m_ts[:] = interpol(self.thflux, self.zh, 30.15)
        self.TKE_30m_ts[:] = interpol_tke(self.u2, self.v2, self.w2, self.z, self.zh, 30.15)

        self.uw_38m_ts[:]  = interpol(self.uflux,  self.zh, 37.51)
        self.vw_38m_ts[:]  = interpol(self.vflux,  self.zh, 37.51)
        self.wth_38m_ts[:] = interpol(self.thflux, self.zh, 37.51)
        self.TKE_38m_ts[:] = interpol_tke(self.u2, self.v2, self.w2, self.z, self.zh, 37.51)

        nc_ts.close()

    def write_profiles(self, average=False):
        # Analysis times
        times = np.arange(0, self.t.max()+0.1, 600)

        # Read list with required variables:
        var_file  = np.genfromtxt('variables_ps.txt', delimiter=',', dtype='str', autostrip=True)
        var_name  = var_file[:,0]
        var_level = var_file[:,1]
        var_lname = var_file[:,2]
        var_unit  = var_file[:,3]

        # Create new NetCDF file
        if(average):
            nc_ps = Dataset('gabls4_meanprofile_les_%s_%s_%s_stage3_%s.nc'%(institute, model, experiment, version), 'w')
        else:
            nc_ps = Dataset('gabls4_instprofile_les_%s_%s_%s_stage3_%s.nc'%(institute, model, experiment, version), 'w')
      
        # Add global attributes:
        add_global_attr(nc_ps)
       
        ## Create dimension:
        dim_t   = nc_ps.createDimension('time', times.size)
        dim_z   = nc_ps.createDimension('levf', self.z.size)
        dim_zh  = nc_ps.createDimension('levh', self.zh.size)
  
        # Create variables:
        self.time_ps = nc_ps.createVariable('time', 'f4', ('time'), fill_value=default_fill_value)
        self.time_ps.setncattr('units', 's')
        self.time_ps.setncattr('long_name', 'seconds since start of experiment')

        for i in range(np.size(var_name)):
            dim = 'levf' if var_level[i] == 'f' else 'levh'
            tmp = nc_ps.createVariable(var_name[i], 'f4', ('time', dim), fill_value=default_fill_value)
            tmp[:,:] = default_fill_value
            tmp.setncattr('units', var_unit[i])
            tmp.setncattr('long_name', var_lname[i])
            setattr(self, var_name[i]+'_ps', tmp)

        # Calculate w-variance and TKE at full level
        w2_full = 0.5 * (self.w2[:,1:]+self.w2[:,:-1])
        TKE_res = 0.5 * (self.u2[:,:] + self.v2[:,:] + w2_full[:,:])

        # Fill them (only the ones which are available, other get fill_value automatically):
        for t in range(times.size):
            tt = key_nearest(self.t, times[t])
            t0 = key_nearest(self.t, times[t]-300) # 600s averaging window
            t1 = key_nearest(self.t, times[t]+300)

            # Same for averaged and non-averaged profiles
            self.zf_ps  [t,:] = self.z[:]
            self.zh_ps  [t,:] = self.zh[:]
            self.ugeo_ps[t,:] = ug           # from gabls4s3init.py
            self.vgeo_ps[t,:] = vg           # from gabls4s3init.py
            self.time_ps[t]   = self.t[tt]

            # Calculate mean profile, or write instantaneous profile. t=0 is never averaged
            if(average and t!=0):
                self.th_ps     [t,:] = np.apply_over_axes(np.mean, self.th    [t0:t1+1,:], 0)[0,:]
                self.u_ps      [t,:] = np.apply_over_axes(np.mean, self.u     [t0:t1+1,:], 0)[0,:] 
                self.v_ps      [t,:] = np.apply_over_axes(np.mean, self.v     [t0:t1+1,:], 0)[0,:]
                self.wth_sbg_ps[t,:] = np.apply_over_axes(np.mean, self.thdiff[t0:t1+1,:], 0)[0,:]
                self.wth_res_ps[t,:] = np.apply_over_axes(np.mean, self.thflux[t0:t1+1,:], 0)[0,:]
                self.uw_sbg_ps [t,:] = np.apply_over_axes(np.mean, self.udiff [t0:t1+1,:], 0)[0,:]
                self.uw_res_ps [t,:] = np.apply_over_axes(np.mean, self.uflux [t0:t1+1,:], 0)[0,:]
                self.vw_sbg_ps [t,:] = np.apply_over_axes(np.mean, self.vdiff [t0:t1+1,:], 0)[0,:]
                self.vw_res_ps [t,:] = np.apply_over_axes(np.mean, self.vflux [t0:t1+1,:], 0)[0,:]
                self.uu_res_ps [t,:] = np.apply_over_axes(np.mean, self.u2    [t0:t1+1,:], 0)[0,:]
                self.vv_res_ps [t,:] = np.apply_over_axes(np.mean, self.v2    [t0:t1+1,:], 0)[0,:]
                self.ww_res_ps [t,:] = np.apply_over_axes(np.mean, self.w2    [t0:t1+1,:], 0)[0,:]
                self.th2_res_ps[t,:] = np.apply_over_axes(np.mean, self.th2   [t0:t1+1,:], 0)[0,:]
                self.TKE_res_ps[t,:] = np.apply_over_axes(np.mean, TKE_res    [t0:t1+1,:], 0)[0,:]
            else:
                self.th_ps     [t,:] = self.th    [tt,:]
                self.u_ps      [t,:] = self.u     [tt,:] 
                self.v_ps      [t,:] = self.v     [tt,:]
                self.wth_sbg_ps[t,:] = self.thdiff[tt,:]
                self.wth_res_ps[t,:] = self.thflux[tt,:]
                self.uw_sbg_ps [t,:] = self.udiff [tt,:]
                self.uw_res_ps [t,:] = self.uflux [tt,:]
                self.vw_sbg_ps [t,:] = self.vdiff [tt,:]
                self.vw_res_ps [t,:] = self.vflux [tt,:]
                self.uu_res_ps [t,:] = self.u2    [tt,:]
                self.vv_res_ps [t,:] = self.v2    [tt,:]
                self.ww_res_ps [t,:] = self.w2    [tt,:]
                self.th2_res_ps[t,:] = self.th2   [tt,:]
                self.TKE_res_ps[t,:] = TKE_res    [tt,:]

        nc_ps.close()

def convert_3d(times, itot=608, jtot=608, ktot=304, ksave=246):
    # Read grid properties from grid.0000000
    n   = itot*jtot*ktot
    fin = open("grid.{:07d}".format(0),"rb")
    raw = fin.read(itot*8)
    x   = np.array(st.unpack('{0}{1}d'.format('<', itot), raw))
    raw = fin.read(itot*8)
    xh  = np.array(st.unpack('{0}{1}d'.format('<', itot), raw))
    raw = fin.read(jtot*8)
    y   = np.array(st.unpack('{0}{1}d'.format('<', jtot), raw))
    raw = fin.read(jtot*8)
    yh  = np.array(st.unpack('{0}{1}d'.format('<', jtot), raw))
    raw = fin.read(ktot*8)
    z   = np.array(st.unpack('{0}{1}d'.format('<', ktot), raw))
    raw = fin.read(ktot*8)
    zh  = np.array(st.unpack('{0}{1}d'.format('<', ktot), raw))
    fin.close()

    for time in times:
        nc = Dataset('gabls4_3D_les_%02i_%s_%s_%s_stage3_%s.nc'%(time/3600., institute, model, experiment, version), 'w')

        add_global_attr(nc)

        # Create dimensions
        dim_t  = nc.createDimension('time', 1)
        dim_z  = nc.createDimension('lev',  ksave)
        dim_zh = nc.createDimension('levh', ksave)
        dim_y  = nc.createDimension('lat',  jtot)
        dim_yh = nc.createDimension('lath', jtot)
        dim_x  = nc.createDimension('lon',  itot)
        dim_xh = nc.createDimension('lonh', itot)

        # Create grid and time variables:
        var_t  = nc.createVariable('time', 'f4', 'time')
        var_t.setncattr('units', 's')
        var_t.setncattr('long_name', 'seconds since start of experiment')

        var_x  = nc.createVariable('lon',  'f4', 'lon')
        var_x.setncattr('units', 'm')
        var_x.setncattr('long_name', 'x-location cell center')

        var_y  = nc.createVariable('lat',  'f4', 'lat')
        var_y.setncattr('units', 'm')
        var_y.setncattr('long_name', 'y-location cell center')
        
        var_z  = nc.createVariable('lev', 'f4', 'lev')
        var_z.setncattr('units', 'm')
        var_z.setncattr('long_name', 'z-location cell center')
 
        var_xh = nc.createVariable('lonh', 'f4', 'lonh')
        var_xh.setncattr('units', 'm')
        var_xh.setncattr('long_name', 'x-location cell edge')

        var_yh = nc.createVariable('lath', 'f4', 'lath')
        var_yh.setncattr('units', 'm')
        var_yh.setncattr('long_name', 'y-location cell edge')

        var_zh = nc.createVariable('levh', 'f4', 'levh')
        var_zh.setncattr('units', 'm')
        var_zh.setncattr('long_name', 'z-location cell edge')

        # Create 3D variables:
        var_u  = nc.createVariable('u',     'f4', ('time','lev','lat','lonh',))
        var_u.setncattr('units', 'm/s')
        var_u.setncattr('long_name', 'zonal wind')

        var_v  = nc.createVariable('v',     'f4', ('time','lev','lath','lon',))
        var_v.setncattr('units', 'm/s')
        var_v.setncattr('long_name', 'meridional wind')

        var_w  = nc.createVariable('w',     'f4', ('time','levh','lat','lon',))
        var_w.setncattr('units', 'm/s')
        var_w.setncattr('long_name', 'vertical velocity')

        var_th = nc.createVariable('theta', 'f4', ('time','lev','lat','lon',))
        var_w.setncattr('units', 'K')
        var_w.setncattr('long_name', 'potential temperature')

        # Write grid and time variables:
        var_t[:]  = time
        var_x[:]  = x
        var_y[:]  = y
        var_z[:]  = z[:ksave]
        var_zh[:] = zh[:ksave]

        # Write 3D fields:
        for var,field in zip(['u','v','w','th'],[var_u, var_v, var_w, var_th]):
            fin = open("%s.%07i"%(var,time),"rb")
            for k in range(ksave):
                raw = fin.read(itot*jtot*8)
                tmp = np.array(st.unpack('{0}{1}d'.format('<', itot*jtot), raw))
                field[0,k,:,:] = tmp.reshape((jtot, itot))
            fin.close()

        nc.close()

if(True):
    # Convert time series and profile statistics
    r1 = read_stat('gabls4s3.default.nc')
    r1.write_time_series()
    r1.write_profiles(average=False)
    r1.write_profiles(average=True)

# Convert 3D files
if(True):
    times = np.array([5,7,9,11,13,15,17,19,21,23])*3600.
    convert_3d(times)
