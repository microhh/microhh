import numpy as np
import pylab as pl
import netCDF4 as nc4

pl.close('all')

# Read the stage 3 driver file
class read_driver:

    def __init__(self):
        nc = nc4.Dataset('SCM_LES_STAGE3_10h.nc', 'r')
        self.t = nc.variables["time"][:]
        self.z = nc.variables["height"][:][::-1]
        self.p = nc.variables["pf"][:][::-1]

        # initial profiles
        self.th = nc.variables["theta"][:][::-1]
        self.T  = nc.variables["t"][:][::-1]
        self.u  = nc.variables["u"][:][::-1]
        self.v  = nc.variables["v"][:][::-1]

        # time varying forcings
        self.ug   = nc.variables["Ug"][0,:][::-1]    # u geo wind = constant in time
        self.vg   = nc.variables["Vg"][0,:][::-1]    # v geo wind = consant in time
        
        self.Ts   = nc.variables["Tg"][:]
        self.ps  = nc.variables["psurf"].getValue()
        self.z0m = nc.variables["z0m"].getValue()

        # Calculate theta_s
        self.ths = self.Ts / (self.ps / 1.e5)**(287.04/1005.)

        # Extrapolate everything linearly to the surface (z=0)
        self.extrapolate_to_surface("th")
        self.extrapolate_to_surface("u")
        self.extrapolate_to_surface("v")
        self.extrapolate_to_surface("ug")
        self.extrapolate_to_surface("vg")
        self.z = np.insert(self.z, 0, 0)

    def extrapolate_to_surface(self, var):
        arr = getattr(self, var)
        ddz = (arr[1]-arr[0])/(self.z[1]-self.z[0])
        val = arr[0] - ddz*self.z[0]
        setattr(self, var, np.insert(arr, 0, val))


float_type = 'f8'
s3 = read_driver()

# MicroHH settings
zsize = 150.
ktot  = 150

# Create vertical grid
dz    = zsize/ktot
z     = np.arange(0.5*dz, zsize, dz)

# Create vertical profiles
th = np.interp(z, s3.z, s3.th)
u  = np.interp(z, s3.z, s3.u )
v  = np.interp(z, s3.z, s3.v )
ug = np.interp(z, s3.z, s3.ug)
vg = np.interp(z, s3.z, s3.vg)

# Save all the input data to NetCDF
nc_file = nc4.Dataset('gabls4s3_nbl_input.nc', mode='w', datamodel='NETCDF4', clobber=True)

nc_file.createDimension('z', ktot)
nc_z = nc_file.createVariable('z', float_type, ('z'))
nc_z[:] = z[:]

# Create a group called 'init' for the initial profiles.
nc_group_init = nc_file.createGroup('init')

nc_th = nc_group_init.createVariable('th',    float_type, ('z'))
nc_u  = nc_group_init.createVariable('u',     float_type, ('z'))
nc_v  = nc_group_init.createVariable('v',     float_type, ('z'))
nc_ug = nc_group_init.createVariable('u_geo', float_type, ('z'))
nc_vg = nc_group_init.createVariable('v_geo', float_type, ('z'))

nc_th[:] = th[:]
nc_u [:] = u [:]
nc_v [:] = v [:]
nc_ug[:] = ug[:]
nc_vg[:] = vg[:]

# Create a group called 'timedep' for the timedep.
nc_group_timedep = nc_file.createGroup('timedep')
nc_group_timedep.createDimension('time_surface', s3.t.size)

nc_time_surface = nc_group_timedep.createVariable('time_surface', float_type, ('time_surface'))
nc_th_sbot = nc_group_timedep.createVariable('th_sbot', float_type, ('time_surface'))

nc_time_surface[:] = s3.t[:]
nc_th_sbot[:] = s3.ths[:]

nc_file.close()
