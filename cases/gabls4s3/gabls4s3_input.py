import numpy as np
import pylab as pl
import netCDF4 as nc4

pl.close('all')

# Create stretched grid
class Grid:
    def __init__(self, kmax, nloc1, nbuf1,dz1, dz2):
        self.kmax  = kmax
        dn         = 1./kmax
        n          = np.linspace(dn, 1.-dn, kmax)
        nloc1     *= dn
        nbuf1     *= dn
        dzdn1      = dz1/dn
        dzdn2      = dz2/dn

        dzdn       = dzdn1 + 0.5*(dzdn2-dzdn1)*(1. + np.tanh((n-nloc1)/nbuf1))
        self.dz    = dzdn*dn

        self.kmax  = kmax
        self.z     = np.zeros(self.dz.size)
        stretch    = np.zeros(self.dz.size)

        self.z[0]  = 0.5*self.dz[0]
        stretch[0] = 1.

        for k in range(1, self.kmax):
              self.z [k] = self.z[k-1] + 0.5*(self.dz[k-1]+self.dz[k])
              stretch[k] = self.dz[k]/self.dz[k-1]

        zsize = self.z[kmax-1] + 0.5*self.dz[kmax-1]
        #print('kmax=%i, zsize=%f'%(kmax,zsize))

# Read the stage 3 driver file
class read_driver:
    def __init__(self):
        nc = nc4.Dataset('SCM_LES_STAGE3.nc', 'r')
        self.t = nc.variables['time'][:]
        self.z = nc.variables['height'][:][::-1]
        self.p = nc.variables['pf'][:][::-1]

        # initial profiles
        self.th = nc.variables['theta'][:][::-1]
        self.T  = nc.variables['t'][:][::-1]
        self.q  = nc.variables['qv'][:][::-1] # = zero
        self.u  = nc.variables['u'][:][::-1]
        self.v  = nc.variables['v'][:][::-1]

        # time varying forcings
        self.ug   = nc.variables['Ug'][0,:][::-1] # u geo wind = constant in time
        self.vg   = nc.variables['Vg'][0,:][::-1] # v geo wind = consant in time
        self.advT = nc.variables['hadvT'][0,:][::-1] # = zero
        self.advq = nc.variables['hadvQ'][0,:][::-1] # = zero
        self.Ts   = nc.variables['Tg'][:]

        self.ps  = nc.variables['psurf'].getValue()
        self.z0m = nc.variables['z0m'].getValue()

        # Calculate theta_s
        self.ths = self.Ts / (self.ps / 1.e5)**(287.04/1005.)

outname = 'gabls4s3'
float_type = 'f8'

s3 = read_driver()

# Large domain (~1 km high):
g20l = Grid(288, 250, 20, 2, 12) # dz = 2 m
g10l = Grid(512, 470, 30, 1, 12) # dz = 1 m

# Restart domain (~200 m high):
g20s  = Grid(128, 65,  10, 2.0,  8) # dz = 2 m
g10s  = Grid(192, 135, 20, 1.0,  8) # dz = 1 m
g05s  = Grid(320, 245, 30, 0.5,  8) # dz = 0.5 m
g02s  = Grid(512, 440, 40, 0.25, 8) # dz = 0.25 m
g02s2 = Grid(480, 410, 30, 0.25, 8) # dz = 0.25 m # mistral

# Switch between vertical grids:
grid = g20l

# Interpolate GABLS4 data to LES grid
th = np.interp(grid.z, s3.z, s3.th)
u  = np.interp(grid.z, s3.z, s3.u)
v  = np.interp(grid.z, s3.z, s3.v)
ug = np.interp(grid.z, s3.z, s3.ug)
vg = np.interp(grid.z, s3.z, s3.vg)

# Save all the input data to NetCDF
nc_file = nc4.Dataset('gabls4s3_input.nc', mode='w', datamodel='NETCDF4', clobber=True)

nc_file.createDimension('z', grid.kmax)
nc_z = nc_file.createVariable('z', float_type, ('z'))
nc_z[:] = grid.z[:]

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
