import matplotlib.pyplot as pl
import numpy as np
import netCDF4 as nc

from microhh_tools import Read_namelist

if __name__ == '__main__':

    float_type = 'f8'
    int_type = 'i4'

    # Read the namelist
    nl = Read_namelist()

    # Grid information
    zsize = nl['grid']['zsize']
    ktot  = nl['grid']['ktot']
    dz = zsize / ktot

    # Define input atmospheric fields
    z   = np.linspace(0.5*dz, zsize-0.5*dz, ktot)
    thl = np.zeros(np.size(z))
    qt  = np.zeros(np.size(z))
    u   = np.zeros(np.size(z))
    ug  = np.zeros(np.size(z))

    zi = 500.
    thl[z<zi ] = 290
    thl[z>=zi] = 290 + (z[z>=zi]-zi)*0.006

    qt[z<zi ] = 8e-3
    qt[z>=zi] = 8e-3 - (z[z>=zi]-zi)*0.000001

    u[:]  = 5.
    ug[:] = 5.

    # Define soil fields
    ktot_soil = nl['soil']['ktot']
    z_soil = np.array([-1, -0.5, -0.25, -0.1])
    theta_soil = 0.2+0.1*np.arange(ktot_soil)
    temp_soil = 289+np.arange(ktot_soil)
    soil_index = np.zeros(ktot_soil, dtype=int)

    # Save all the input data to NetCDF
    nc_file = nc.Dataset('lsm_dev_input.nc', mode='w', datamodel='NETCDF4')

    nc_file.createDimension('z', ktot)
    nc_file.createDimension('z_soil', ktot_soil)

    nc_z = nc_file.createVariable('z', float_type, ('z'))
    nc_z_soil = nc_file.createVariable("z_soil", float_type, ("z_soil"))

    nc_z[:] = z[:]
    nc_z_soil[:]     = z_soil[:]

    # Create a group called "init" for the initial profiles.
    nc_group_init = nc_file.createGroup('init')

    nc_thl = nc_group_init.createVariable('thl', float_type, ('z'))
    nc_qt  = nc_group_init.createVariable('qt' , float_type, ('z'))
    nc_u   = nc_group_init.createVariable('u'  , float_type, ('z'))
    nc_ug  = nc_group_init.createVariable('ug' , float_type, ('z'))

    nc_thl[:] = thl[:]
    nc_qt [:] = qt [:]
    nc_u  [:] = u  [:]
    nc_ug [:] = ug [:]

    # Create soil group
    nc_group_soil = nc_file.createGroup("soil")

    nc_soil_index = nc_group_soil.createVariable("soil_index", int_type, ("z_soil"))
    nc_theta_soil = nc_group_soil.createVariable("theta_soil", float_type, ("z_soil"))
    nc_temp_soil  = nc_group_soil.createVariable("t_soil", float_type, ("z_soil"))

    nc_soil_index[:] = soil_index[:]
    nc_theta_soil[:] = theta_soil[:]
    nc_temp_soil[:]  = temp_soil[:]

    nc_file.close()
