import numpy as np
import netCDF4 as nc
from scipy.special import erf

# Available in microhh/python
import microhh_tools as mht

float_type = 'f8'


def pad_profile(z, dz, hc=10.0, pai=0.479, sigma=0.5):
    """
    PAD profile with smooth canopy top using error function.
    Profile is scaled such that the integral of PAD matches the input PAI.
    """
    a_rel = 0.5 * (1 - erf((z - hc) / sigma))
    a_base = pai / (np.sum(a_rel) * dz)
    return a_rel * a_base


ini = mht.Read_namelist('canopy.ini')

# Initial vertical profiles.
zsize = ini['grid']['zsize']
ktot = ini['grid']['ktot']

dz = zsize / ktot
z = np.arange(dz/2, zsize, dz)

u  = np.zeros(ktot) + 6
v  = np.zeros(ktot)


# Plant area density at full levels. Code interpolates to the `u`, `v`, and `w` locations internally.
hc = 10.        # Canopy height (m)
pai = 0.479     # Integrated plant area density PAI = int(a * dz)
sigma = 1       # Width transition are canopy top (m)

pad = pad_profile(z, dz, hc, pai, sigma)


# Create case_input.nc
nc_file = nc.Dataset('canopy_input.nc', mode='w', datamodel='NETCDF4', clobber=True)

nc_file.createDimension('z', ktot)

nc_z = nc_file.createVariable('z', float_type, ('z'))

nc_group_init = nc_file.createGroup('init');
nc_u = nc_group_init.createVariable('u', float_type, ('z'))
nc_v = nc_group_init.createVariable('v', float_type, ('z'))
nc_pad = nc_group_init.createVariable('pad', float_type, ('z'))

nc_z[:] = z[:]

nc_u[:] = u[:]
nc_v[:] = v[:]
nc_pad[:] = pad[:]

nc_file.close()
