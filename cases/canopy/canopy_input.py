import numpy as np
import netCDF4 as nc
from scipy.special import erf

# Available in microhh/python
import microhh_tools as mht

float_type = np.float32


def pad_profile(z, dz, hc=10.0, pai=0.479, sigma=0.5):
    """
    PAD profile with smooth canopy top using error function.
    Profile is scaled such that the integral of PAD matches the input PAI.
    """
    a_rel = 0.5 * (1 - erf((z - hc) / sigma))
    a_base = pai / (np.sum(a_rel) * dz)
    return a_rel * a_base


ini = mht.Read_namelist('canopy.ini')

# Switch between a horizontally homogeneous canopy (1D `pad` profile in the NetCDF input),
# or a 3D `pad` field with a block of canopy in the center of the domain (binary `pad.0000000`).
sw_3d_pad = ini['canopy']['sw_3d_pad']

# Initial vertical profiles.
itot = ini['grid']['itot']
jtot = ini['grid']['jtot']
ktot = ini['grid']['ktot']

xsize = ini['grid']['xsize']
ysize = ini['grid']['ysize']
zsize = ini['grid']['zsize']

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


if sw_3d_pad:
    # Same vertical profile, but only in a block of 1/3 of the domain width in the center of the domain.
    dx = xsize / itot
    dy = ysize / jtot

    x = np.arange(dx/2, xsize, dx)
    y = np.arange(dy/2, ysize, dy)

    mask_x = (x > xsize/3) & (x < 2*xsize/3)
    mask_y = (y > ysize/3) & (y < 2*ysize/3)
    mask_xy = mask_y[:,None] & mask_x[None,:]

    # Vertical extent of the 3D field: up to the highest non-zero `pad`.
    ktot_canopy = np.nonzero(pad)[0].max() + 1
    mht.replace_namelist_value('ktot_canopy', ktot_canopy, group='canopy', namelist_file='canopy.ini')

    pad_3d = np.zeros((ktot_canopy, jtot, itot), dtype=float_type)
    pad_3d[:,:,:] = pad[:ktot_canopy,None,None] * mask_xy[None,:,:]

    # NOTE: `float_type` has to match the precision of the MicroHH build.
    pad_3d.tofile('pad.0000000')
