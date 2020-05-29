import matplotlib.pyplot as pl
import netCDF4 as nc4
import numpy as np

pl.close('all'); pl.ion()

def grid_stretched(kmax, dz0, alpha):
    """
    Simple linearly stretched grid
    """
    dz = dz0 * (1 + alpha)**np.arange(kmax)
    zh = np.zeros(kmax+1)
    zh[1:] = np.cumsum(dz)
    z = 0.5 * (zh[1:] + zh[:-1])
    zsize = zh[-1]
    return z, zsize

def interp_z(array, z_in, z_out):
    """
    Interpolate array (with or without time dimension) to new vertical grid
    """
    if array.ndim == 1:
        return np.interp(z_out, z_in, array)
    else:
        nt = array.shape[0]
        array_out = np.zeros((nt, z_out.size))
        for i in range(nt):
            array_out[i] = np.interp(z_out, z_in, array[i,:])
        return array_out

def has_time_dim(dims):
    """
    Check if any of the dimensions in `dims` is a time dimension
    """
    has_time = False
    for dim in dims:
        if 'time' in dim:
            has_time = True
    return has_time


# Input file (`cabauw_input_hr.nc`) runs from 00:00 to 00:00 utc.
# Period (indices) used for experiment:
t0 = 6
t1 = 19

# Simple linearly stretched vertical grid:
ktot = 256
dz0 = 20
z, zsize = grid_stretched(ktot, dz0, 0.008)

# Hi-res input file from LS2D
f_hr  = nc4.Dataset('cabauw_input_hr.nc', 'r')

# New (interpolated) input file for MicroHH
f_new = nc4.Dataset('cabauw_input.nc', 'w')

# Set/copy dimensions/groups/...
f_new.createDimension('z', ktot)
f_new.createVariable('z', 'f8', ('z'))
f_new.variables['z'][:] = z

for group in f_hr.groups:
    nc_g = f_new.createGroup(group)

    for dim in f_hr[group].dimensions:
        if 'time' in dim:
            nc_g.createDimension(dim, t1-t0)
        else:
            nc_g.createDimension(dim, f_hr[group].dimensions[dim].size)

    for var in f_hr[group].variables:
        nc_var_in = f_hr[group].variables[var]
        dims = nc_var_in.dimensions
        nc_var = nc_g.createVariable(var, 'f8', dims)

        islice = np.s_[t0:t1] if has_time_dim(dims) else np.s_[:]

        if 'time' in var:
            nc_var[:] = nc_var_in[islice] - nc_var_in[t0]
        else:
            if 'z' not in dims or group == 'soil':
                nc_var[:] = nc_var_in[islice]
            else:
                nc_var[:] = interp_z(nc_var_in[islice], f_hr.variables['z'][:], z)

f_new.close()
