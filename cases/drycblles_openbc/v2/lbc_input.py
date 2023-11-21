import numpy as np
import xarray as xr

def lbc_input(
        fields,
        x, y, z,
        xh, yh, zh,
        time, ngc, dtype=np.float64):
    """
    Help function to define an Xarray Dataset with the 
    correct dimensions/coordinates/fields for the
    open BCs / nesting in MicroHH.
    """

    # Checks.
    #if 'w' in fields:
    #    raise Exception('You can not specify vertical velocity BCs!')

    nt = time.size
    itot = x.size
    jtot = y.size
    ktot = z.size

    # Dimension sizes.
    dims = {
        'time': nt,
        'x': itot + 2*ngc,
        'xh': itot + 2*ngc,
        'xgw': ngc,
        'xge': ngc,
        'xhgw': ngc+1,
        'xhge': ngc,
        'y': jtot + 2*ngc,
        'yh': jtot + 2*ngc,
        'ygn': ngc,
        'ygs': ngc,
        'yhgs': ngc+1,
        'yhgn': ngc,
        'z': ktot}

    # Coordinates.
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    # Pad x,y dimensions with ghost cells.
    xp = np.zeros(x.size+2*ngc)
    xhp = np.zeros(x.size+2*ngc)

    yp = np.zeros(y.size+2*ngc)
    yhp = np.zeros(y.size+2*ngc)

    xp[ngc:-ngc] = x
    xhp[ngc:-ngc] = xh

    yp[ngc:-ngc] = y
    yhp[ngc:-ngc] = yh

    for i in range(ngc):
        xp[i] = x[0] - (ngc-i)*dx
        xhp[i] = xh[0] - (ngc-i)*dx

        yp[i] = y[0] - (ngc-i)*dy
        yhp[i] = yh[0] - (ngc-i)*dy
        
        xp[itot+ngc+i] = x[-1] + (i+1)*dx
        xhp[itot+ngc+i] = xh[-1] + (i+1)*dx

        yp[jtot+ngc+i] = y[-1] + (i+1)*dy
        yhp[jtot+ngc+i] = yh[-1] + (i+1)*dy
    
    # Define coordinates.
    coords = {
        'time': time,
        'x': xp,
        'xh': xhp,
        'y': yp,
        'yh': yhp,
        'z': z,
        'zh': zh,
        'xgw': xp[:ngc],
        'xge': xp[itot+ngc:],
        'xhgw': xhp[:ngc+1],
        'xhge': xhp[itot+ngc:],
        'ygs': yp[:ngc],
        'ygn': yp[jtot+ngc:],
        'yhgs': yhp[:ngc+1],
        'yhgn': yhp[jtot+ngc:]}

    # Create Xarray dataset.
    ds = xr.Dataset(coords=coords)

    def get_dim_size(dim_in):
        out = []
        for dim in dim_in:
            out.append(coords[dim].size)
        return out

    def add_var(name, dims):
        dim_size = get_dim_size(dims)
        ds[name] = (dims, np.zeros(dim_size, dtype=dtype))

    for fld in fields:
        if fld not in ('u','v','w'):
            add_var(f'{fld}_west', ('time', 'z', 'y', 'xgw'))
            add_var(f'{fld}_east', ('time', 'z', 'y', 'xge'))
            add_var(f'{fld}_south', ('time', 'z', 'ygs', 'x'))
            add_var(f'{fld}_north', ('time', 'z', 'ygn', 'x'))

    if 'u' in fields:
        add_var('u_west', ('time', 'z', 'y', 'xhgw'))
        add_var('u_east', ('time', 'z', 'y', 'xhge'))
        add_var('u_south', ('time', 'z', 'ygs', 'xh'))
        add_var('u_north', ('time', 'z', 'ygn', 'xh'))

    if 'v' in fields:
        add_var('v_west', ('time', 'z', 'yh', 'xgw'))
        add_var('v_east', ('time', 'z', 'yh', 'xge'))
        add_var('v_south', ('time', 'z', 'yhgs', 'x'))
        add_var('v_north', ('time', 'z', 'yhgn', 'x'))

    if 'w' in fields:
        add_var('w_west', ('time', 'zh', 'y', 'xgw'))
        add_var('w_east', ('time', 'zh', 'y', 'xge'))
        add_var('w_south', ('time', 'zh', 'ygs', 'x'))
        add_var('w_north', ('time', 'zh', 'ygn', 'x'))

    return ds



if __name__ == '__main__':
    # Just for testing...

    xsize = 3200
    ysize = 2400
    zsize = 1600

    itot = 32
    jtot = 24
    ktot = 16

    dx = xsize / itot
    dy = ysize / jtot
    dz = zsize / ktot

    x = np.arange(dx/2, xsize, dx)
    y = np.arange(dy/2, ysize, dy)
    z = np.arange(dz/2, zsize, dz)

    xh = np.arange(0, xsize, dx)
    yh = np.arange(0, ysize, dy)
    zh = np.arange(0, zsize, dz)

    time = np.array([0, 10800], np.float64)
    ngc = 2

    lbc = lbc_input(
            ['u', 'v', 's'],
            x, y, z, xh, yh, zh, time, ngc)

    # `lbc` is now an Xarray Dataset with all the correct
    # dimensions, coordinates, fields, etc.
    lbc['u_west'][:] = 2.

    lbc.to_netcdf('test_lbc_input.nc')
