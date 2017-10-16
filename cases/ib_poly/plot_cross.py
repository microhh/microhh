import numpy as np
import matplotlib.pylab as pl
import matplotlib
import netCDF4 as nc4

pl.close('all')

class Readnc:
    def __init__(self, file_name):
        f = nc4.Dataset(file_name, 'r')
        for v in f.variables:
            setattr(self, v, f.variables[v][:])

def plot_polys(poly):
    ax=pl.gca()
    for p in poly:
        patch = matplotlib.patches.Polygon(np.vstack((p['x'], p['y'])).transpose(),\
                                           facecolor='k', edgecolor='k', alpha=0.8, zorder=10)
        ax.add_patch(patch)

def mask(field, m):
    field = np.ma.masked_array(field)
    field.mask = m
    return field


if __name__ == "__main__":
    import microhh_tools as mht
    import ib_tools as ibt
    import ib_tools_cython as ibtc

    # Read the namelist, and grid info:
    nl = mht.Read_namelist()
    gr = mht.Read_grid(nl['grid']['itot'], nl['grid']['jtot'], nl['grid']['ktot'], nl['grid']['zsize'])

    # Read polyons
    poly = ibt.read_polygon_unstructured('poly_coordinates.txt')

    # Flag in/outside IB
    inside_s = ibtc.flag_inside_ib_2d(gr.x,  gr.y,  poly)
    m = np.zeros((nl['grid']['jtot'], nl['grid']['itot']), dtype=bool)
    m[inside_s != -1] = True

    # Read cross-sections
    uxy = Readnc('u.xy.nc')
    vxy = Readnc('v.xy.nc')
    wxy = Readnc('w.xy.nc')
    uxz = Readnc('u.xz.nc')

    pl.figure()
    pl.subplot(111, aspect='equal')
    pl.pcolormesh(uxy.xh, uxy.y, mask(uxy.u[-1,0,:,:], m), vmin=-0.4, vmax=0.8, cmap=pl.cm.magma)
    pl.xlabel('x')
    pl.ylabel('y')
    pl.colorbar()

    pl.figure()
    pl.subplot(111, aspect='equal')
    pl.streamplot(uxy.xh, uxy.y, mask(uxy.u[-1,0,:,:], m), vxy.v[-1,0,:,:], density=4, linewidth=0.5)
    plot_polys(poly)
