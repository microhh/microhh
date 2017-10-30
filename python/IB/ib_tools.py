import numpy as np
import matplotlib.pylab as pl
import time
import sys

pl.close('all')

class Timer:
    def __init__(self):
        self.t0 = time.time()

    def reset(self):
        self.t1 = time.time()
        dt = self.t1-self.t0
        print('elapsed = {0:.1f} sec'.format(dt))
        self.t0 = time.time()


def is_left(xp, yp, x0, y0, x1, y1):
    """
    Check whether point (xp,yp) is left of line segment ((x0,y0) to (x1,y1))
    returns:  >0 if left of line, 0 if on line, <0 if right of line
    """

    return (x1-x0) * (yp-y0) - (xp-x0) * (y1-y0)


def is_inside(xp, yp, x_set, y_set):
    """
    Given location (xp,yp) and set of line segments (x_set, y_set), determine
    whether (xp,yp) is inside (or on) polygon.
    """

    # First simple check on bounds
    if (xp < x_set.min() or xp > x_set.max() or yp < y_set.min() or yp > y_set.max()):
        return False

    wn = 0
    for i in range(x_set.size-1):
        # Second check: see if point exactly on line segment:
        if (is_left(xp, yp, x_set[i], y_set[i], x_set[i+1], y_set[i+1]) == 0):
            return True

        # Calculate winding number
        if (y_set[i] <= yp):
            if (y_set[i+1] > yp):
                if (is_left(xp, yp, x_set[i], y_set[i], x_set[i+1], y_set[i+1]) > 0):
                    wn += 1
        else:
            if (y_set[i+1] <= yp):
                if (is_left(xp, yp, x_set[i], y_set[i], x_set[i+1], y_set[i+1]) < 0):
                    wn -= 1

    if (wn == 0):
        return False
    else:
        return True


def is_inside_vec(xp, yp, x_set, y_set):
    """
    Vectorized version of is_inside()
    """

    left   = is_left(xp, yp, x_set[:-1], y_set[:-1], x_set[1:], y_set[1:])

    # Return true if exactly on line
    if (np.any(left == 0)):
        return True

    ydiff1 = y_set[:-1]-yp
    ydiff2 = y_set[1:] -yp

    d1 = ((ydiff1 <= 0) & (ydiff2 >  0) & (left > 0)) *  1
    d2 = ((ydiff1 >  0) & (ydiff2 <= 0) & (left < 0)) * -1

    return (np.sum(d1) + np.sum(d2)) != 0


def nearest_on_line(xp, yp, x0, y0, x1, y1):
    """
    Given location (xp,yp) and line segment from {x0,y0} to {x1,y1},
    find closest location of (xp,yp) on line segment.
    """

    dist_l_2 = (x1-x0)**2 + (y1-y0)**2                      # Length line segment squared
    dot_l0_p_l0_l1 = (xp-x0)*(x1-x0) + (yp-y0)*(y1-y0)      # Dot product vector p0-p,p0-p1
    t  = dot_l0_p_l0_l1 / dist_l_2                          # Normalized distance
    t  = min(1,max(0,t))                                    # Limit distance
    xb = x0+(x1-x0)*t                                       # Closest location on line segment
    yb = y0+(y1-y0)*t                                       # Closest location on line segment
    d  = ((xp-xb)**2 + (yp-yb)**2)**0.5                     # Distance to {xp, yp}

    return xb,yb,d


def nearest_on_polygon(xp, yp, zp, x_set, y_set, z):
    """
    Given location (xp,yp) and set of line segments, find closest location
    on any of the line segments.
    """

    # 1. Check if set is closed; if not, close the set
    if x_set[-1] != x_set[0]:
        x_set = np.append(x_set, x_set[0])
        y_set = np.append(y_set, y_set[0])

    # 2. Find the nearest location on the set of line segments (polygon) or roof
    d_min = np.abs(zp-z)
    xb    = xp
    yb    = yp
    zb    = z
    for i in range(x_set.size-1):
        x_tmp, y_tmp, d_tmp = nearest_on_line(xp, yp, x_set[i], y_set[i], x_set[i+1], y_set[i+1])
        if (d_tmp < d_min):
            d_min = d_tmp
            xb    = x_tmp
            yb    = y_tmp
            zb    = zp
    return xb,yb,zb


def abs_dist(x1, y1, z1, x2, y2, z2):
    return ((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)**0.5


def read_polygon_unstructured(file_name, xoffs=0, yoffs=0):
    """
    Read the polygon input for MicroHH, for unstructured input (each polygon can
    have a different number of vertices).
    Each object is defined on a new line, with e.g.: z, x1, y1, x2, y2, etc.
    """
    print('Reading {}'.format(file_name))

    polygons = []

    # Parse data line-by-line as the input per object can differ in size
    f = open(file_name, 'r')
    for l in f.readlines():

        # Skip comments
        if '#' not in l:

            # Split line on comma, and read height and x,y coordinates
            data = l.split(',')
            z = float(data[0])
            x = np.array(data[1::2]).astype(float) + xoffs
            y = np.array(data[2::2]).astype(float) + yoffs

            # Check if polygon is closed, if not close it
            if (x[0] != x[-1]) or (y[0] != y[-1]):
                x = np.append(x,x[0])
                y = np.append(y,y[0])

            # Add to list of polygons
            polygons.append({'z':z, 'x':x, 'y':y})

    f.close()

    return polygons


def flag_inside_ib_2d(x, y, poly):
    """
    Find grid points which are inside the IB
    all grid points are flagged -1 (outside IB), or otherwise
    given the ID of the polygon
    """

    # Initialise all values at -1 (outside IB)
    inside = np.ones((y.size, x.size), dtype=int) * -1

    # Painfully slow, but required for strangly shaped polygons
    for id in range(len(poly)):
        print('Working on object {0:4d} of {1:4d}'.format(id+1, len(poly)))
        for i in range(x.size):
            for j in range(y.size):
                if (is_inside(x[i], y[j],  poly[id]['x'], poly[id]['y'])):
                    inside[j,i] = id

    return inside


def flag_ghost_cells_2d(inside_ib):
    """
    Find the ghost cells, i.e. a grid point which is inside the IB,
    which has a neighbour grid point outside the IB (where inside_ib = -1)
    """

    ghost = np.ones_like(inside_ib, dtype=int) * -1
    itot  = ghost.shape[1]
    jtot  = ghost.shape[0]

    ij_in = np.where(inside_ib != -1)
    for n in range(ij_in[0].size):
        # Index in 2d grid
        j = ij_in[0][n]
        i = ij_in[1][n]
        o = inside_ib[j,i]

        # Check check
        assert o != -1

        # Don't allow IB near domain edge, otherwhise we can't
        # search for e.g. interpolation points
        if (i < 2 or i > itot-3 or \
            j < 2 or j > jtot-3):
            sys.exit('NEIN!, object too close to edge of domain')
        else:
            if (inside_ib[j,i-1] == -1 or inside_ib[j,i+1] == -1 or \
                inside_ib[j-1,i] == -1 or inside_ib[j+1,i] == -1):
                    # Give ghost cell value of object id:
                    ghost[j,i] = o

    return ghost


def find_interpolation_points(xi, yi, zi, x, y, z, i, j, k, ib_mask, n_idw, kmax):
    """
    Find the "n_idw" nearest interpolation points outside the IB
    """

    # Limit vertical search stencil near surface
    dk0 = -4
    if (k+dk0 < 0):
        dk0=-k

    # Lists holding indexes and distance
    i_n = []
    j_n = []
    k_n = []
    d   = []

    dx = x[1]-x[0]
    dy = y[1]-y[0]
    dz = z[k+1]-z[k]
    mind = np.min([dx,dy,dz])

    # Loop over a +/- 2 grid point stencil to
    # find potential interpolation points
    for di in range(-2, 3):
        for dj in range(-2, 3):
            for dk in range(dk0, 5):
                if (ib_mask[j+dj,i+di] == -1 or k+dk >= kmax):
                        dist = abs_dist(xi, yi, zi, x[i+di], y[j+dj], z[k+dk])
                        i_n.append(i+di)
                        j_n.append(j+dj)
                        k_n.append(k+dk)
                        d  .append(dist)

    # Sort on distance, and clip to the requested amount of interpolation points
    inds = np.array(d).argsort()
    d_tmp = np.array(d)[inds]

    # Exclude grid points which are too close to the boundary
    too_close = np.where(d_tmp < 0.1*mind)[0].size

    d    = np.array(d)  [inds][too_close:n_idw+too_close]
    i_n  = np.array(i_n)[inds][too_close:n_idw+too_close]
    j_n  = np.array(j_n)[inds][too_close:n_idw+too_close]
    k_n  = np.array(k_n)[inds][too_close:n_idw+too_close]

    return i_n, j_n, k_n


def calc_ghost_cell_properties(ghost_cells, inside_ib, polygons, x, y, z, n_idw, label='', kstart=0):
    print('Calculation ghost cell properties for: {}'.format(label))

    # List to hold all the ghost cell info
    ghost_list = []

    # Loop over all points inside IB
    ij_in = np.where(inside_ib != -1)
    for n in range(ij_in[0].size):
        j = ij_in[0][n]
        i = ij_in[1][n]
        o = inside_ib[j,i]    # Object id

        # Check check
        assert o != -1

        # Lateral wall: pick nearest location on vertical wall
        if (ghost_cells[j,i] != -1):
            xb, yb, zb =  nearest_on_polygon(x[i], y[j], 99999, polygons[o]['x'], polygons[o]['y'], -99999)

            # First vertical index above object:
            kmax  = np.abs(z-polygons[o]['z']).argmin()
            if z[kmax] < polygons[o]['z']:
                kmax += 1

            # Loop over object height
            for k in range(kstart, kmax):

                # Interpolation location:
                xi = 2*xb - x[i]
                yi = 2*yb - y[j]
                zi = z[k]

                # Nearest grid point to interpolation location
                ii = np.abs(x-xi).argmin()
                jj = np.abs(y-yi).argmin()

                # Lateral wall, so boundary location height == grid point height
                zb = z[k]

                # Find "n_idw" nearest grid points outside IB
                i_n, j_n, k_n = find_interpolation_points(xi, yi, zi, x, y, z, ii, jj, k, inside_ib, n_idw, kmax)

                # Check whether there are enough interpolation points
                assert i_n.size == n_idw

                # Add to list of ghost cells as a dict
                ghost_list.append({'i':i, 'j':j, 'k':k,
                                   'xb':xb, 'yb':yb, 'zb':zb,
                                   'xi':xi, 'yi':yi, 'zi':zi,
                                   'in':i_n, 'jn':j_n, 'kn':k_n})
        # Top boundary
        else:
            k = kmax-1

            # Top wall, so x,y location is directly above ghost cell
            xb = x[i]
            yb = y[j]
            zb = polygons[o]['z']

            # Interpolation location:
            xi = x[i]
            yi = y[j]
            zi = 2*zb-z[k]

            # Nearest grid point to interpolation location
            kk = np.abs(z-zi).argmin()

            # Find "n_idw" nearest grid points outside IB
            i_n, j_n, k_n = find_interpolation_points(xi, yi, zi, x, y, z, i, j, kk, inside_ib, n_idw, kmax)

            # Check whether there are enough interpolation points
            assert i_n.size == n_idw

            # Add to list of ghost cells as a dict
            ghost_list.append({'i':i, 'j':j, 'k':k,
                               'xb':xb, 'yb':yb, 'zb':zb,
                               'xi':xi, 'yi':yi, 'zi':zi,
                               'in':i_n, 'jn':j_n, 'kn':k_n})

    return ghost_list


def write_output(ghost_cells, file_name):
    print('Writing {}.....'.format(file_name))
    n_neighbours = ghost_cells[0]['in'].size

    f = open(file_name, 'w')
    f.write('{0:^4s} {1:^4s} {2:^4s} {3:^20s} {4:^20s} {5:^20s} {6:^4s}'.format('i','j','k','xb','yb','zb','nn'))
    for n in range(n_neighbours):
        f.write(' {0:>2s}{1:<2d} {2:>2s}{3:<2d} {4:>2s}{5:<2d}'.format('i',n,'j',n,'k',n))
    f.write('\n')

    for g in ghost_cells:
        f.write('{0:04d} {1:04d} {2:04d} {3:1.14E} {4:1.14E} {5:1.14E} {6:04d}'\
            .format(g['i'], g['j'], g['k'], g['xb'], g['yb'], g['zb'], n_neighbours))

        for n in range(n_neighbours):
            f.write(' {0:04d} {1:04d} {2:04d}'.format(g['in'][n], g['jn'][n], g['kn'][n]))
        f.write('\n')

    f.close()


if __name__ == '__main__':
    x = np.array([0,1,1,0,0])
    y = np.array([0,0,1,1,0])

    print(True,  is_inside(0.5, 0.5, x, y))   # in
    print(True,  is_inside(0.0, 0.5, x, y))   # boundary x0
    print(True,  is_inside(1.0, 0.5, x, y))   # boundary x1
    print(True,  is_inside(0.5, 0.0, x, y))   # boundary y0
    print(True,  is_inside(0.5, 1.0, x, y))   # boundary y1
    print(False, is_inside(1.5, 0.5, x, y))   # out
