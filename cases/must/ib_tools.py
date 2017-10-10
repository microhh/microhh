import numpy as np
import matplotlib.pylab as pl
import time
import sys

pl.close('all')

def is_left(xp, yp, x0, y0, x1, y1):
    """
    Check whether point (xp,yp) is left of line segment ((x0,y0) to (x1,y1))
    returns:  >0 if left of line, 0 if on line, <0 if right of line
    """

    return (x1-x0) * (yp-y0) - (xp-x0) * (y1-y0)


def is_inside(xp, yp, x_set, y_set):
    """
    Given location (xp,yp) and set of line segments (x_set, y_set), determine
    whether (xp,yp) is inside polygon.
    """

    wn = 0
    for i in range(x_set.size-1):
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

    # First simple check
    if (xp < x_set.min() or xp > x_set.max() or yp < y_set.min() or yp > y_set.max()):
        return False

    wn = np.ones(x_set.size-1)*1000

    left   = is_left(xp, yp, x_set[:-1], y_set[:-1], x_set[1:], y_set[1:])
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


def read_polygon_structured(file_name, xoffs=0, yoffs=0):
    """
    Read the polygon input for MicroHH, for structured input (each polygon has same
    number of vertices).
    Each object is defined on a new line, with e.g.: z, x1, y1, x2, y2, etc.
    """
    print('Reading {}'.format(file_name))

    # Open file with Numpy, and read the dimensions.
    f = np.loadtxt(file_name, delimiter=',')
    n_vert = int((f[0,:].size-1)/2)
    n_obj  = f[:,0].size
    print('Found {} objects with {} vertices each'.format(n_obj, n_vert))

    # Check if polygons are closed. For the structured read_polygon(), either all polygons have to be closed or not...
    x = f[0,:][1::2]
    y = f[0,:][2::2]

    print(x,y)

    #if (x[0] != x[-1] and y[0] != y[-1]):
    #    print('Polygons not closed, will close them...')
    #    close_poly = True
    #else:
    #    close_poly = False

    ## Load all the polygons and store in list with x,y combination dicts
    #polygons = []
    #for n in range(n_obj):
    #    x = f[n,:][0::2] + xoffs
    #    y = f[n,:][1::2] + yoffs

    #    if close_poly:
    #        x = np.append(x, x[0])
    #        y = np.append(y, y[0])

    #    polygons.append({'x':x, 'y':y})

    #return polygons


def find_ghost_cells_2d(inside_ib):
    """
    Find the ghost cells, i.e. a grid point which is inside the IB,
    which has a neighbour grid point outside the IB (where inside_ib = -1)
    """

    ghost = np.ones_like(inside_ib, dtype=int) * -1

    ij_in = np.where(inside_ib != -1)
    for n in range(ij_in[0].size):
        # Index in 2d grid
        j = ij_in[0][n]
        i = ij_in[1][n]

        # Don't allow IB near domain edge, otherwhise we can't 
        # search for e.g. interpolation points
        if (i < 2 or i > gr.x.size-3 or \
            j < 2 or j > gr.y.size-3):
            sys.exit('NEIN!, container too close to edge of domain')
        else:
            if (inside_ib[j,i-1] == -1 or inside_ib[j,i+1] == -1 or \
                inside_ib[j-1,i] == -1 or inside_ib[j+1,i] == -1):
                    # Give ghost cell same value (object id):
                    ghost[j,i] = inside_ib[j,i]

    return ghost

def find_interpolation_points(xi, yi, zi, x, y, z, ib_mask, n_idw):
    """
    Find the "n_idw" nearest interpolation points outside the IB
    """

    # Limit vertical search stencil near surface
    dk0 = -4; if (k+dk0 < 0): dk0=-k

    # Lists holding indexes and distance
    i_n = []
    j_n = []
    k_n = []
    d   = []

    for di in range(-2, 3):
        for dj in range(2, 3):
            for dk in range(dk0, 5):
                if (ib_mast == -1):
                    dist = abs_distance(tmp_ghost['xi'], tmp_ghost['yi'], tmp_ghost['zi'],x[i+di], y[j+dj], z[k+dk])
                    i_n.append(i+di)
                    j_n.append(j+dj)
                    k_n.append(k+dk)
                    d  .append(dist)

    # Sort on distance, and clip to the requested amount of interpolation points
    inds = np.array(d).argsort()
    d    = np.array(d)  [inds][:n_idw]
    i_n  = np.array(i_n)[inds][:n_idw]
    j_n  = np.array(j_n)[inds][:n_idw]
    k_n  = np.array(k_n)[inds][:n_idw]

    return i_n, j_n, k_n

def calc_ghost_cell_properties(ghost_cells, inside_ib, polygons, x, y, z, kmax, n_idw):
    # List to hold all the ghost cell info
    ghost_list = []

    # Loop over all points inside IB
    ij_in = np.where(inside_ib != -1)
    for n in range(ij_in[0].size):
        j = ij_in[0][n]
        i = ij_in[1][n]
        o = ghost_cells[j,i]    # Object id

        # Lateral wall: pick nearest location on vertical wall
        if (ghost_cells[j,i] != -1):
            xb, yb, zb =  nearest_on_polygon(x[i], y[j], 99999, polygons[o]['x'], polygons[o]['y'], -99999)

            # Loop over object height
            for k in range(kmax):

                # Interpolation location:
                xi = 2*xb - x[i]
                yi = 2*yb - y[j]
                zi = z[k]

                # Find "n_idw" nearest grid points outside IB
                i_n, j_n, k_n = find_interpolation_points(xi, yi, zi, x, y, z, inside_ib, n_idw)

                # Add to dict with ghost cells
                ghost_list.append({'i':i, 'j':j, 'k':k,
                                   'xb':xb, 'yb':yb, 'zb':zb,
                                   'xi':xi, 'yi':yi, 'zi':zi,
                                   'in':i_n, 'jn':j_n, 'kn':k_n})

    return ghost_list


if __name__ == "__main__":
    import microhh_tools as mht

    z_container = 2.54  # container height (m)

    # Read the namelist, and grid info
    nl = mht.Read_namelist()
    gr = mht.Read_grid(nl['grid']['itot'], nl['grid']['jtot'], nl['grid']['ktot'], nl['grid']['zsize'])

    # Offset inside model domain compared to input location containers
    x_offset = 25
    y_offset = 25

    # Read the polygon input
    poly = read_polygon_structured('container_corners.txt', x_offset, y_offset)

    """

    # Subset for testing
    poly = poly[:10]

    # Find grid points inside IB (2D)
    # Scalar location is the same {x,y} as w location on a 2D plane
    # Skip this part if executed from IPython with run -i ... (interactive mode),
    # as this is pretty expensive / slow :-(
    if "inside_u_2d" not in locals():

        # Flags for outside (-1) or inside (>= 0 ->container ID) IB
        inside_u_2d = np.ones((gr.y.size, gr.x.size), dtype=int) * -1
        inside_v_2d = np.ones((gr.y.size, gr.x.size), dtype=int) * -1
        inside_w_2d = np.ones((gr.y.size, gr.x.size), dtype=int) * -1

        for n in range(len(poly)):
            print('Working on object {0:4d} of {1:4d}'.format(n+1, len(poly)))
            for i in range(gr.x.size):
                for j in range(gr.y.size):
                    # u-location
                    if (is_inside_vec(gr.xh[i], gr.y[j],  poly[n]['x'], poly[n]['y'])):
                        inside_u_2d[j,i] = n
                    # v-location
                    if (is_inside_vec(gr.x[i],  gr.yh[j], poly[n]['x'], poly[n]['y'])):
                        inside_v_2d[j,i] = n
                    # w and s-location
                    if (is_inside_vec(gr.x[i],  gr.y[j],  poly[n]['x'], poly[n]['y'])):
                        inside_w_2d[j,i] = n

    # Find the lateral ghost cells
    ghost_u_2d = find_ghost_cells_2d(inside_u_2d)
    ghost_v_2d = find_ghost_cells_2d(inside_v_2d)
    ghost_w_2d = find_ghost_cells_2d(inside_w_2d)

    # Find x,y location nearest walls
    wall_u_x, wall_u_y = find_nearest_wall_2d(ghost_u_2d, poly, gr.xh, gr.y)
    wall_v_x, wall_v_y = find_nearest_wall_2d(ghost_v_2d, poly, gr.x, gr.yh)
    wall_w_x, wall_w_y = find_nearest_wall_2d(ghost_w_2d, poly, gr.x, gr.y)

    # From 2d to 3d....
    # Find first grid levels above containers
    kmax_block  = np.abs(gr.z -z_container).argmin()   # Full grid level
    kmaxh_block = np.abs(gr.zh-z_container).argmin()   # Half grid level

    if gr.z[kmax_block] < z_container:
        kmax_block += 1
    if gr.zh[kmaxh_block] < z_container:
        kmaxh_block += 1

    # Calculate all ghost cell properties (nearest location wall, interpolation points, ....)
    ghost_u = calc_ghost_cell_properties(ghost_u_2d, inside_u_2d, poly, gr.xh, gr.y, gr.z, kmax_block)

    """


    





    pl.figure()
    pl.subplot(221, aspect='equal')
    pl.title('Inside IB / container ID', loc='left')
    pl.pcolormesh(gr.x, gr.y, inside_w_2d)
    pl.colorbar()

    pl.subplot(222, aspect='equal')
    pl.title('Ghost cells', loc='left')
    pl.pcolormesh(gr.x, gr.y, ghost_w_2d)
    pl.colorbar()

    pl.subplot(223, aspect='equal')
    pl.title('Wall_x', loc='left')
    pl.pcolormesh(gr.x, gr.y, wall_u_x)
    pl.colorbar()

    pl.subplot(224, aspect='equal')
    pl.title('Wall_y', loc='left')
    pl.pcolormesh(gr.x, gr.y, wall_u_y)
    pl.colorbar()




    #""" Test test """
    #xysize = 1600
    #ijtot  = 16
    #dxy    = xysize/ijtot

    ## Simplified grid at full and half (h) locations
    #x  = np.arange(0.5*dxy, ijtot*dxy, dxy)
    #y  = np.arange(0.5*dxy, ijtot*dxy, dxy)
    #xh = np.arange(0, ijtot*dxy+0.1, dxy)
    #yh = np.arange(0, ijtot*dxy+0.1, dxy)

    ## Create some buildings...
    #set_x1 = np.array([300, 1000, 800, 1300, 400,  300])
    #set_y1 = np.array([300, 500,  700, 1400, 1000, 300])
    #set_x2 = np.array([1200, 1400, 1400, 1200, 1200])
    #set_y2 = np.array([300,  300,  500,  500,  300])

    ## Plot, and set ticks equal to grid half locations
    #pl.figure()
    #ax = pl.subplot(111)
    #ax.set_aspect(1)
    #ax.set_xticks(xh)
    #ax.set_yticks(yh)
    #ax.set_xlim(0,xysize)
    #ax.set_ylim(0,xysize)
    #ax.grid(linestyle='dotted')

    ## Loop over the different buildings:
    #for set_x,set_y in zip([set_x1,set_x2],[set_y1,set_y2]):
    #
    #    # Plot the building
    #    pl.plot   (set_x, set_y, 'r-', linewidth=1.5, dashes=[4,2])
    #    pl.scatter(set_x, set_y, marker='+', s=50, facecolor='r', edgecolor='r')
    #
    #    # For all grid points; check wether it is inside a block
    #    gp    = np.zeros((ijtot,ijtot), dtype=bool)
    #    for i in range(1,ijtot-1):
    #        for j in range(1,ijtot-1):
    #            if(is_inside(x[i], y[j], set_x, set_y)):
    #                gp[i,j] = True
    #
    #    # Loop over the points which are inside a block
    #    ij_in = np.where(gp)
    #    for n in range(ij_in[0].size):
    #        i = ij_in[0][n]
    #        j = ij_in[1][n]

    #        # Check if it is a ghost cell:
    #        if (gp[i-1,j] == False or gp[i+1,j] == False or gp[i,j-1] == False or gp[i,j+1] == False):
    #            # Find nearest location on the wall:
    #            xb,yb,zb = nearest_on_polygon(x[i], y[j], 0., set_x, set_y, z=100000.)

    #            pl.scatter(x[i], y[j], s=30, facecolor='r')
    #            pl.plot([x[i],xb],[y[j],yb], 'k-')
    #        else:
    #            pl.scatter(x[i], y[j], s=30, facecolor='g')
