import numpy as np
import matplotlib.pylab as pl
import time

pl.close('all')

def is_left(xp, yp, x0, y0, x1, y1):
    """
    Check whether (xp,yp) is left of line segment ((x0,y0) to (x1,y1))
    returns:  >0 if left of line, 0 if on line, <0 if right of line
    Source: http://geomalgorithms.com/a03-_inclusion.html
    """

    return (x1-x0) * (yp-y0) - (xp-x0) * (y1-y0) 

def is_inside(xp, yp, x_set, y_set):
    """
    Given location (xp,yp) and set of line segments (x_set, y_set), determine
    whether (xp,yp) is inside polygon. 
    Source: http://geomalgorithms.com/a03-_inclusion.html
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






if __name__ == "__main__":
    """ Test test """
    xysize = 1600
    ijtot  = 16
    dxy    = xysize/ijtot
   
    # Simplified grid at full and half (h) locations 
    x  = np.arange(0.5*dxy, ijtot*dxy, dxy)
    y  = np.arange(0.5*dxy, ijtot*dxy, dxy)
    xh = np.arange(0, ijtot*dxy+0.1, dxy)
    yh = np.arange(0, ijtot*dxy+0.1, dxy)

    # Create some buildings... 
    set_x1 = np.array([300, 1000, 800, 1300, 400,  300])
    set_y1 = np.array([300, 500,  700, 1400, 1000, 300])
    set_x2 = np.array([1200, 1400, 1400, 1200, 1200])
    set_y2 = np.array([300,  300,  500,  500,  300])

    # Plot, and set ticks equal to grid half locations
    pl.figure()
    ax = pl.subplot(111)
    ax.set_aspect(1)
    ax.set_xticks(xh)
    ax.set_yticks(yh)
    ax.set_xlim(0,xysize)
    ax.set_ylim(0,xysize)
    ax.grid(linestyle='dotted')

    # Loop over the different buildings:
    for set_x,set_y in zip([set_x1,set_x2],[set_y1,set_y2]):
       
        # Plot the building
        pl.plot   (set_x, set_y, 'r-', linewidth=1.5, dashes=[4,2])
        pl.scatter(set_x, set_y, marker='+', s=50, facecolor='r', edgecolor='r')
        
        # For all grid points; check wether it is inside a block
        gp    = np.zeros((ijtot,ijtot), dtype=bool)
        for i in range(1,ijtot-1):
            for j in range(1,ijtot-1):
                if(is_inside(x[i], y[j], set_x, set_y)):
                    gp[i,j] = True
      
        # Loop over the points which are inside a block 
        ij_in = np.where(gp)
        for n in range(ij_in[0].size):
            i = ij_in[0][n]
            j = ij_in[1][n]

            # Check if it is a ghost cell:
            if (gp[i-1,j] == False or gp[i+1,j] == False or gp[i,j-1] == False or gp[i,j+1] == False):
                # Find nearest location on the wall:
                xb,yb,zb = nearest_on_polygon(x[i], y[j], 0., set_x, set_y, z=100000.)

                pl.scatter(x[i], y[j], s=30, facecolor='r')
                pl.plot([x[i],xb],[y[j],yb], 'k-')
            else:
                pl.scatter(x[i], y[j], s=30, facecolor='g')
