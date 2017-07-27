import matplotlib
matplotlib.use('TkAgg')

import numpy as np
import matplotlib.pylab as pl

pl.close('all')

# Custom functions for IB & MicroHH
import ib_tools as it
import microhh_tools as mht

def is_inside_block(x, y, x0, x1, y0, y1):
    if x > x0 and x < x1 and y > y0 and y < y1:
        return True
    else:
        return False

# Read namelist and grid:
nl = mht.Read_namelist()
gr = mht.Read_grid(nl['grid']['itot'], nl['grid']['jtot'], nl['grid']['ktot'], nl['grid']['zsize'])

n_idw = 5 # number of interpolation points in inverse distance weighted interpolation

# Find correct (approximate) location of block with size 1x1
#i0 = np.abs(gr.xh - 3.0).argmin()
#i1 = np.abs(gr.xh - 4.0).argmin()
#j0 = np.abs(gr.yh - 2.7).argmin()
#j1 = np.abs(gr.yh - 3.7).argmin()
#
## Define boundaries block:
#x0 = gr.xh[i0] + 0.05 * gr.dx
#x1 = gr.xh[i1] - 0.05 * gr.dx
#y0 = gr.yh[j0] + 0.05 * gr.dy
#y1 = gr.yh[j1] - 0.05 * gr.dy

x0 = 3
x1 = 4
y0 = 2.7
y1 = 3.7

print('block size = {0:.2f} x {1:.2f}, x0={2:.2f}, x1={3:.2f}, y0={4:.2f}, y1={5:.2f}'.format(x1-x0, y1-y0, x0, x1, y0, y1))

xblock = np.array([x0, x1, x1, x0, x0])
yblock = np.array([y0, y0, y1, y1, y0])
zblock = 1

# Write to namelist
mht.replace_namelist_value('x0_block', x0)
mht.replace_namelist_value('x1_block', x1)
mht.replace_namelist_value('y0_block', y0)
mht.replace_namelist_value('y1_block', y1)
mht.replace_namelist_value('z1_block', zblock)

def find_ghost_cells(x, y, z, kstart=0):
    print('Searching ghost cells.....')

    # Lists holding all the ghost cell info
    ghost_list = []

    # Find all grid points which are inside the block:
    inside_block_2d = np.zeros((nl['grid']['itot'], nl['grid']['jtot']), dtype=bool)
    for i in range(nl['grid']['itot']):
        for j in range(nl['grid']['jtot']):
            inside_block_2d[i,j] = is_inside_block(x[i], y[j], x0, x1, y0, y1)

    # Find first grid level above block
    kmax_block = np.abs(z-zblock).argmin()
    if z[kmax_block] < zblock:
        kmax_block += 1

    # Find all the ghost cells, and give them a flag for east/west/... walls
    ghost_cells = np.zeros((nl['grid']['itot'], nl['grid']['jtot'], kmax_block), dtype=int)

    Walls = {'west':1, 'east':2, 'south':4, 'north':8, 'top':16}

    is_ib = np.where(inside_block_2d)
    for n in range(is_ib[0].size):
        i = is_ib[0][n]
        j = is_ib[1][n]

        if not inside_block_2d[i-1,j]:
            ghost_cells[i,j,kstart:kmax_block] += Walls['west']
        if not inside_block_2d[i+1,j]:
            ghost_cells[i,j,kstart:kmax_block] += Walls['east']
        if not inside_block_2d[i,j-1]:
            ghost_cells[i,j,kstart:kmax_block] += Walls['south']
        if not inside_block_2d[i,j+1]:
            ghost_cells[i,j,kstart:kmax_block] += Walls['north']

        ghost_cells[i,j,kmax_block-1] += Walls['top']

    is_ghost = np.where(ghost_cells > 0)
    for n in range(is_ghost[0].size):
        # Indices of ghost cell
        i = is_ghost[0][n]
        j = is_ghost[1][n]
        k = is_ghost[2][n]

        # "Nearest" location on wall (biased...)
        if ghost_cells[i,j,k] & Walls['west'] > 0:
            xb = x0
            yb = y[j]
            zb = z[k]
        elif ghost_cells[i,j,k] & Walls['east'] > 0:
            xb = x1
            yb = y[j]
            zb = z[k]
        elif ghost_cells[i,j,k] & Walls['south'] > 0:
            xb = x[i]
            yb = y0
            zb = z[k]
        elif ghost_cells[i,j,k] & Walls['north'] > 0:
            xb = x[i]
            yb = y1
            zb = z[k]
        elif ghost_cells[i,j,k] & Walls['top'] > 0:
            xb = x[i]
            yb = y[j]
            zb = zblock

        # Interpolation location in fluid domain
        xi = 2*xb - x[i];
        yi = 2*yb - y[j];
        zi = 2*zb - z[k];

        # Nearest grid point to interpolation location:
        ii = np.abs(x-xi).argmin()
        jj = np.abs(x-yi).argmin()
        kk = np.abs(z-zi).argmin()

        # Find n nearest fluid points (outside IB)
        i_n = []
        j_n = []
        k_n = []
        d   = []
        for di in range(-2,3):
            for dj in range(-2,3):
                for dk in range(-5,6):       # more points in the vertical because of stretched grid!
                    if (not inside_block_2d[ii+di,jj+dj] or z[kk+dk] > zblock):
                        dist = it.abs_dist(xi, yi, zi, x[ii+di], y[jj+dj], z[kk+dk])
                        i_n.append(ii+di)
                        j_n.append(jj+dj)
                        k_n.append(kk+dk)
                        d  .append(dist)

        # Sort on distance, and clip to the required number of interpolation points
        inds = np.array(d).argsort()
        d    = np.array(d)  [inds][:n_idw]
        i_n  = np.array(i_n)[inds][:n_idw]
        j_n  = np.array(j_n)[inds][:n_idw]
        k_n  = np.array(k_n)[inds][:n_idw]

        # Add to dict with ghost cells
        ghost_list.append({'i':i, 'j':j, 'k':k,
                           'xb':xb, 'yb':yb, 'zb':zb,
                           'xi':xi, 'yi':yi, 'zi':zi,
                           'in':i_n, 'jn':j_n, 'kn':k_n})

    return ghost_list


    ## Empty list which holds all the ghost cell info
    #ghost_cells = []
    #count = 0
    #
    ## Loop over the points which are inside a block
    #ij_in = np.where(gp)
    #for n in range(ij_in[0].size):
    #    i = ij_in[0][n]
    #    j = ij_in[1][n]

    #    for k in range(kmax_block):
    #        # Check if it is a ghost cell:
    #        if (gp[i-1,j] == False or gp[i+1,j] == False or gp[i,j-1] == False or gp[i,j+1] == False or k == kmax_block-1):
    #
    #            # Find nearest location on the wall:
    #            xb,yb,zb = it.nearest_on_polygon(x[i], y[j], z[k], xblock, yblock, z=zblock)

    #            # Interpolation location:
    #            xi = 2*xb - x[i];
    #            yi = 2*yb - y[j];
    #            zi = 2*zb - z[k];

    #            # Nearest grid point to interpolation location:
    #            ii = np.abs(x-xi).argmin()
    #            jj = np.abs(x-yi).argmin()
    #            kk = np.abs(z-zi).argmin()

    #            # Find n nearest fluid points (outside IB)
    #            i_n = []
    #            j_n = []
    #            k_n = []
    #            d   = []
    #            for di in range(-3,4):
    #                for dj in range(-3,4):
    #                    for dk in range(-5,6):       # more points in the vertical because of stretched grid!
    #                        if (gp[ii+di,jj+dj] == False or z[kk+dk] > zblock):
    #                            dist = it.abs_dist(xi, yi, zi, x[ii+di], y[jj+dj], z[kk+dk])
    #                            i_n.append(ii+di)
    #                            j_n.append(jj+dj)
    #                            k_n.append(kk+dk)
    #                            d  .append(dist)
    #
    #            # Sort on distance
    #            inds = np.array(d).argsort()
    #            d    = np.array(d)  [inds][:n_idw]
    #            i_n  = np.array(i_n)[inds][:n_idw]
    #            j_n  = np.array(j_n)[inds][:n_idw]
    #            k_n  = np.array(k_n)[inds][:n_idw]

    #            # Add to list of ghost cells
    #            ghost_cells.append({'i':i, 'j':j, 'k':k,
    #                                'xb':xb, 'yb':yb, 'zb':zb,
    #                                'xi':xi, 'yi':yi, 'zi':zi,
    #                                'in':i_n, 'jn':j_n, 'kn':k_n})

    #            #if k == 42:
    #            #    print(count)
    #            count += 1

    #return ghost_cells


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

# Skip if script is executed with -i in IPython
if "ghost_cells_u" not in locals():
    ghost_cells_u = find_ghost_cells(gr.xh, gr.y,  gr.z )
    ghost_cells_v = find_ghost_cells(gr.x,  gr.yh, gr.z )
    ghost_cells_w = find_ghost_cells(gr.x,  gr.y,  gr.zh, kstart=1)
    ghost_cells_s = find_ghost_cells(gr.x,  gr.y,  gr.z)

    write_output(ghost_cells_u, 'u.ib_input')
    write_output(ghost_cells_v, 'v.ib_input')
    write_output(ghost_cells_w, 'w.ib_input')
    write_output(ghost_cells_s, 's.ib_input')


if (False):
    from mpl_toolkits.mplot3d import Axes3D

    x = gr.xh
    y = gr.y
    z = gr.z
    ghost_cells = ghost_cells_u

    pl.figure()
    ax=pl.subplot(111, projection='3d')
    ax.set_xticks(gr.xh)
    ax.set_yticks(gr.yh)
    ax.set_zticks(gr.zh)

    for zz in [0,zblock]:
        ax.plot3D([x0,x1],[y0,y0],[zz,zz],'k-')
        ax.plot3D([x0,x1],[y1,y1],[zz,zz],'k-')
        ax.plot3D([x0,x0],[y0,y1],[zz,zz],'k-')
        ax.plot3D([x1,x1],[y0,y1],[zz,zz],'k-')
    ax.plot3D([x0,x0],[y0,y0],[0,zz],'k-')
    ax.plot3D([x1,x1],[y0,y0],[0,zz],'k-')
    ax.plot3D([x0,x0],[y1,y1],[0,zz],'k-')
    ax.plot3D([x1,x1],[y1,y1],[0,zz],'k-')

    # Plot ghost cells
    ii = np.random.randint(low=0, high=len(ghost_cells), size=5)

    for i in ii:
        g  = ghost_cells[i]
        xg = x[g['i']]
        yg = y[g['j']]
        zg = z[g['k']]
        ax.scatter(xg, yg, zg, s=30, facecolor='r')
        pl.plot([g['xi'],xg], [g['yi'],yg], [g['zi'],zg], 'k-')

        for n in range(g['in'].size):
            xn = x[g['in'][n]]
            yn = y[g['jn'][n]]
            zn = z[g['kn'][n]]
            ax.plot([g['xi'],xn],[g['yi'],yn],[g['zi'],zn], 'r:')
            ax.scatter(xn, yn, zn, marker='x')
