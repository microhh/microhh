import numpy as np
import matplotlib.pylab as pl

# Import MicroHH and IB specific tools
import microhh_tools as mht
import ib_tools as ibt
import ib_tools_cython as ibtc

# Read the namelist, and grid info:
nl = mht.Read_namelist()
gr = mht.Create_grid('must.prof')

# Offset inside model domain compared to input location objects:
x_offset = 20
y_offset = 35

# Read the polygon input as a list of dictionaries, with each dict having x,y,z keys:
poly = ibt.read_polygon_unstructured('container_corners.txt', x_offset, y_offset)

#poly = poly[:5]

# Find all grid points inside IB
inside_u_2d = ibtc.flag_inside_ib_2d(gr.xh, gr.y,  poly)
inside_v_2d = ibtc.flag_inside_ib_2d(gr.x,  gr.yh, poly)
inside_w_2d = ibtc.flag_inside_ib_2d(gr.x,  gr.y,  poly)

# Find the lateral ghost cells
ghost_u_2d = ibt.flag_ghost_cells_2d(inside_u_2d)
ghost_v_2d = ibt.flag_ghost_cells_2d(inside_v_2d)
ghost_w_2d = ibt.flag_ghost_cells_2d(inside_w_2d)

# Calculate all ghost cell properties (nearest location wall, interpolation points, ....)
ghost_u = ibt.calc_ghost_cell_properties(ghost_u_2d, inside_u_2d, poly, gr.xh, gr.y,  gr.z,  nl['IB']['n_idw'], 'u')
ghost_v = ibt.calc_ghost_cell_properties(ghost_v_2d, inside_v_2d, poly, gr.x,  gr.yh, gr.z,  nl['IB']['n_idw'], 'v')
ghost_w = ibt.calc_ghost_cell_properties(ghost_w_2d, inside_w_2d, poly, gr.x,  gr.y,  gr.zh, nl['IB']['n_idw'], 'w', kstart=1)
ghost_s = ibt.calc_ghost_cell_properties(ghost_w_2d, inside_w_2d, poly, gr.x,  gr.y,  gr.z,  nl['IB']['n_idw'], 's')

# Some debugging output
print('Found {} u-ghost cells'.format(len(ghost_u)))
print('Found {} v-ghost cells'.format(len(ghost_v)))
print('Found {} w-ghost cells'.format(len(ghost_w)))
print('Found {} s-ghost cells'.format(len(ghost_s)))

# Write input data for MicroHH
ibt.write_output(ghost_u, 'u.ib_input')
ibt.write_output(ghost_v, 'v.ib_input')
ibt.write_output(ghost_w, 'w.ib_input')
ibt.write_output(ghost_s, 's.ib_input')


# Some plottings...
if (False):
    def scatter_obj(arr, xloc, yloc):
        ij_in = np.where(arr != -1)
        for n in range(ij_in[0].size):
            j = ij_in[0][n]
            i = ij_in[1][n]
            pl.scatter(xloc[i], yloc[j], marker='.', s=5, color='g')

    def plot_poly(poly):
        for p in poly:
            pl.plot(p['x'], p['y'], 'k-')

    def plot_ghost(ghost, xloc, yloc):
        for g in ghost:
            if g['k'] == 0:
                i = g['i']
                j = g['j']
                pl.scatter(xloc[i], yloc[j], marker='o', s=15, color='r')
                pl.scatter(g['xi'], g['yi'], marker='o', s=15, color='m')
                pl.scatter(g['xb'], g['yb'], marker='o', s=15, color='b')
                pl.plot([xloc[i],g['xi']], [yloc[j],g['yi']], linewidth=1, color='k', dashes=[4,2])

                for i in range(g['in'].size):
                    pl.scatter(xloc[g['in'][i]], yloc[g['jn'][i]], marker='o', s=15, color='g')
                    pl.plot([g['xi'], xloc[g['in'][i]]], [g['yi'], yloc[g['jn'][i]]], color='k', dashes=[1,1])


    pl.figure()
    pl.subplot(121, aspect='equal')
    pl.title('Inside IB', loc='left')
    scatter_obj(inside_u_2d, gr.xh, gr.y)
    plot_poly(poly)

    pl.subplot(122, aspect='equal')
    pl.title('Ghost cells', loc='left')
    scatter_obj(ghost_u_2d, gr.xh, gr.y)
    plot_poly(poly)

    pl.figure()
    pl.subplot(111, aspect='equal')
    plot_poly(poly)
    plot_ghost(ghost_u, gr.xh, gr.y)
