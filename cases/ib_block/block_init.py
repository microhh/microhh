import numpy as np
import matplotlib.pylab as pl

pl.close('all')

# Custom functions for IB & MicroHH
import ib_tools as ibt
import microhh_tools as mht

# Read namelist and grid:
nl = mht.Read_namelist()
gr = mht.Create_grid('block.prof')

# Location block
x0 = 3
x1 = 4
y0 = 2.7
y1 = 3.7
z1 = 1

# Coordinates block
x = np.array([x0, x1, x1, x0, x0])
y = np.array([y0, y0, y1, y1, y0])
poly = [{'z':z1, 'x':x, 'y':y}]

# Find all grid points inside IB
inside_u_2d = ibt.flag_inside_ib_2d(gr.xh, gr.y,  poly)
inside_v_2d = ibt.flag_inside_ib_2d(gr.x,  gr.yh, poly)
inside_w_2d = ibt.flag_inside_ib_2d(gr.x,  gr.y,  poly)

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
