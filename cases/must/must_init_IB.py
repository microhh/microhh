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
x_offset = 450
y_offset = 150

# Read the polygon input as a list of dictionaries, with each dict having x,y,z keys:
poly = ibt.read_polygon_unstructured('poly_coordinates_0deg.txt', x_offset, y_offset)

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

