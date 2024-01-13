#
# MicroHH
# Copyright (c) 2011-2023 Chiel van Heerwaarden
# Copyright (c) 2011-2023 Thijs Heus
# Copyright (c) 2014-2023 Bart van Stratum
# 
# This file is part of MicroHH
# 
# MicroHH is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# MicroHH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
#

import numpy as np

# Available in `microhh/python/`.
import helpers as hlp
import microhh_lbc_tools as mlt

ktot = 144

# Inner domain:
npx_in = 48  # = 144/3
npy_in = 64  # = 16*192/48

xsize_in = 500_000
ysize_in = 300_000

itot_in = 4800
jtot_in = 2880

dx_in = xsize_in / itot_in
dy_in = ysize_in / jtot_in

print('In, gridpoints/core=', int(itot_in*jtot_in*ktot/npx_in/npy_in))

hlp.check_grid_decomposition(itot_in, jtot_in, ktot, npx_in, npy_in)

# Outer domain:
npx_out = 48  # = 144/3
npy_out = 32  # = 8*192/48

itot_out = int(itot_in/2)
jtot_out = int(jtot_in/2)

dx_out = 3*dx_in
dy_out = 3*dy_in

xsize_out = itot_out * dx_out
ysize_out = jtot_out * dy_out

print('Out, gridpoints/core=', int(itot_out*jtot_out*ktot/npx_out/npy_out))

hlp.check_grid_decomposition(itot_out, jtot_out, ktot, npx_out, npy_out)

"""
Horizontal projection.
"""
central_lon = -57.7
central_lat = 13.3

hgrid_in = mlt.Projection(
        xsize_in, ysize_in,
        itot_in, jtot_in,
        central_lon, central_lat,
        anchor='center',
        proj_str='+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs +type=crs')


"""
Vertical grid.
"""
heights = [0, 4000, 10000]
factors = [1.01, 1.02]
vgrid = hlp.Grid_stretched_manual(ktot, 20, heights, factors)

"""
import matplotlib.pyplot as pl
pl.close('all')
pl.figure()
pl.plot(vgrid.dz, vgrid.z, '-x')
pl.xlabel(r'$\Delta z$ (m)')
pl.ylabel(r'$z$ (m)')
"""
