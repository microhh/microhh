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

import matplotlib.pyplot as pl
import numpy as np

# Available in `microhh/python/`.
import helpers as hlp
import microhh_lbc_tools as mlt

from global_settings import domain_name


"""
Vertical grid.
"""
ktot = 144
dz0 = 20
heights = [0, 4000, 10000]
factors = [1.01, 1.02]
vgrid = hlp.Grid_stretched_manual(ktot, dz0, heights, factors)

#ktot = 192
#dz0 = 20
#heights = [0, 4000, 10000]
#factors = [1.005, 1.015]
#vgrid = hlp.Grid_stretched_manual(ktot, dz0, heights, factors)

n_ghost = 3
n_buffer = 5


if domain_name == 'develop':
    """
    Small development domain.
    """

    central_lon = -57.7
    central_lat = 13.3

    npx_inner = 2
    npy_inner = 4
    
    xsize_inner = 64*100
    ysize_inner = 64*100
    
    itot_inner = 64
    jtot_inner = 64

    dx_inner = xsize_inner / itot_inner
    dy_inner = ysize_inner / jtot_inner

    istart_in_parent = 8
    jstart_in_parent = 8

    npx_outer = 2
    npy_outer = 4
    
    itot_outer = 64
    jtot_outer = 64
    
    dx_outer = 3*dx_inner
    dy_outer = 3*dy_inner


    xsize_outer = itot_outer * dx_outer
    ysize_outer = jtot_outer * dy_outer


elif domain_name == 'mip':
    """
    Domain used for MIP results.
    """

    central_lon = -57.7
    central_lat = 13.3

    npx_inner = 48  # = 144/3
    npy_inner = 64  # = 16*192/48
    
    xsize_inner = 500_000
    ysize_inner = 300_000
    
    itot_inner = 4800
    jtot_inner = 2880
    
    dx_inner = xsize_inner / itot_inner
    dy_inner = ysize_inner / jtot_inner

    istart_in_parent = 150
    jstart_in_parent = 150
    
    # Outer domain:
    npx_outer = 48  # = 144/3
    npy_outer = 32  # = 8*192/48
    
    itot_outer = int(itot_inner/2)
    jtot_outer = int(jtot_inner/2)
    
    dx_outer = 4*dx_inner
    dy_outer = 4*dy_inner
    
    xsize_outer = itot_outer * dx_outer
    ysize_outer = jtot_outer * dy_outer
    

else:
    raise Exception('Invalid domain.')


#"""
#Checks!
#"""
#print('Inner, gridpoints/core=', int(itot_inner*jtot_inner*ktot/npx_inner/npy_inner))
#hlp.check_grid_decomposition(itot_inner, jtot_inner, ktot, npx_inner, npy_inner)
#
#print('Out, gridpoints/core=', int(itot_outer*jtot_outer*ktot/npx_outer/npy_outer))
#hlp.check_grid_decomposition(itot_outer, jtot_outer, ktot, npx_outer, npy_outer)


"""
Horizontal projection.
"""
proj_str = '+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs +type=crs'

hgrid_inner = mlt.Projection(
        xsize_inner, ysize_inner,
        itot_inner, jtot_inner,
        central_lon, central_lat,
        anchor='center',
        proj_str=proj_str)

w_spacing = dx_outer * istart_in_parent
s_spacing = dx_outer * jstart_in_parent
sw_lon, sw_lat = hgrid_inner.to_lonlat(-w_spacing, -s_spacing)

hgrid_outer = mlt.Projection(
        xsize_outer, ysize_outer,
        itot_outer, jtot_outer,
        sw_lon, sw_lat,
        anchor='southwest',
        proj_str=proj_str)


"""
Horizontal grid with ghost cells.
"""
def add_ghost_cells(hgrid_in, n_ghost):

    new_xsize = hgrid_in.xsize + 2*n_ghost*hgrid_in.dx
    new_ysize = hgrid_in.ysize + 2*n_ghost*hgrid_in.dy

    new_itot = hgrid_in.itot + 2*n_ghost
    new_jtot = hgrid_in.jtot + 2*n_ghost

    return mlt.Projection(
            new_xsize, new_ysize,
            new_itot, new_jtot,
            hgrid_in.central_lon,
            hgrid_in.central_lat,
            anchor='center',
            proj_str=proj_str)

hgrid_inner_pad = add_ghost_cells(hgrid_inner, n_ghost)
hgrid_outer_pad = add_ghost_cells(hgrid_outer, n_ghost)


#"""
#Plot!
#"""
#pl.close('all')
#
#pl.figure()
#pl.plot(vgrid.dz, vgrid.z, '-x')
#pl.xlabel(r'$\Delta z$ (m)')
#pl.ylabel(r'$z$ (m)')
#
#pl.figure()
#pl.plot(hgrid_inner.bbox_lon, hgrid_inner.bbox_lat, color='tab:red', label=f'$\Delta={hgrid_inner.dx:.1f} m$')
#pl.plot(hgrid_outer.bbox_lon, hgrid_outer.bbox_lat, color='tab:blue', label=f'$\Delta={hgrid_outer.dx:.1f} m$')
#pl.plot(hgrid_inner_pad.bbox_lon, hgrid_inner_pad.bbox_lat, '--', color='tab:red')
#pl.plot(hgrid_outer_pad.bbox_lon, hgrid_outer_pad.bbox_lat, '--', color='tab:blue')
#pl.legend()
