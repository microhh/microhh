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

import ls2d

# Make sure `microhhpy` is in the Python path!
from microhhpy.spatial import Domain, plot_domains, calc_vertical_grid_2nd, refine_grid_for_nesting


"""
Define vertical grid. This is used in several scripts, so define it globally.
"""
"""
ktot = 144
dz0 = 20
heights = [0, 4000, 10000]
factors = [1.01, 1.02]
vgrid = ls2d.grid.Grid_stretched_manual(ktot, dz0, heights, factors)
vgrid.plot()
"""

ktot = 128
dz0 = 20
heights = [0, 3000, 5000, 10000]
factors = [1.01, 1.03, 1.1]
vgrid = ls2d.grid.Grid_stretched_manual(ktot, dz0, heights, factors)
#vgrid.plot()

gd_outer = calc_vertical_grid_2nd(vgrid.z, vgrid.zsize)
gd_inner = gd_outer


"""
# Define low-resolution vertical grid, and use grid refinement
# to calculate the high-resolution grid.
# This step is necessary to ensure that the half level heights match.
ktot = 64
dz0 = 40
heights = [0, 3000, 5000, 10000]
factors = [1.02, 1.05, 1.1]
vgrid_out = ls2d.grid.Grid_stretched_manual(ktot, dz0, heights, factors)
#vgrid_out.plot()

# Convert (LS)2D grid to MicroHH compatible grid.
gd_outer = calc_vertical_grid_2nd(vgrid_out.z, vgrid_out.zsize)
z,zh = refine_grid_for_nesting(gd_outer['z'], gd_outer['zh'], ratio=2)
gd_inner = calc_vertical_grid_2nd(z, zh[-1])
#gd_inner = gd_outer
"""

# Define buffer height globally; needed by multiple scripts.
zstart_buffer = 0.75 * gd_outer['zsize']


#proj_str = '+proj=utm +zone=20 +ellps=WGS84 +towgs84=0,0,0 +units=m +no_defs +type=crs'
proj_str = '+proj=tmerc +lat_0=13.3 +lon_0=-57.7 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'


def get_develop_domain():
    """
    Test domain over BCO, to allow testing column stats.
    """
    outer_dom = Domain(
        xsize=32000,
        ysize=16000,
        itot=128,
        jtot=64,
        n_ghost=3,
        n_sponge=5,
        lbc_freq=3600,
        lon=-57.7,
        lat=13.3,
        anchor='center',
        proj_str=proj_str,
    )

    outer_dom.npx = 2
    outer_dom.npy = 4

    inner_dom = Domain(
        xsize=16000,
        ysize=8000,
        itot=128,
        jtot=64,
        n_ghost=3,
        n_sponge=0,
        lbc_freq=60,
        parent=outer_dom,
        xstart_in_parent=2000,
        ystart_in_parent=2000
    )

    outer_dom.child = inner_dom

    inner_dom.npx = 2
    inner_dom.npy = 4

    #plot_domains([outer_dom, inner_dom], use_projection=True)

    return outer_dom, inner_dom


def get_full_domain_100m():
    """
    First test nested run; 300 m outer + 100 m inner.
    """

    # Outer:
    # npx=  64, npy=  48 -> cores=3072 | 3264 x 1920 x  128 | pts/core=261120.0
    # or:
    # npx=  64, npy=  96 -> cores=6144 | 3264 x 1920 x  128 | pts/core=130560.0

    # Inner:
    # npx=  64, npy=  96 -> cores=6144 | 5568 x 3264 x  128 | pts/core=378624.0

    #proj_str = '+proj=tmerc +lat_0=13.3 +lon_0=-57.7 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
    proj_str = '+proj=lcc +lat_1=12.0 +lat_2=16.5 +lat_0=14.3 +lon_0=-55.7 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

    # Dummy, just for plotting.
    # This is the domain requested by the intercomparison.
    ref_dom = Domain(
        xsize=500_000,
        ysize=300_000,
        itot=128,
        jtot=64,
        n_ghost=3,
        n_sponge=5,
        lbc_freq=3600,
        lon=-57.7,  # Don't change
        lat=13.3,   # Don't change
        anchor='center',
        proj_str=proj_str,
    )

    outer_dom = Domain(
        xsize=3264*300,
        ysize=1920*300,
        itot=3264,
        jtot=1920,
        n_ghost=3,
        n_sponge=10,
        lbc_freq=3600,
        lon=-55.7,
        lat=14.3,
        anchor='center',
        proj_str=proj_str,
    )

    outer_dom.npx = 64
    outer_dom.npy = 96

    inner_dom = Domain(
        xsize=5568*100,
        ysize=3264*100,
        itot=5568,
        jtot=3264,
        n_ghost=3,
        n_sponge=0,
        lbc_freq=60,
        parent=outer_dom,
        xstart_in_parent=15000,
        ystart_in_parent=18000)

    inner_dom.npx = 64
    inner_dom.npy = 96

    outer_dom.child = inner_dom

#    domains = [outer_dom, inner_dom, ref_dom]
#    labels = []
#    labels.append(rf'Outer: {outer_dom.xsize/1000} x {outer_dom.ysize/1000} km$^2$ @ $\Delta$={outer_dom.dx:.0f} m')
#    labels.append(rf'Inner: {inner_dom.xsize/1000} x {inner_dom.ysize/1000} km$^2$ @ $\Delta$={inner_dom.dx:.0f} m')
#    labels.append(rf'MIP-ref: {ref_dom.xsize/1000} x {ref_dom.ysize/1000} km$^2$')
#    plot_domains(domains, use_projection=True, labels=labels)

    return outer_dom, inner_dom


def get_quarter_domain_100m():
    """
    Half domain size (in x+y, so quarter area) at full resolution.
    """
    # Outer:
    # npx=  32, npy=  24 -> cores= 768 | 1632 x  960 x  128 | pts/core=261120.0
    # or:
    # npx=  32, npy=  48 -> cores=1536 | 1632 x  960 x  128 | pts/core=130560.0

    # Inner:
    # npx=  32, npy=  48 -> cores=1536 | 2784 x 1632 x  128 | pts/core=378624.0

    # Dummy, just for plotting.
    # This is the domain requested by the intercomparison.
    ref_dom = Domain(
        xsize=500_000,
        ysize=300_000,
        itot=128,
        jtot=64,
        n_ghost=3,
        n_sponge=5,
        lbc_freq=3600,
        lon=-57.7,  # Don't change
        lat=13.3,   # Don't change
        anchor='center',
        proj_str=proj_str,
    )

    outer_dom = Domain(
        xsize=1632*300,
        ysize=960*300,
        itot=1632,
        jtot=960,
        n_ghost=3,
        n_sponge=10,
        lbc_freq=3600,
        lon=-57.7,
        lat=13.3,
        anchor='center',
        proj_str=proj_str,
    )

    outer_dom.npx = 32
    outer_dom.npy = 48

    inner_dom = Domain(
        xsize=2784*100,
        ysize=1632*100,
        itot=2784,
        jtot=1632,
        n_ghost=3,
        n_sponge=0,
        lbc_freq=60,
        parent=outer_dom,
        xstart_in_parent=12000,
        ystart_in_parent=12000)

    inner_dom.npx = 32
    inner_dom.npy = 48

    outer_dom.child = inner_dom

#    domains = [outer_dom, inner_dom, ref_dom]
#    labels = []
#    labels.append(rf'Outer: {outer_dom.xsize/1000} x {outer_dom.ysize/1000} km$^2$ @ $\Delta$={outer_dom.dx:.0f} m')
#    labels.append(rf'Inner: {inner_dom.xsize/1000} x {inner_dom.ysize/1000} km$^2$ @ $\Delta$={inner_dom.dx:.0f} m')
#    labels.append(rf'MIP-ref: {ref_dom.xsize/1000} x {ref_dom.ysize/1000} km$^2$')
#    plot_domains(domains, use_projection=True, labels=labels)

    return outer_dom, inner_dom


if __name__ == '__main__':
    """
    Just for plotting vertical grid and/or domain.
    """

    #outer_dom, inner_dom = get_develop_domain()
    outer_dom, inner_dom = get_full_domain_100m()
    #outer_dom, inner_dom = get_quarter_domain_100m()
