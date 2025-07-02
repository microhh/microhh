import numpy as np

import ls2d

# Make sure `microhhpy` is in the Python path!
from microhhpy.spatial import Domain, plot_domains, calc_vertical_grid_2nd

# Help scripts from current directory.
import global_settings
import helpers as hlp

"""
Define vertical grid. This is used in several scripts, so define it globally.
"""
ktot = 144
dz0 = 20
heights = [0, 4000, 10000]
factors = [1.01, 1.02]
vgrid =hlp.Grid_stretched_manual(ktot, dz0, heights, factors)
#vgrid = ls2d.grid.Grid_stretched_manual(ktot, dz0, heights, factors)

#ktot = 128
#dz0 = 20
#heights = [0, 3000, 5000, 10000]
#factors = [1.01, 1.03, 1.08]
#vgrid = ls2d.grid.Grid_stretched_manual(ktot, dz0, heights, factors)

# Define buffer height globally; needed by multiple scripts.
zstart_buffer = 0.75 * vgrid.zsize


"""
Define horizontal grid with projection.
"""
proj_str = '+proj=utm +zone=20 +ellps=WGS84 +towgs84=0,0,0 +units=m +no_defs +type=crs'

if global_settings.case == 'develop':

    outer_dom = Domain(
        xsize=32000,
        ysize=16000,
        itot=128,
        jtot=64,
        n_ghost=3,
        n_sponge=5,
        lon=-57.7,
        lat=13.3,
        anchor='center',
        proj_str=proj_str,
    )

    outer_dom.npx = 2
    outer_dom.npy = 4

    inner_dom = None

    #inner_dom = Domain(
    #    xsize=3200,
    #    ysize=1600,
    #    itot=128,
    #    jtot=64,
    #    n_ghost=3,
    #    n_sponge=5,
    #    parent=outer_dom,
    #    center_in_parent=True
    #)


elif global_settings.case == 'test':

    # Snellius test #1.
    #outer_dom = Domain(
    #    xsize=1536*150,
    #    ysize=768*150,
    #    itot=1536,
    #    jtot=768,
    #    n_ghost=3,
    #    n_sponge=5,
    #    lon=-57.7,
    #    lat=13.3,
    #    anchor='center',
    #    proj_str=proj_str,
    #)

    # Snellius test #2.
    # Full domain, ~200 m resolution.
    outer_dom = Domain(
        xsize=500_000,
        ysize=300_000,
        itot=2592,
        jtot=1536,
        n_ghost=3,
        n_sponge=5,
        lon=-57.7,
        lat=13.3,
        anchor='center',
        proj_str=proj_str,
    )

    outer_dom.npx = 48
    outer_dom.npy = 96

    inner_dom = None

elif global_settings.case == 'large':

    # Dummy, just for plotting.
    ref_dom = Domain(
        xsize=500_000,
        ysize=300_000,
        itot=128,
        jtot=64,
        n_ghost=3,
        n_sponge=5,
        lon=-57.7,
        lat=13.3,
        anchor='center',
        proj_str=proj_str,
    )

    outer_dom = Domain(
        xsize=1_000_000,
        ysize=600_000,
        itot=1920,
        jtot=1152,
        n_ghost=3,
        n_sponge=5,
        lon=-55.7,
        lat=14.3,
        anchor='center',
        proj_str=proj_str,
    )

    outer_dom.npx = 48
    outer_dom.npy = 64

    inner_dom = None


#plot_domains([outer_dom, ref_dom], use_projection=True)
