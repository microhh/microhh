import numpy as np

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
vgrid = hlp.Grid_stretched_manual(ktot, dz0, heights, factors)

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


if global_settings.case == '500x300':

    outer_dom = Domain(
        xsize=700_000,
        ysize=500_000,
        itot=2400,
        jtot=1440,
        n_ghost=3,
        n_sponge=5,
        lon=-57.7,
        lat=13.3,
        anchor='center',
        proj_str=proj_str,
    )

    #inner_dom = Domain(
    #    xsize=2400*200,
    #    ysize=1600,
    #    itot=64,
    #    jtot=32,
    #    n_ghost=3,
    #    n_sponge=5,
    #    parent=outer_dom,
    #    center_in_parent=True
    #)


#plot_domains([outer_dom], use_projection=True)
