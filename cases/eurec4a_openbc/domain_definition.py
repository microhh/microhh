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


"""
Define horizontal grid with projection.
"""
if global_settings.case == 'develop':

    outer_dom = Domain(
        xsize=6400,
        ysize=3200,
        itot=64,
        jtot=32,
        n_ghost=3,
        n_sponge=5,
        lon=-57.7,
        lat=13.3,
        anchor='center',
        proj_str='+proj=utm +zone=20 +datum=WGS84 +units=m +no_defs',
    )

    inner_dom = Domain(
        xsize=3200,
        ysize=1600,
        itot=64,
        jtot=32,
        n_ghost=3,
        n_sponge=5,
        parent=outer_dom,
        center_in_parent=True
    )

#plot_domains([outer_proj, inner_proj])