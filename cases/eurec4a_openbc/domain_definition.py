import numpy as np

import ls2d

# Make sure `microhhpy` is in the Python path!
from microhhpy.spatial import Domain, plot_domains, calc_vertical_grid_2nd

"""
Define vertical grid. This is used in several scripts, so define it globally.
"""
#ktot = 144
#dz0 = 20
#heights = [0, 4000, 10000]
#factors = [1.01, 1.02]
##vgrid =hlp.Grid_stretched_manual(ktot, dz0, heights, factors)
#vgrid = ls2d.grid.Grid_stretched_manual(ktot, dz0, heights, factors)
#vgrid.plot()

ktot = 128
dz0 = 20
heights = [0, 3000, 5000, 10000]
factors = [1.01, 1.03, 1.08]
vgrid = ls2d.grid.Grid_stretched_manual(ktot, dz0, heights, factors)
#vgrid.plot()

# Define buffer height globally; needed by multiple scripts.
zstart_buffer = 0.75 * vgrid.zsize

proj_str = '+proj=utm +zone=20 +ellps=WGS84 +towgs84=0,0,0 +units=m +no_defs +type=crs'


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
        lon=-59.432,
        lat=13.165,
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
        n_sponge=3,
        lbc_freq=60,
        parent=outer_dom,
        xstart_in_parent=500,
        ystart_in_parent=500
    )

    outer_dom.child = inner_dom

    inner_dom.npx = 2
    inner_dom.npy = 4

    #plot_domains([outer_dom, inner_dom], use_projection=True)

    return outer_dom, inner_dom


def get_test_domain():

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

    return outer_dom, inner_dom


def get_large_domain():

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

    return outer_dom, inner_dom


def get_production_domain_200m():
    """
    First test nested run; 600 m outer + 200 m inner.
    """

    # npx=  32, npy=  48 -> cores=1536 | 1728 x  960 x  128 | pts/core=138240.0
    # npx=  64, npy=  48 -> cores=3072 | 2688 x 1728 x  128 | pts/core=193536.0

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
        lon=-57.7,
        lat=13.3,
        anchor='center',
        proj_str=proj_str,
    )

    outer_dom = Domain(
        xsize=1728*600,
        ysize=960*600,
        itot=1728,
        jtot=960,
        n_ghost=3,
        n_sponge=6,
        lbc_freq=3600,
        lon=-55.5,
        lat=14.2,
        anchor='center',
        proj_str=proj_str,
    )

    outer_dom.npx = 32
    outer_dom.npy = 48

    inner_dom = Domain(
        xsize=2688*200,
        ysize=1728*200,
        itot=2688,
        jtot=1728,
        n_ghost=3,
        n_sponge=3,
        lbc_freq=60,
        parent=outer_dom,
        xstart_in_parent=25200,
        ystart_in_parent=25200)

    inner_dom.npx = 64
    inner_dom.npy = 48

    outer_dom.child = inner_dom

    #domains = [outer_dom, inner_dom, ref_dom]
    #labels = []
    #labels.append(rf'Outer: {outer_dom.xsize/1000} x {outer_dom.ysize/1000} km$^2$ @ $\Delta$={outer_dom.dx:.0f} m')
    #labels.append(rf'Inner: {inner_dom.xsize/1000} x {inner_dom.ysize/1000} km$^2$ @ $\Delta$={inner_dom.dx:.0f} m')
    #labels.append(rf'MIP-ref: {ref_dom.xsize/1000} x {ref_dom.ysize/1000} km$^2$')
    #plot_domains(domains, use_projection=True, labels=labels)

    return outer_dom, inner_dom


def get_production_domain_100m():
    """
    First test nested run; 300 m outer + 100 m inner.
    """

    # npx=  64, npy=  48 -> cores=3072 | 3264 x 1920 x  128 | pts/core=261120.0
    # npx=  64, npy=  96 -> cores=6144 | 5760 x 3456 x  128 | pts/core=414720.0

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
        n_sponge=6,
        lbc_freq=3600,
        lon=-55.7,
        lat=14.2,
        anchor='center',
        proj_str=proj_str,
    )

    outer_dom.npx = 64
    outer_dom.npy = 48

    inner_dom = Domain(
        xsize=5760*100,
        ysize=3456*100,
        itot=5760,
        jtot=3456,
        n_ghost=3,
        n_sponge=3,
        lbc_freq=60,
        parent=outer_dom,
        xstart_in_parent=21000,
        ystart_in_parent=25200)

    inner_dom.npx = 64
    inner_dom.npy = 96

    outer_dom.child = inner_dom

    #domains = [outer_dom, inner_dom, ref_dom]
    #labels = []
    #labels.append(rf'Outer: {outer_dom.xsize/1000} x {outer_dom.ysize/1000} km$^2$ @ $\Delta$={outer_dom.dx:.0f} m')
    #labels.append(rf'Inner: {inner_dom.xsize/1000} x {inner_dom.ysize/1000} km$^2$ @ $\Delta$={inner_dom.dx:.0f} m')
    #labels.append(rf'MIP-ref: {ref_dom.xsize/1000} x {ref_dom.ysize/1000} km$^2$')
    #plot_domains(domains, use_projection=True, labels=labels)

    return outer_dom, inner_dom


if __name__ == '__main__':
    """
    Just for plotting vertical grid and/or domain.
    """

    outer_dom, inner_dom = get_develop_domain()
    #outer_dom, inner_dom = get_test_domain()
    #outer_dom, inner_dom = get_large_domain()
    #outer_dom, inner_dom = get_production_domain_200m()
    #outer_dom, inner_dom = get_production_domain_100m()
