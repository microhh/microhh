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

from datetime import datetime
import numpy as np
import ls2d

from microhhpy.spatial import Domain, plot_domains, calc_vertical_grid_2nd

float_type = np.float64

# Data paths.
#syst = 'eddy'
syst = 'snellius'

if syst == 'eddy':
    env = dict(
        work_path = 'develop_case',
        microhh_path = '/home/bart/meteo/models/microhh/',
        microhh_bin = '/home/bart/meteo/models/microhh/build_dp_cpumpi/microhh',
        gpt_veerman_path = '/home/bart/meteo/models/coefficients_veerman',
        cosmo_path = '/home/scratch2/bart/eurec4a_cosmo/',
        ls2d_era5_path = '/home/scratch1/bart/LS2D_ERA5',
        ls2d_cams_path = '/home/scratch1/bart/LS2D_CAMS')

elif  syst == 'snellius':
    env = dict(
        work_path = '/gpfs/work2/0/nwo21036/bart/eurec4a_full',
        work_path2 = '/scratch-shared/stratum2/eurec4a_full',
        microhh_path = '/home/stratum2/meteo/models/microhh',
        microhh_bin = '/home/stratum2/meteo/models/microhh/build_dp_cpumpi',
        gpt_veerman_path = '/home/stratum2/meteo/models/coefficients_veerman',
        cosmo_path = '/gpfs/work3/0/lesmodels/eurec4a',
        ls2d_era5_path = '/gpfs/work3/0/lesmodels/ls2d_era5',
        ls2d_cams_path = '/gpfs/work3/0/lesmodels/ls2d_cams')

host_model = 'COSMO'   # COSMO or ERA5


"""
Define simulation period.
COSMO has data from 2020-02-01 00:00:00 to 2020-02-12 00:00:00.
"""
start_date = datetime(year=2020, month=2, day=1,  hour=0)
end_date   = datetime(year=2020, month=2, day=12, hour=0)


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
_vgrid = ls2d.grid.Grid_stretched_manual(ktot, dz0, heights, factors)
#vgrid.plot()

# Vertical grid consistent with MicroHH definition.
gd = calc_vertical_grid_2nd(_vgrid.z, _vgrid.zsize)

zstart_buffer = 0.75 * gd['zsize']

#proj_str = '+proj=utm +zone=20 +ellps=WGS84 +towgs84=0,0,0 +units=m +no_defs +type=crs'
#proj_str = '+proj=tmerc +lat_0=13.3 +lon_0=-57.7 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
proj_str = '+proj=lcc +lat_1=12.0 +lat_2=16.5 +lat_0=14.3 +lon_0=-55.7 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'


# Requested output domain MIP:
ref_domain = Domain(
    xsize=500_000,
    ysize=300_000,
    itot=128,
    jtot=64,
    lon=-57.7,  # Don't change
    lat=13.3,   # Don't change
    anchor='center',
    proj_str=proj_str
)


def get_develop_domain(plot=False):
    """
    Test domain over BCO, to allow testing column stats.
    """

    dom0 = Domain(
        xsize=64*400,
        ysize=32*400,
        itot=64,
        jtot=32,
        n_ghost=3,
        n_sponge=5,
        lbc_freq=3600,
        lon=-57.7,
        lat=13.3,
        anchor='center',
        proj_str=proj_str,
        work_dir=f'{env["work_path"]}/dom0'
    )

    dom0.npx = 2
    dom0.npy = 4

    dom1 = Domain(
        xsize=64*200,
        ysize=32*200,
        itot=64,
        jtot=32,
        n_ghost=3,
        n_sponge=0,
        lbc_freq=60,
        parent=dom0,
        xstart_in_parent=800,
        ystart_in_parent=800,
        work_dir=f'{env["work_path"]}/dom1'
    )

    dom1.npx = 2
    dom1.npy = 4

    dom2 = Domain(
        xsize=64*100,
        ysize=32*100,
        itot=64,
        jtot=32,
        n_ghost=3,
        n_sponge=0,
        lbc_freq=60,
        parent=dom1,
        xstart_in_parent=400,
        ystart_in_parent=400,
        work_dir=f'{env["work_path"]}/dom2'
    )

    dom2.npx = 2
    dom2.npy = 4

    dom0.child = dom1
    dom1.child = dom2

    if plot:
        plot_domains([dom0, dom1, dom2], use_projection=True)

    return [dom0, dom1, dom2]


def get_quarter_domain_100m(plot=False):
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
        work_dir=f'{env["work_path"]}/dom0'
    )

    dom0 = Domain(
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
        work_dir=f'{env["work_path"]}/dom1'
    )

    dom0.npx = 32
    dom0.npy = 48

    dom1 = Domain(
        xsize=2784*100,
        ysize=1632*100,
        itot=2784,
        jtot=1632,
        n_ghost=3,
        n_sponge=0,
        lbc_freq=60,
        parent=dom0,
        xstart_in_parent=12000,
        ystart_in_parent=12000,
        work_dir=f'{env["work_path"]}/dom2'
    )

    dom1.npx = 32
    dom1.npy = 48

    dom0.child = dom1

    if plot:
        domains = [dom0, dom1, ref_dom]
        labels = []
        labels.append(rf'Outer: {dom0.xsize/1000} x {dom0.ysize/1000} km$^2$ @ $\Delta$={dom0.dx:.0f} m')
        labels.append(rf'Inner: {dom1.xsize/1000} x {dom1.ysize/1000} km$^2$ @ $\Delta$={dom1.dx:.0f} m')
        labels.append(rf'MIP-ref: {ref_dom.xsize/1000} x {ref_dom.ysize/1000} km$^2$')
        plot_domains(domains, use_projection=True, labels=labels)

    return dom0, dom1, None


def get_full_domain_300_100(plot=False):
    """
    First test nested run; 300 m outer + 100 m inner.
    """

    # Outer:
    # npx=  64, npy=  48 -> cores=3072 | 3264 x 1920 x  128 | pts/core=261120.0
    # or:
    # npx=  64, npy=  96 -> cores=6144 | 3264 x 1920 x  128 | pts/core=130560.0

    # Inner:
    # npx=  64, npy=  96 -> cores=6144 | 5568 x 3264 x  128 | pts/core=378624.0

    # Dummy, just for plotting.
    # This is the domain requested by the intercomparison.
    ref_dom = Domain(
        xsize=500_000,
        ysize=300_000,
        itot=128,
        jtot=64,
        lon=-57.7,  # Don't change
        lat=13.3,   # Don't change
        anchor='center',
        proj_str=proj_str
    )

    dom0 = Domain(
        xsize=3264*300,
        ysize=1920*300,
        itot=3264,
        jtot=1920,
        n_ghost=3,
        n_sponge=10,
        lbc_freq=3600,
        lon=-55.7,
        lat=14.35,
        anchor='center',
        proj_str=proj_str,
        work_dir=f'{env["work_path"]}/dom0'
    )

    dom0.npx = 64
    dom0.npy = 96

    dom1 = Domain(
        xsize=5568*100,
        ysize=3264*100,
        itot=5568,
        jtot=3264,
        n_ghost=3,
        n_sponge=0,
        lbc_freq=60,
        parent=dom0,
        xstart_in_parent=16800,
        ystart_in_parent=16800,
        work_dir=f'{env["work_path"]}/dom1'
    )

    dom1.npx = 64
    dom1.npy = 96

    d33 = 100 / 3.

    dom2 = Domain(
        xsize=1536*d33,
        ysize=768*d33,
        itot=1536,
        jtot=768,
        n_ghost=3,
        n_sponge=0,
        lbc_freq=60,
        parent=dom1,
        xstart_in_parent=59600,
        ystart_in_parent=132400,
        work_dir=f'{env["work_path"]}/dom2'
    )

    dom2.npx = 64
    dom2.npy = 96

    dom0.child = dom1
    dom1.child = dom2

    if plot:
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs

        domains = [dom0, dom1, dom2, ref_dom]
        labels = []
        labels.append(rf'Outer: {dom0.xsize/1000} x {dom0.ysize/1000} km$^2$ @ $\Delta$={dom0.dx:.0f} m')
        labels.append(rf'Inner: {dom1.xsize/1000} x {dom1.ysize/1000} km$^2$ @ $\Delta$={dom1.dx:.0f} m')
        labels.append(rf'HR: {dom2.xsize/1000} x {dom2.ysize/1000} km$^2$ @ $\Delta$={dom2.dx:.0f} m')
        labels.append(rf'MIP-ref: {ref_dom.xsize/1000} x {ref_dom.ysize/1000} km$^2$')
        plot_domains(domains, use_projection=True, labels=labels)
        plt.scatter(-59.432, 13.165, marker='+', transform=ccrs.PlateCarree())

    return dom0, dom1, None



def get_full_domain_400_100(plot=False):
    """
    400 m outer + 100 m inner.
    """

    # Outer:
    # npx=  64, npy=  48 -> cores=3072 | 3264 x 1920 x  128 | pts/core=261120.0
    # or:
    # npx=  64, npy=  96 -> cores=6144 | 3264 x 1920 x  128 | pts/core=130560.0

    # Inner:
    # npx=  64, npy=  96 -> cores=6144 | 5568 x 3264 x  128 | pts/core=378624.0

    # Dummy, just for plotting.
    # This is the domain requested by the intercomparison.
    ref_dom = Domain(
        xsize=500_000,
        ysize=300_000,
        itot=128,
        jtot=64,
        lon=-57.7,  # Don't change
        lat=13.3,   # Don't change
        anchor='center',
        proj_str=proj_str
    )

    dom0 = Domain(
        xsize=2496*400,
        ysize=1536*400,
        itot=2496,
        jtot=1536,
        n_ghost=3,
        n_sponge=10,
        lbc_freq=3600,
        lon=-55.6,
        lat=14.55,
        anchor='center',
        proj_str=proj_str,
        work_dir=f'{env["work_path"]}/dom0'
    )

    dom0.npx = 64
    dom0.npy = 96

    dom1 = Domain(
        xsize=5568*100,
        ysize=3264*100,
        itot=5568,
        jtot=3264,
        n_ghost=3,
        n_sponge=0,
        lbc_freq=60,
        parent=dom0,
        xstart_in_parent=16800,
        ystart_in_parent=15200,
        work_dir=f'{env["work_path"]}/dom1'
    )

    dom1.npx = 64
    dom1.npy = 96

    d33 = 100 / 3.

    dom2 = Domain(
        xsize=1536*d33,
        ysize=768*d33,
        itot=1536,
        jtot=768,
        n_ghost=3,
        n_sponge=0,
        lbc_freq=60,
        parent=dom1,
        xstart_in_parent=59600,
        ystart_in_parent=132400,
        work_dir=f'{env["work_path2"]}/dom2'
    )

    dom2.npx = 32
    dom2.npy = 48

    dom0.child = dom1
    dom1.child = dom2

    if plot:
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs

        domains = [dom0, dom1, dom2, ref_dom]
        labels = []
        labels.append(rf'Outer: {dom0.xsize/1000} x {dom0.ysize/1000} km$^2$ @ $\Delta$={dom0.dx:.0f} m')
        labels.append(rf'Inner: {dom1.xsize/1000} x {dom1.ysize/1000} km$^2$ @ $\Delta$={dom1.dx:.0f} m')
        labels.append(rf'HR: {dom2.xsize/1000} x {dom2.ysize/1000} km$^2$ @ $\Delta$={dom2.dx:.0f} m')
        labels.append(rf'MIP-ref: {ref_dom.xsize/1000} x {ref_dom.ysize/1000} km$^2$')
        plot_domains(domains, use_projection=True, labels=labels)
        plt.scatter(-59.432, 13.165, marker='+', transform=ccrs.PlateCarree())

    return dom0, dom1, dom2


def get_domain(name):
    plot = False
    if name == 'develop':
        return get_develop_domain(plot=plot)
    elif name == 'quarter_100m':
        return get_quarter_domain_100m(plot=plot)
    elif name == 'full_300_100':
        return get_full_domain_300_100(plot=plot)
    elif name == 'full_400_100':
        return get_full_domain_400_100(plot=plot)


if __name__ == '__main__':
    """
    Just for testing/plotting.
    """
    #domains = get_develop_domain(plot=False)
    #domains = get_quarter_domain_100m()
    domains = get_full_domain_400_100(plot=True)
