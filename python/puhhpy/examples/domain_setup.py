#
#  MicroHH
#  Copyright (c) 2011-2024 Chiel van Heerwaarden
#  Copyright (c) 2011-2024 Thijs Heus
#  Copyright (c) 2014-2024 Bart van Stratum
#
#  This file is part of MicroHH
#
#  MicroHH is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MicroHH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
#


# Standard library

# Third-party.

# Local library
from puhhpy.spatial import Domain, plot_domains


def domain_setup_without_proj():

    """
    Specify domains.
    """
    d0 = Domain(
         xsize=3200,
         ysize=3200,
         itot=32,
         jtot=32
         )

    d1 = Domain(
        xsize=1600,
        ysize=1600,
        itot=32,
        jtot=32,
        parent=d0,
        center_in_parent=True
        )

    d2 = Domain(
        xsize=800,
        ysize=800,
        itot=32,
        jtot=32,
        parent=d1,
        xstart_in_parent=200,
        ystart_in_parent=200
        )

    d0.child = d1
    d1.child = d2


    """
    Plot domains.
    """
    plot_domains([d0, d1, d2])


def domain_setup_with_proj():

    """
    Specify domains.
    """
    d0 = Domain(
         xsize=128000,
         ysize=128000,
         itot=64,
         jtot=64,
         lon=4.92,
         lat=51.97,
         anchor='center',
         proj_str='+proj=utm +zone=31 +ellps=intl +towgs84=-87,-98,-121,0,0,0,0 +units=m +no_defs +type=crs'
         )

    d1 = Domain(
        xsize=64000,
        ysize=64000,
        itot=64,
        jtot=64,
        parent=d0,
        center_in_parent=True
        )

    d0.child = d1


    """
    Plot domains.
    """
    plot_domains([d0, d1], use_projection=True, scatter_lonlat=True)


    """
    Using the projection.
    Each domain has `.proj` instance that contains the lat/lon coordinates
    (see docstring of `puhhpy.spatial.Projection`), and can do x/y <-> lon/lat transforms.
    """
    x,y = d0.proj.to_xy(4.92, 51.97)
    lon,lat = d0.proj.to_lonlat(x, y)

    print(f'Coordinates ({lon:.2f}, {lat:.2f}) degrees = ({x:.2f}, {y:.2f}) m in LES.')


if __name__ == '__main__':

    # Demonstration without a spatial projection.
    # For idealized cases.
    domain_setup_without_proj()

    # Demonstration with spatial projection.
    # Needed for e.g. nesting in ERA5, generation
    # of land-use maps, etc.
    domain_setup_with_proj()
