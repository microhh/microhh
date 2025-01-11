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
