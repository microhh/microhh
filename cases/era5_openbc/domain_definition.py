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

d0 = Domain(
     xsize=51200,
     ysize=51200,
     itot=128,
     jtot=128,
     n_ghost=3,
     lon=5.34,
     lat=53.4,
     anchor='center',
     proj_str='+proj=utm +zone=31 +ellps=intl +towgs84=-87,-98,-121,0,0,0,0 +units=m +no_defs +type=crs'
     )

d1 = Domain(
    xsize=25600,
    ysize=25600,
    itot=128,
    jtot=128,
    n_ghost=3,
    parent=d0,
    center_in_parent=True
    )

domains = [d0, d1]

for i in range(len(domains)-1):
    domains[i].child = domains[i+1]


if __name__ == '__main__':

    plot_domains([d0, d1], use_projection=True, scatter_lonlat=False)