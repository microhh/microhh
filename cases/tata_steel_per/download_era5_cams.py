#
# This file is part of LS2D.
#
# Copyright (c) 2017-2025 Wageningen University & Research
# Author: Bart van Stratum (WUR)
#
# LS2D is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LS2D is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with LS2D.  If not, see <http://www.gnu.org/licenses/>.
#


from datetime import datetime
import ls2d

from global_settings import ls2d_settings, cams_eac4_variables, cams_egg4_variables

# Download ERA5 meteorology.
ls2d.download_era5(ls2d_settings, exit_when_waiting=False)

# Download CAMS chemical species.
ls2d.download_cams(ls2d_settings, variables=cams_eac4_variables, grid=0.25)

# CAMS-EGG only available < 2021..
#ls2d.download_cams(ls2d_settings, variables=cams_egg4_variables, grid=0.25)