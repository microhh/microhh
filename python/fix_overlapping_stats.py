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
"""
Fix overlapping time in MicroHH statistics, when simulation hit/crashed on 
wall-clock time, and was restarted at an earlier time.
Requires NCO (`ncks` and `ncrcat` commands).

Example usage:
    `python fix_overlapping_stats.py drycblles.default.0000000.nc drycblles.default.0003600.nc`
    `python fix_overlapping_stats.py drycblles.default.*.nc`

This results in a file:
    `drycblles.default.nc`
in which all files are merged and the overlapping times are removed.

NOTE: none of the original statistics files are modified/removed, so this is a safe operation.
"""

import xarray as xr
import numpy as np
import subprocess
import datetime
import shutil
import sys
import os

if len(sys.argv) < 3:
    raise Exception('Provide at least two statistics files as arguments!')

files = sys.argv[1:]
n_files = len(files)

print(f'Merging {n_files} files...')

# Output file name, e.g. case_name.default.nc
name = files[0].split('.')[0]
mode = files[0].split('.')[1]
file_out = f'{name}.{mode}.nc'

# List with NetCDF files with duplicate times removed.
files_fixed = []
to_remove = []

for n in range(n_files-1):
    print(f'Fixing {files[n]} + {files[n+1]}...')

    # Read NetCDF files to determine overlapping times.
    nc1 = xr.open_dataset(files[n],   decode_times=False)
    nc2 = xr.open_dataset(files[n+1], decode_times=False)

    i0 = int(np.abs(nc2.time - nc1.time[-1]).argmin())+1
    i1 = nc2.time.size-1

    if i0 > 1:
        # Fix second NetCDF file use ncks. Replace with Xarray?
        file_fixed = f'{name}.{mode}_{n+1}.nc'
        subprocess.run(['ncks', '-d', f'time,{i0},{i1}', files[n+1], file_fixed])
        files_fixed.append(file_fixed)
        to_remove.append(file_fixed)
    else:
        files_fixed.append(files[n+1])

print(files_fixed)

# Merge statistics files with `ncrcat`.
subprocess.run(['ncrcat', files[0], *files_fixed, file_out])

# Check.
nc = xr.open_dataset(file_out, decode_times=False)
dt = np.diff(nc.time.values)

if np.all(np.isclose(dt, dt[0])):
    print('Merge succes...!')
else:
    print('Merge failed...!')

# Cleanup temporary files.
for f in to_remove:
    os.remove(f)
