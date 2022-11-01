"""
Fix overlapping time in MicroHH statistics, when simulation hit/crashed on 
wall-clock time, and was restarted at an earlier time.
Requires NCO (`ncks` and `ncrcat` commands).

Example usage:
    `python fix_overlapping_stats.py drycblles.default.0000000.nc drycblles.default.0003600.nc`

This results in a file:
    `drycblles.default.0003600.nc_fixed.nc`,
from which the overlapping times with `drycblles.default.0000000.nc` are removed,
and a final merged file called:
    `drycblles.default.nc`.

NOTE: none of the original statistics files are modified/removed,
so this is a safe operation...
"""

import xarray as xr
import numpy as np
import subprocess
import shutil
import sys

f1 = sys.argv[1]
f2 = sys.argv[2]

name = f1.split('.')[0]
mode = f1.split('.')[1]
f2_fixed = '{}_fixed.nc'.format(f2)
out = '{}.{}.nc'.format(name, mode)

# 2. Read NetCDF files to determine overlapping times.
nc1 = xr.open_dataset(f1, decode_times=False)
nc2 = xr.open_dataset(f2, decode_times=False)

i0 = int(np.abs(nc2.time - nc1.time[-1]).argmin())+1
i1 = nc2.time.size-1

# 3. Fix second NetCDF file use ncks.
subprocess.run(['ncks', '-d', 'time,{},{}'.format(i0, i1), f2, f2_fixed])

# 4. Merge statistics files with `ncrcat`.
subprocess.run(['ncrcat', f1, f2_fixed, out])

# 5. Check:
nc3 = xr.open_dataset(out, decode_times=False)
dt = np.diff(nc3.time.values)

if np.all(np.isclose(dt, dt[0])):
    print('Merge succes...!')
else:
    print('Merge failed...!')
