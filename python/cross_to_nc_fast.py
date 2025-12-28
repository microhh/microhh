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

import os
import microhh_tools as mht  # available in microhh/python directory
import argparse
import collections
import glob
import numpy as np
from multiprocessing import Pool, set_start_method
import platform

if platform.system() == 'Darwin':
    set_start_method('fork')

def convert_to_nc(tasks):

    # Loop over the (variable, mode) combinations
    for variable, mode in tasks:
        print(f'Converting {variable}-{mode}...')

        try:
            iotime = int(round(starttime / 10**iotimeprec))

            # Check if variable is surface cross-section.
            if os.path.isfile("{0}.xy.000.{1:07d}".format(variable, iotime)):
                if mode != 'xy':
                    continue
                at_surface = True
            else:
                at_surface = False

            # Get the half level (000=full, 100=u, 010=v, etc) index.
            filename = "{0}.{1}.nc".format(variable, mode)
            halflevel = '000'

            if not at_surface:
                if indexes is None:
                    indexes_local,halflevel = mht.get_cross_indices(variable, mode)
                else:
                    indexes_local = indexes

                    files = glob.glob("{0:}.{1}.*.{2:05d}.{3:07d}".format(
                            variable, mode, indexes_local[0], iotime))

                    if len(files) == 0:
                        raise Exception('Cannot find any cross-section')

                    halflevel = files[0].split('.')[-3]

            dim = collections.OrderedDict()
            dim['time'] = []
            dim['z'] = range(ktot)
            dim['y'] = range(jtot)
            dim['x'] = range(itot)

            if at_surface:
                dim.pop('z')
                indexes_local = [-1]
            elif mode == 'xy':
                dim.update({'z': []})
            elif mode == 'xz':
                dim.update({'y': []})
            elif mode == 'yz':
                dim.update({'x': []})

            if halflevel[0] == '1':
                dim['xh'] = dim.pop('x')
            if halflevel[1] == '1':
                dim['yh'] = dim.pop('y')
            if halflevel[2] == '1':
                dim['zh'] = dim.pop('z')

            ncfile = mht.Create_ncfile(
                grid, filename, variable, dim, precision, compression)

            for key, val in dim.items():
                if key == 'time':
                    continue
                elif val == []:
                    ncfile.dimvar[key][:] = grid.dim[key][indexes_local]

            for t, time in enumerate(np.arange(starttime, endtime + sampletime, sampletime)):
                for k in range(len(indexes_local)):
                    index = indexes_local[k]
                    otime = int(
                        round(
                            (time) / 10**iotimeprec))
                    if at_surface:
                        f_in = "{0}.{1}.{2}.{3:07d}".format(
                            variable, mode, halflevel, otime)
                    else:
                        f_in = "{0:}.{1}.{2}.{3:05d}.{4:07d}".format(
                            variable, mode, halflevel, index, otime)
                    try:
                        fin = np.fromfile(f_in, dtype)
                    except Exception as ex:
                        #print (ex)
                        break

                    #print(
                    #    "Processing %8s, time=%7i, index=%4i" %
                    #    (variable, otime, index))

                    ncfile.dimvar['time'][t] = time

                    if at_surface:
                        ncfile.var[t, :, :] = fin
                    elif mode == 'xy':
                        ncfile.var[t, k, :, :] = fin
                    elif mode == 'xz':
                        ncfile.var[t, :, k, :] = fin
                    elif mode == 'yz':
                        ncfile.var[t, :, :, k] = fin

            ncfile.close()

        except Exception as ex:
            print(ex)
            print("Failed to create %s" % filename)


"""
Parse command line arguments.
"""
cross_modes = ['xy', 'xz', 'yz']
parser = argparse.ArgumentParser(
    description='Convert MicroHH binary cross-sections to netCDF4 files.')

parser.add_argument(
    '-m',
    '--modes',
    nargs='*',
    help='mode of the cross section',
    choices=cross_modes)

parser.add_argument('-f', '--filename', help='ini file name')
parser.add_argument('-d', '--directory', help='directory')
parser.add_argument('-v', '--vars', nargs='*', help='variable names')
parser.add_argument('-x', '--index', nargs='*', help='indices', type=int)
parser.add_argument('-t0', '--starttime', help='first time step to be parsed')
parser.add_argument('-t1', '--endtime', help='last time step to be parsed')

parser.add_argument(
    '-tstep',
    '--sampletime',
    help='time interval to be parsed')

parser.add_argument(
    '-p',
    '--precision',
    help='precision',
    choices=[
        'single',
         'double'])

parser.add_argument(
    '-n',
    '--nprocs',
    help='Number of processes',
    type=int,
    default=1)

parser.add_argument(
    '-c',
    '--nocompression',
    help='do not compress the netcdf file',
    action='store_true')

parser.add_argument(
    '-o',
    '--order',
    help='order',
    choices=[
        2, 4], type = int)


"""
Read .ini file to set missing arguments to default values.
"""
args = parser.parse_args()

if args.directory is not None:
    os.chdir(args.directory)

modes = args.modes
indexes = args.index

nl = mht.Read_namelist(args.filename)
itot = nl['grid']['itot']
jtot = nl['grid']['jtot']
ktot = nl['grid']['ktot']

def default(arg, default):
    return arg if arg is not None else default

starttime  = int(default(args.starttime,  nl['time']['starttime']))
endtime    = int(default(args.endtime,    nl['time']['endtime']))
sampletime = int(default(args.sampletime, nl['cross']['sampletime']))

variables  = default(args.vars, nl['cross']['crosslist'])
nprocs     = default(args.nprocs, len(variables))
order      = default(args.order, nl['grid'].get('swspatialorder', 2))

# In case variables is a single string, convert to list.
variables = [variables] if not isinstance(variables, list) else variables

if args.modes is None:
    modes = list(nl['cross'].keys() & cross_modes)

    # Check if there are paths in the cross-list
    if 'xy' not in modes:
        for v in nl['cross']['crosslist']:
            if 'path' in v:
                modes.append('xy')
                break
else:
    modes = args.modes

iotimeprec = nl['time'].get('iotimeprec', 0)
precision = args.precision
compression = not(args.nocompression)

# End option parsing
grid = mht.Read_grid(itot, jtot, ktot, order=order)
dtype = np.dtype(grid.en + grid.prec)

# Create all (variable, mode) combinations to parallelize over
tasks = [(var, mode) for var in variables for mode in modes]
chunks = [tasks[i::nprocs] for i in range(nprocs)]

pool = Pool(processes=nprocs)
pool.imap_unordered(convert_to_nc, chunks)
pool.close()
pool.join()
