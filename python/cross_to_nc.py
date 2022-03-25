#
#  MicroHH
#  Copyright (c) 2011-2020 Chiel van Heerwaarden
#  Copyright (c) 2011-2020 Thijs Heus
#  Copyright (c) 2014-2020 Bart van Stratum
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

def convert_to_nc(variables):
    # Loop over the different variables and crosssections
    for variable in variables:
        for mode in modes:
            try:
                otime = int(round(starttime / 10**iotimeprec))
                if os.path.isfile("{0}.xy.{1:07d}".format(variable, otime)):
                    if mode != 'xy':
                        continue
                    at_surface = True
                else:
                    at_surface = False

                filename = "{0}.{1}.nc".format(variable, mode)
                if not at_surface:
                    if indexes is None:
                        indexes_local = mht.get_cross_indices(variable, mode)
                    else:
                        indexes_local = indexes

                #dim = {'time' : range(niter), 'z' : range(ktot), 'y' : range(jtot), 'x' : range(itot)}
                dim = collections.OrderedDict()
                dim['time'] = range(niter)
                dim['z'] = range(ktot)
                dim['y'] = range(jtot)
                dim['x'] = range(itot)

                if at_surface:
                    dim.pop('z')
                    n = itot * jtot
                    indexes_local = [-1]
                elif mode == 'xy':
                    dim.update({'z': indexes_local})
                    n = itot * jtot
                elif mode == 'xz':
                    dim.update({'y': indexes_local})
                    n = itot * ktot
                elif mode == 'yz':
                    dim.update({'x': indexes_local})
                    n = ktot * jtot

                if variable == 'u':
                    dim['xh'] = dim.pop('x')
                if variable == 'v':
                    dim['yh'] = dim.pop('y')
                if variable == 'w':
                    dim['zh'] = dim.pop('z')

                ncfile = mht.Create_ncfile(
                    grid, filename, variable, dim, precision, compression)

                for t in range(niter):
                    for k in range(len(indexes_local)):
                        index = indexes_local[k]
                        otime = int(
                            round(
                                (starttime + t * sampletime) / 10**iotimeprec))
                        if at_surface:
                            f_in = "{0}.{1}.{2:07d}".format(
                                variable, mode, otime)
                        else:
                            f_in = "{0:}.{1}.{2:05d}.{3:07d}".format(
                                variable, mode, index, otime)
                        try:
                            fin = mht.Read_binary(grid, f_in)
                        except Exception as ex:
                            print (ex)
                            raise Exception(
                                'Stopping: cannot find file {}'.format(f_in))

                        print(
                            "Processing %8s, time=%7i, index=%4i" %
                            (variable, otime, index))

                        ncfile.dimvar['time'][t] = otime * 10**iotimeprec

                        if at_surface:
                            ncfile.var[t, :, :] = fin.read(n)
                        elif mode == 'xy':
                            ncfile.var[t, k, :, :] = fin.read(n)
                        elif mode == 'xz':
                            ncfile.var[t, :, k, :] = fin.read(n)
                        elif mode == 'yz':
                            ncfile.var[t, :, :, k] = fin.read(n)

                        fin.close()
                ncfile.close()

            except Exception as ex:
                print(ex)
                print("Failed to create %s" % filename)



# Parse command line and namelist options
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

args = parser.parse_args()

if args.directory is not None:
    os.chdir(args.directory)

modes = args.modes
indexes = args.index

nl = mht.Read_namelist(args.filename)
itot = nl['grid']['itot']
jtot = nl['grid']['jtot']
ktot = nl['grid']['ktot']

starttime = float(
    args.starttime) if args.starttime is not None else nl['time']['starttime']
endtime = float(
    args.endtime) if args.endtime is not None else nl['time']['endtime']
sampletime = float(
    args.sampletime) if args.sampletime is not None else nl['cross']['sampletime']

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

if 'iotimeprec' in nl['time']:
    iotimeprec = nl['time']['iotimeprec']
else:
    iotimeprec = 0.

variables = args.vars if args.vars is not None else nl['cross']['crosslist']

# In case variables is a single string, convert to list.
variables = [ variables ] if not isinstance(variables, list) else variables

precision = args.precision
nprocs = args.nprocs if args.nprocs is not None else len(variables)
compression = not(args.nocompression)
# End option parsing

# Calculate the number of iterations
for time in np.arange(starttime, endtime, sampletime):
    otime = int(round(time / 10**iotimeprec))
    if not glob.glob('*.{0:07d}'.format(otime)):
        endtime = time - sampletime
        break
niter = int((endtime - starttime) / sampletime + 1)

grid = mht.Read_grid(itot, jtot, ktot)

chunks = [variables[i::nprocs] for i in range(nprocs)]

pool = Pool(processes=nprocs)

pool.imap_unordered(convert_to_nc, chunks)

pool.close()
pool.join()
