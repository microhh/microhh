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

import netCDF4 as nc
import numpy as np
import struct as st
import glob
import re
import subprocess
import importlib
import shutil
import os
import sys
import filecmp
import timeit
import csv
import copy
import datetime
import itertools
import inspect
from copy import deepcopy

# -------------------------
# General help functions
# -------------------------

def _int_or_float_or_str(value):
    """ Helper function: convert a string to int/float/str """
    try:
        if ('.' in value):
            return float(value)
        else:
            return int(float(value))
    except BaseException:
        return value.rstrip()


def _convert_value(value):
    """ Helper function: convert namelist value or list """
    if ',' in value:
        value = value.split(',')
        return [_int_or_float_or_str(val) for val in value]
    else:
        return _int_or_float_or_str(value)


def _find_namelist_file():
    """ Helper function: automatically find the .ini file in the current directory """
    namelist_file = glob.glob('*.ini')
    if len(namelist_file) == 0:
        raise RuntimeError(
            'Can\'t find any .ini files in the current directory!')
    if len(namelist_file) > 1:
        print('There are multiple .ini files:')
        for i,f in enumerate(namelist_file):
            print('{}: {}'.format(i, f))
        print('Which one do you want to use: ', end='')
        i = int(input())
        return namelist_file[i]
    else:
        return namelist_file[0]

# -------------------------
# Classes and functions to read and write MicroHH things
# -------------------------

class Read_namelist:
    """ Reads a MicroHH .ini file to memory
        All available variables are accessible as e.g.:
            nl = Read_namelist()    # with no name specified, it searches for a .ini file in the current dir
            itot = nl['grid']['itot']
            enttime = nl['time']['endtime']
    """

    def __init__(self,  namelist_file=None, ducktype=True):
        if (namelist_file is None):
            namelist_file = _find_namelist_file()

        self.groups = {}   # Dictionary holding all the data
        with open(namelist_file) as f:
            for line in f:
                lstrip = line.strip()
                if (len(lstrip) > 0 and lstrip[0] != "#"):
                    if lstrip[0] == '[' and lstrip[-1] == ']':
                        curr_group_name = lstrip[1:-1]
                        self.groups[curr_group_name] = {}
                    elif ("=" in line):
                        var_name = lstrip.split('=')[0]
                        value    = lstrip.split('=')[1]

                        if ducktype:
                            value = _convert_value(value)

                        self.groups[curr_group_name][var_name] = value

    def __getitem__(self, name):
        """
        Get group dictionary with `nl['group_name']` syntax
        """
        if name in self.groups.keys():
            return self.groups[name]
        else:
            raise RuntimeError(
                'Can\'t find group \"{}\" in .ini file'.format(name))


    def __repr__(self):
        """
        Print list of availabe groups
        """
        return 'Available groups:\n{}'.format(', '.join(self.groups.keys()))


    def set_value(self, group, variable, value):
        """
        Set value in namelist file/dict, if the group or
        variable does not exist, it is newly defined
        """
        if group not in self.groups:
            self.groups[group] = {}
        self.groups[group][variable] = value


    def save(self, namelist_file, allow_overwrite=False):
        """
        Write namelist from (nested) dictionary back to .ini file
        """
        if os.path.exists(namelist_file) and not allow_overwrite:
            raise Exception('.ini file \"{}\" already exists!'.format(namelist_file))

        with open(namelist_file, 'w') as f:
            for group in self.groups:
                f.write('[{}]\n'.format(group))
                for variable, value in self.groups[group].items():
                    if isinstance(value, list):
                        value = ','.join(str(v) for v in value)
                    elif isinstance(value, bool):
                        value = '1' if value else '0'
                    f.write('{}={}\n'.format(variable, value))
                f.write('\n')


def replace_namelist_value(item, new_value, group=None, namelist_file=None):
    """ Replace a item value in an existing namelist """
    if namelist_file is None:
        namelist_file = _find_namelist_file()
    with open(namelist_file, "r") as source:
        lines = source.readlines()
    with open(namelist_file, "w") as source:
        current_group = None
        has_replaced = False

        for line in lines:
            lstrip = line.strip()

            if len(lstrip)>0 and lstrip[0] == '[' and lstrip[-1] == ']':
                current_group = lstrip[1:-1]

            if group is None or group==current_group:
                source.write(re.sub(r'({}).*'.format(item), r'\1={}'.format(new_value), line))
                has_replaced = True
            else:
                source.write(line)

        if (not has_replaced):
            raise RuntimeError(
                'There is no item \"{0}\" in group \"{1}\" in .ini file'.format(item, group))


def determine_ntasks():
    namelist = Read_namelist()['master']

    npx = namelist['npx'] if 'npx' in namelist.keys() else 1
    npy = namelist['npy'] if 'npy' in namelist.keys() else 1

    return npx * npy


class Read_statistics:
    """ Read all the NetCDF statistics
        Example:
        f = Read_statistics('drycblles.default.0000000.nc')
        print(f) prints a list with the available variables
        The data can be accessed as either f['th'] or f.th, which returns the numpy array with data
        The variable names can be accessed as f.names['th'], the units as f.units['th'], the dimensions as f.dimensions['th']
        This allows you to automatically format axis labels as e.g.:
        pl.xlabel("{0:} ({1:})".format(f.names['th'], f.units['th']))
        """

    def __init__(self, stat_file):
        f = nc.Dataset(stat_file, 'r')

        # Dictionaries which hold the variable names, units, etc.
        self.data = {}
        self.units = {}
        self.names = {}
        self.dimensions = {}

        # For each variable in the NetCDF file, read all the content and info
        for var in f.variables:
            self.data[var] = f.variables[var].__array__()
            self.units[var] = f.variables[var].units
            self.names[var] = f.variables[var].long_name
            self.dimensions[var] = f.variables[var].dimensions

        f.close()

    def __getitem__(self, name):
        if name in self.data.keys():
            return self.data[name]
        else:
            raise RuntimeError(
                'Can\'t find variable \"{}\" in statistics file'.format(name))

    def __getattr__(self, name):
        if name in self.data.keys():
            return self.data[name]
        else:
            raise RuntimeError(
                'Can\'t find variable \"{}\" in statistics file'.format(name))

    def __repr__(self):
        return 'Available variables:\n{}'.format(', '.join(self.names.keys()))


class Read_grid:
    """ Read the grid file from MicroHH.
        If no file name is provided, grid.0000000 from the current directory is read """

    def __init__(self, itot, jtot, ktot, filename=None):
        self.en = '<' if sys.byteorder == 'little' else '>'
        filename = 'grid.0000000' if filename is None else filename
        self.TF = round(os.path.getsize(filename) /
                        (2 * itot + 2 * jtot + 2 * ktot))
        if self.TF == 8:
            self.prec = 'd'
        else:
            self.prec = 'f'

        self.fin = open(filename, 'rb')

        self.dim = {}

        self.dim['zh'] = np.zeros(ktot+1)

        self.dim['x'] = self.read(itot)
        self.dim['xh'] = self.read(itot)
        self.dim['y'] = self.read(jtot)
        self.dim['yh'] = self.read(jtot)
        self.dim['z'] = self.read(ktot)
        self.dim['zh'][:-1] = self.read(ktot)

        self.dim['zh'][-1] = self.dim['z'][-1] + 2*(self.dim['z'][-1] - self.dim['zh'][-2])

        self.fin.close()
        del self.fin

    def read(self, n):
        return np.array(
            st.unpack(
                '{0}{1}{2}'.format(
                    self.en, n, self.prec), self.fin.read(
                    n * self.TF)))


class Read_binary:
    """ Read a binary file from MicroHH. """

    def __init__(self, grid, filename):
        self.en = grid.en
        self.prec = grid.prec
        self.TF = grid.TF

        try:
            self.file = open(filename, 'rb')
        except BaseException:
            raise Exception('Cannot find file {}'.format(filename))

    def close(self):
        self.file.close()

    def read(self, n):
        return np.array(
            st.unpack(
                '{0}{1}{2}'.format(
                    self.en, n, self.prec), self.file.read(
                    n * self.TF)))


class Create_ncfile():
    def __init__(
            self,
            grid,
            filename,
            varname,
            dimensions,
            precision='',
            compression=True):
        self.ncfile = nc.Dataset(filename, "w", clobber=True)
        if not precision:
            precision = 'f{}'.format(grid.TF)
        elif precision == 'single':
            precision = 'f4'
        else:
            precision = 'f8'

        half_level_vars = [
            'w',
            'sw_flux_dn', 'sw_flux_dn_dir', 'sw_flux_up',
            'sw_flux_dn_clear', 'sw_flux_dn_dir_clear', 'sw_flux_up_clear',
            'lw_flux_dn', 'lw_flux_up'
            'lw_flux_dn_clear', 'lw_flux_up_clear']

        if(varname == 'u'):
            try:
                dimensions['xh'] = dimensions.pop('x')
            except KeyError:
                pass
        if(varname == 'v'):
            try:
                dimensions['yh'] = dimensions.pop('y')
            except KeyError:
                pass
        if(varname in half_level_vars):
            try:
                dimensions['zh'] = dimensions.pop('z')
            except KeyError:
                pass

        # create dimensions in netCDF file
        self.dim = {}
        self.dimvar = {}
        for key, value in dimensions.items():
            self.dim[key] = self.ncfile.createDimension(key, len(value))
            self.dimvar[key] = self.ncfile.createVariable(
                key, precision, (key))
            if key != 'time':
                self.dimvar[key][:] = grid.dim[key][value]

        self.var = self.ncfile.createVariable(
            varname, precision, tuple(
                self.sortdims(
                    dimensions.keys())), zlib=compression)

    def sync(self):
        self.ncfile.sync()

    def close(self):
        self.ncfile.close()

    def sortdims(self, lst=[]):
        ordered_dims = ['time', 'z', 'zh', 'y', 'yh', 'x', 'xh']
        lst_out = [value for value in ordered_dims if value in lst]
        return lst_out


def get_cross_indices(variable, mode):
    """ Find the cross-section indices given a variable name and mode (in 'xy','xz','yz') """
    if mode not in ['xy', 'xz', 'yz']:
        raise ValueError('\"mode\" should be in {\"xy\", \"xz\", \"yz\"}')

    # Get a list of all the cross-section files
    files = glob.glob('{}.{}.*.*'.format(variable, mode))
    if len(files) == 0:
        raise Exception('Cannot find any cross-section')

    # Get a list with all the cross-section files for one time
    time = files[0].split('.')[-1]
    files = glob.glob('{}.{}.*.{}'.format(variable, mode, time))

    # Get the indices
    indices = sorted([int(f.split('.')[-2]) for f in files])
    return indices


_opts = {
    'blue': '\033[94m',
    'green': '\033[92m',
    'purple': '\033[95m',
    'red': '\033[91m',
    'yellow': '\033[93m',
    'bf': '\033[1m',
    'ul': '\033[4m',
    'end': '\033[0m'
}


def print_header(message, time=True):
    """
    Format of print statements indicating new main routine
    """
    if time:
        now = datetime.datetime.now()
        print(
            '[{}] {}{}{}'.format(
                now.strftime('%d-%m: %H:%M'),
                _opts['green'],
                message,
                _opts['end']))
    else:
        print(
            '{}{}{}{}'.format(
                _opts['green'],
                _opts['bf'],
                message,
                _opts['end']))


def print_message(message):
    """
    Format of print statements
    """
    print(' - {}'.format(message))


def print_warning(message):
    """
    Format of print warnings
    """
    print(
        '{}{}WARNING:{} {}'.format(
            _opts['yellow'],
            _opts['bf'],
            _opts['end'],
            message))


def print_error(message):
    """
    Format of print errors
    """
    print(
        '{}{}ERROR:{} {}'.format(
            _opts['red'],
            _opts['bf'],
            _opts['end'],
            message))


def merge_options(options, options_to_add):
    """
    Merge dictionaries of dicts with run options.
    """
    for group in options_to_add:
        if group in options:
            options[group].update(options_to_add[group])
        else:
            options[group] = copy.deepcopy(options_to_add[group])


def run_scripts(scripts):
    def exec_function(lib, function, *args):
        rc = getattr(lib, function)(*args)

        if (rc is not None) and (rc != 0):
            raise Exception(
                '{}: {}() returned {}'.format(
                    script, function, rc))

    if scripts is not None:
        # Loop over, and execute all functions
        for script, functions in scripts.items():
            if (script == __file__):
                lib = sys.modules[__name__]
            else:
                # Module name = script name minus the `.py`
                module = script.replace('.py', '')
                # The full module name is relative to the source file, with
                # dots instead of slashes
                full_module = os.path.relpath(
                    os.getcwd(), sys.path[0]).replace(
                    '/', '.') + '.' + module
                # Import module; this executes all code that is not in
                # classes/functions
                if full_module not in sys.modules:
                    lib = importlib.import_module(full_module)
                else:
                    importlib.reload(sys.modules[full_module])

            # If any specific routines are specified, run them
            if functions is not None:
                for function in functions:
                    args = function[1:]
                    exec_function(lib, function[0], *args)


def restart_pre(origin, timestr):
    fnames = glob.glob('../' + origin + '/*_input.nc')
    fnames += glob.glob('../' + origin + '/grid.0000000')
    fnames += glob.glob('../' + origin + '/fftwplan.0000000')
    fnames += glob.glob('../' + origin + '/*.' + timestr)
    for file in fnames:
        shutil.copy(file, '.')


def compare_bitwise(f1, f2):
    # Compare with Python's `filecmp`
    cmp_python = filecmp.cmp(f1, f2)

    # Backup check with OS `cmp`
    sp = subprocess.Popen(
            'cmp {} {}'.format(f1, f2),
            executable='/bin/bash',
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    out, err = sp.communicate()
    sp.wait()
    cmp_os = not sp.returncode

    return cmp_python, cmp_os


def restart_post(origin, timestr):
    file_names = glob.glob('*.' + timestr)
    not_identical = False
    for file_name in file_names:
        cmp_python, cmp_os = compare_bitwise('../' + origin + '/' + file_name, file_name)

        if not cmp_python and not cmp_os:
            not_identical = True
            print_warning('{} is not identical (python+OS)'.format(file_name))
        elif not cmp_python:
            not_identical = True
            print_warning('{} is not identical (python)'.format(file_name))
        elif not cmp_os:
            not_identical = True
            print_warning('{} is not identical (OS)'.format(file_name))

    if not_identical:
        raise Warning('One or more restart files are not identical.')


def compare(origin, file, starttime=-1, vars={}):
    nc_new = nc.Dataset(file, mode="r")
    nc_old = nc.Dataset('../' + origin + '/' + file, mode="r")

    blacklist = ['iter']
    rtol = 1e-3
    atol = 1e-8
    if len(vars) == 0:
        for key in nc_new.variables.keys():
            if key not in blacklist:
                vars.update({key: [rtol, atol]})

    for key, opts in vars.items():
        var_new = np.mean(nc_new.variables[key][starttime:, ...], axis=0)
        var_old = np.mean(nc_old.variables[key][starttime:, ...], axis=0)
        if not np.allclose(
                var_new,
                var_old,
                rtol=opts[0],
                atol=opts[1],
                equal_nan=True):
            with np.errstate(all='ignore'):
                raise Warning('{0} in {1} has a relative error of up to {2:.2%}'.format(
                    key, file, np.max(np.abs((var_new - var_old) / var_old))))


def execute(command):
    sp = subprocess.Popen(
        command,
        executable='/bin/bash',
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    out, err = sp.communicate()
    sp.wait()

    # Write the standard output and errors to logflies
    with open('stdout.log', 'a') as f:
        f.write(out.decode('utf-8'))
    with open('stderr.log', 'a') as f:
        f.write(err.decode('utf-8'))

    if sp.returncode != 0:
        raise Exception(
            '\'{}\' returned \'{}\'.'.format(
                command, sp.returncode))

    return sp.returncode


def run_cases(cases, executable, mode, outputfile=''):
    """
    Function that iterates over a list of cases and runs all of them
    """

    if not os.path.exists(executable):
        raise Exception(
            'ERROR: Executable {} does not exists'.format(executable))

    # Get the absolute path to the executable
    executable_rel = executable
    executable = os.path.abspath(executable)
    rootdir = os.getcwd()

    for case in cases:
        print_header(
            'Running case \'{}\' for executable \'{}\' in dir \'{}\''.format(
                case.name, executable_rel, case.rundir))

        # Move to working directory
        rootdir = os.getcwd()
        rundir = rootdir + '/' + case.casedir + '/' + case.rundir + '/'

        casedir = rootdir + '/' + case.casedir + '/'
        if case.rundir != '':
            try:
                shutil.rmtree(rundir)
            except Exception:
                pass
            os.mkdir(rundir)
            os.chdir(rundir)

            try:
                for fname in case.files:
                    shutil.copy(casedir + fname, rundir)
            except BaseException:
                print_warning(
                    case.name +
                    ': Cannot find {} for copying, skipping case!'.format(
                        casedir +
                        fname))
                os.chdir(rootdir)
                continue
        else:
            case.keep = True

        try:
            # Update .ini file for testing
            ini_file = '{0}.ini'.format(case.name)
            nl = Read_namelist(ini_file, ducktype=False)

            for group, group_dict in case.options.items():
                for variable, value in group_dict.items():
                    nl.set_value(group, variable, value)

            nl.save(ini_file, allow_overwrite=True)

            # Find the number of MPI tasks
            ntasks = determine_ntasks()

            # Create input data, and do other pre-processing
            run_scripts(case.pre)

            for phase in case.phases:
                case.time = timeit.default_timer()
                if mode == 'cpu' or mode == 'gpu':
                    execute('{} {} {}'.format(executable, phase, case.name))
                elif mode == 'cpumpi':
                    execute('mpiexec --oversubscribe -n {} {} {} {}'.format(
                        ntasks, executable, phase, case.name))
                else:
                    raise ValueError('{} is an illegal value for mode'.format(mode))

                case.time = timeit.default_timer() - case.time

            # Run the post-processing steps
            run_scripts(case.post)
            case.success = True

        except Exception as e:
            print(str(e))
            print_error('Case Failed!')
            case.success = False
        else:
            print_message('Success!')

        finally:
            # Go back to root of all cases
            os.chdir(rootdir)

    # Write the output file and remove unnecssary dirs
    if outputfile != '':
        with open(outputfile, 'w') as csv_file:
            write = csv.writer(csv_file)
            write.writerow(['Name', 'Run Dir', 'Success', 'Time', 'Options'])
            for case in cases:
                write.writerow(
                    [case.name, case.rundir, case.success, case.time, case.options])
        csv_file.close()

    for case in cases:
        if case.success and not case.keep:
            rundir = rootdir + '/' + case.name + '/' + case.rundir + '/'
            shutil.rmtree(rundir)

"""
def generator_restart(cases):
    cases_out = []
    for case in cases:
        nl = Read_namelist('{0}/{0}.ini'.format(case.name))
        # Everything relevant is in the time group, so merge that with the
        # overriding options
        options = {'iotimeprec': 0}
        options.update(nl['time'])
        if case.options is not None:
            options.update(case.options)

        iotimeprec = options['iotimeprec']
        endtime = options['endtime']
        savetime = int(endtime / 2)
        endtimestr = '{0:07d}'.format(endtime * 10**(-iotimeprec))
        savetimestr = '{0:07d}'.format(savetime * 10**(-iotimeprec))

        case_init = case
        case_init.rundir = 'init'
        case_init.options.update({'savetime': savetime, 'endtime': endtime})

        case_restart = copy.deepcopy(case)
        case_restart.rundir = 'restart'
        case_restart.phases = ['run']
        case_restart.options.update(
            {'starttime': savetime, 'endtime': endtime})
        case_restart.pre = {__file__: [
            ['restart_pre', case_init.rundir, savetimestr]]}
        case_restart.post = {__file__: [
            ['restart_post', case_init.rundir, endtimestr]]}

        cases_out.append(case_init)
        cases_out.append(case_restart)

    return cases_out
"""

def generator_restart(case, endtime):
    cases_out = []
    nl = Read_namelist('{}/{}.ini'.format(case.casedir, case.name))

    iotimeprec = nl['time']['iotimeprec'] if 'iotimeprec' in nl['time'] else 0
    savetime = endtime/2

    savetime_io = int(round(savetime * 10**(-iotimeprec)))
    endtime_io = int(round(endtime * 10**(-iotimeprec)))

    endtimestr = '{0:07d}'.format(endtime_io)
    savetimestr = '{0:07d}'.format(savetime_io)

    case_init = copy.deepcopy(case)
    case_init.rundir = case.rundir + '_init'

    merge_options(case_init.options, {'time': {'savetime': savetime, 'endtime': endtime}})

    case_restart = copy.deepcopy(case)
    case_restart.rundir = case.rundir + '_restart'
    case_restart.phases = ['run']
    case_restart.pre = {__file__: [
        ['restart_pre', case_init.rundir, savetimestr]]}
    case_restart.post = {__file__: [
        ['restart_post', case_init.rundir, endtimestr]]}

    merge_options(case_restart.options, {'time': {'starttime': savetime, 'savetime': savetime, 'endtime': endtime}})

    cases_out.append(case_init)
    cases_out.append(case_restart)

    return cases_out


def prime_factors(n):
    import math

    result = []
    for i in range(2, int(math.sqrt(n)) + 1):
        # while i divides n , print i ad divide n
        while n % i == 0:
            result.append(i),
            n = n / i

    if (n > 1):
        result.append(int(n))

    return result


def generator_scaling(cases, procs, type='strong', dir='y'):
    cases_out = []
    for case in cases:
        if type == 'weak':
            nl = Read_namelist('{0}/{0}.ini'.format(case.name))
            itot = nl['grid']['itot']
            jtot = nl['grid']['jtot']
            xsize = nl['grid']['xsize']
            ysize = nl['grid']['ysize']

        for proc in procs:
            if dir == 'x':
                option = {'npx': proc}
            elif dir == 'y':
                option = {'npy': proc}
            elif dir == 'xy':
                primes = prime_factors(proc)
                npy = 1
                npx = 1
                for i in range(0, len(primes), 2):
                    npy *= primes[i]
                    if i + 1 < len(primes):
                        npx *= primes[i + 1]
                option = {'npy': npy, 'npx': npx}
            if type == 'weak':
                option.update({'itot': itot * npx,
                               'jtot': jtot * npy,
                               'xsize': xsize * npx,
                               'ysize': ysize * npy})
            new_case = copy.deepcopy(case)
            new_case.options.update(option)
            new_case.rundir = '{0:03d}'.format(proc)
            cases_out.append(new_case)

    return cases_out


"""
def generator_parameter_change(cases, **kwargs):
    cases_out = []
    if len(kwargs) > 0:
        for case in cases:
            key, value = list(kwargs.items())[0]
            for val in value:
                new_case = copy.deepcopy(case)
                new_case.options.update({key : val})
                new_case.rundir += (key + str(val)).replace('.', '')

                cases_out.append(new_case)
        del kwargs[key]
        if len(kwargs) > 0:
            cases_out = generator_parameter_change(cases_out, **kwargs)

    return cases_out
"""

def generator_parameter_permutations(base_case, lists):
    """
    Function to permutate lists of dictionaries to generate cases to run
    """
    cases_out = []

    # Put a single dictionary into a list with one item.
    if type(lists) is dict:
        lists = [lists]

    # Convert the dictionaries into tuples to enable to permutate the list.
    tuple_lists = []
    for l in lists:
        tuple_list = []
        for name, name_dict in l.items():
            tuple_list.append((name, name_dict))
        tuple_lists.append(tuple_list)

    # Create permutation of all lists. Each item contains 1 value of each list.
    lists_permutations = list(itertools.product(*tuple_lists))

    for lp in lists_permutations:
        case = copy.deepcopy(base_case)

        # Construct the directory name from tuple names.
        for name_dict in lp:
            case.rundir += '_' + name_dict[0]

        for name_dict in lp:
            merge_options(case.options, name_dict[1])

        cases_out.append(case)

    return cases_out


class Case:
    """
    Class that contains a case to run with the required runtime settings
    """
    def __init__(
            self,
            name,
            options=[],
            pre={},
            post={},
            phases=['init', 'run'],
            casedir='',
            rundir='default_run',
            files=[],
            keep=True):

        self.name = name       # Case name
        self.options = options # List of options to override
        self.pre = pre         # List of pre-processing python scripts
        self.post = post       # List of post-processing python scripts
        self.phases = phases   # List of the run phases we have to go through
        self.casedir = casedir # Directory of the case; self.name by default
        self.rundir = rundir   # Relative run directory, defaults to `default_run`
        self.files = files     # List of files necessary to run the case
        self.success = None    # Whether the entire case was run succesfully or not
        self.time = None       # Duration of the last phase (usually run)
        self.keep = keep       # Whether to keep the results of succefull simulations afterwards

        # By default; run {name}_input.py in preprocessing phase
        self.pre = pre if pre else {'{}_input.py'.format(name): None}
        self.files = files if files else [
            '{0}.ini'.format(name), '{}_input.py'.format(name)]
        self.casedir = casedir if casedir else name


def run_case(
        case_name, options_in, options_mpi_in,
        executable='microhh', mode='cpu',
        case_dir='.', experiment='local',
        additional_pre_py={}):

    options = deepcopy(options_in)

    if mode == 'cpumpi':
        merge_options(options, options_mpi_in)

    if additional_pre_py:
        # Aarghh
        pre = {'{}_input.py'.format(case_name): None}

        files = [
            '{}_input.py'.format(case_name),
            '{}.ini'.format(case_name)]

        for key, value in additional_pre_py.items():
            pre[key] = value
            files.append(key)

        cases = [
            Case(
                case_name,
                casedir=case_dir,
                rundir=experiment,
                options=options,
                pre=pre,
                files=files)]
    else:
        cases = [
            Case(
                case_name,
                casedir=case_dir,
                rundir=experiment,
                options=options)]

    run_cases(
        cases,
        executable,
        mode,
        outputfile='{}/{}_{}.csv'.format(case_dir, case_name, experiment))

    for case in cases:
        if not case.success:
            return 1
    return 0


def run_permutations(
        case_name, options_in, options_mpi_in, permutations_in,
        executable='microhh', mode='cpu',
        case_dir='.', experiment='local'):

    options = deepcopy(options_in)

    if mode == 'cpumpi':
        merge_options(options, options_mpi_in)

    base_case = Case(
            case_name,
            casedir=case_dir,
            rundir=experiment,
            options=options)

    cases = generator_parameter_permutations(base_case, permutations_in)

    run_cases(
        cases,
        executable,
        mode,
        outputfile='{}/{}_{}.csv'.format(case_dir, case_name, experiment))

    for case in cases:
        if not case.success:
            return 1
    return 0


def run_restart(
        case_name, options_in, options_mpi_in, permutations_in=None,
        executable='microhh', mode='cpu',
        case_dir='.', experiment='local'):

    # Deep copy the small version of the reference case and disable stats.
    options = deepcopy(options_in)

    if mode == 'cpumpi':
        merge_options(options, options_mpi_in)

    if permutations_in is None:
        base_cases = [Case(
                case_name,
                casedir=case_dir,
                rundir=experiment,
                options=options)]
    else:
        base_case = Case(
            case_name,
            casedir=case_dir,
            rundir=experiment,
            options=options)

        base_cases = generator_parameter_permutations(base_case, permutations_in)

    cases = []
    for case in base_cases:
        cases.extend(generator_restart(case, options['time']['endtime']))

    run_cases(
        cases,
        executable,
        mode,
        outputfile='{}/{}_restart_{}.csv'.format(case_dir, case_name, experiment))

    for case in cases:
        if not case.success:
            return 1
    return 0


def copy_radfiles(srcdir=None, destdir=None, gpt='128_112'):
    if srcdir is None:
        srcdir = os.path.dirname(inspect.getabsfile(inspect.currentframe()))+'/../rte-rrtmgp-cpp/rte-rrtmgp/' 
    if destdir is None:
        destdir = os.getcwd()
    if gpt == '128_112':
        shutil.copy(srcdir+'rrtmgp/data/rrtmgp-data-lw-g128-210809.nc', destdir+'/coefficients_lw.nc')
        shutil.copy(srcdir+'rrtmgp/data/rrtmgp-data-sw-g112-210809.nc', destdir+'/coefficients_sw.nc')
    elif gpt == '256_224':
        shutil.copy(srcdir+'rrtmgp/data/rrtmgp-data-lw-g256-210809.nc', destdir+'/coefficients_lw.nc')
        shutil.copy(srcdir+'rrtmgp/data/rrtmgp-data-sw-g224-210809.nc', destdir+'/coefficients_sw.nc')
    else:
        raise ValueError('gpt should be in {\'128_112\', \'256_224\'}')

    shutil.copy(srcdir+'extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc', destdir+'/cloud_coefficients_lw.nc')
    shutil.copy(srcdir+'extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc', destdir+'/cloud_coefficients_sw.nc')

def copy_lsmfiles(srcdir=None, destdir=None):
    if srcdir is None:
        srcdir = os.path.dirname(inspect.getabsfile(inspect.currentframe()))+'/../misc/'
    if destdir is None:
        destdir = os.getcwd()
    shutil.copy(srcdir+'van_genuchten_parameters.nc', destdir)
    