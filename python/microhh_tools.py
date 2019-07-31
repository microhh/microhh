import netCDF4 as nc
import numpy   as np
import struct  as st
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
    except:
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
        raise RuntimeError('Can\'t find any .ini files in the current directory!')
    if len(namelist_file) > 1:
        raise RuntimeError('There are multiple .ini files: {}'.format(namelist_file))
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
            printing e.g. nl['grid'] provides an overview of the available variables in a group
    """
    def __init__(self, namelist_file=None):
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
                        value = _convert_value(lstrip.split('=')[1])
                        self.groups[curr_group_name][var_name] = value

    def __getitem__(self, name):
        if name in self.groups.keys():
            return self.groups[name]
        else:
            raise RuntimeError('Can\'t find group \"{}\" in .ini file'.format(name))

    def __repr__(self):
        return 'Available groups:\n{}'.format(', '.join(self.groups.keys()))


def replace_namelist_value(variable, new_value, namelist_file=None):
    """ Replace a variables value in an existing namelist """
    if namelist_file is None:
        namelist_file = _find_namelist_file()

    with open(namelist_file, "r") as source:
        lines = source.readlines()
    with open(namelist_file, "w") as source:
        for line in lines:
            source.write(re.sub(r'({}).*'.format(variable), r'\1={}'.format(new_value), line))

def determine_mode():
    namelist = Read_namelist()['master']

    npx = namelist['npx'] if 'npx' in namelist.keys() else 1
    npy = namelist['npy'] if 'npy' in namelist.keys() else 1
    mode = 'serial' if npx*npy == 1 else 'parallel'
    return mode, npx*npy


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
        self.data       = {}
        self.units      = {}
        self.names      = {}
        self.dimensions = {}

        # For each variable in the NetCDF file, read all the content and info
        for var in f.variables:
            self.data[var]       = f.variables[var].__array__()
            self.units[var]      = f.variables[var].units
            self.names[var]      = f.variables[var].long_name
            self.dimensions[var] = f.variables[var].dimensions

        f.close()

    def __getitem__(self, name):
        if name in self.data.keys():
            return self.data[name]
        else:
            raise RuntimeError('Can\'t find variable \"{}\" in statistics file'.format(name))

    def __getattr__(self, name):
        if name in self.data.keys():
            return self.data[name]
        else:
            raise RuntimeError('Can\'t find variable \"{}\" in statistics file'.format(name))

    def __repr__(self):
        return 'Available variables:\n{}'.format(', '.join(self.names.keys()))


class Read_grid:
    """ Read the grid file from MicroHH.
        If no file name is provided, grid.0000000 from the current directory is read """
    def __init__(self, itot, jtot, ktot, filename=None):
        self.en = '<' if sys.byteorder == 'little' else '>'
        filename = 'grid.0000000' if filename is None else filename
        self.TF  = round(os.path.getsize(filename)/(2*itot + 2*jtot + 2*ktot))
        if self.TF == 8:
            self.prec = 'd'
        else:
            self.prec = 'f'

        self.fin = open(filename, 'rb')

        self.dim = {}
        self.dim['x']  = self.read(itot)
        self.dim['xh'] = self.read(itot)
        self.dim['y']  = self.read(jtot)
        self.dim['yh'] = self.read(jtot)
        self.dim['z']  = self.read(ktot)
        self.dim['zh'] = self.read(ktot)

        self.fin.close()
        del self.fin

    def read(self, n):
        return np.array(st.unpack('{0}{1}{2}'.format(self.en, n, self.prec), self.fin.read(n*self.TF)))

class Read_binary:
     """ Read a binary file from MicroHH. """
     def __init__(self, grid, filename):
        self.en   = grid.en
        self.prec = grid.prec
        self.TF   = grid.TF

        try:
            self.file = open(filename, 'rb')
        except:
            raise Exception('Cannot find file {}'.format(filename))

     def close(self):
        self.file.close()

     def read(self, n):
        return np.array(st.unpack('{0}{1}{2}'.format(self.en, n, self.prec), self.file.read(n*self.TF)))

class Create_ncfile():
    def __init__(self, grid, filename, varname, dimensions, precision='',compression=True):
        self.ncfile = nc.Dataset(filename, "w", clobber=False)
        if not precision:
            precision = 'f{}'.format(grid.TF)
        elif precision == 'single':
            precision = 'f4'
        else:
            precision = 'f8'
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
        if(varname == 'w'):
            try:
                dimensions['zh'] = dimensions.pop('z')
            except KeyError:
                pass
        # create dimensions in netCDF file
        self.dim = {}
        self.dimvar = {}
        for key, value in dimensions.items():
            self.dim[key]       = self.ncfile.createDimension(key, len(value))
            self.dimvar[key]    = self.ncfile.createVariable(key, precision, (key))
            if key is not 'time':
                self.dimvar[key][:] = grid.dim[key][value]
        self.var = self.ncfile.createVariable(varname, precision, tuple(self.sortdims(dimensions.keys())), zlib=compression)

    def sync(self):
        self.ncfile.sync()

    def close(self):
        self.ncfile.close()
    
    def sortdims(self, lst = []):
        ordered_dims = ['time', 'z', 'zh', 'y', 'yh', 'x', 'xh']
        lst_out = [value for value in ordered_dims if value in lst]
        return lst_out

def get_cross_indices(variable, mode):
    """ Find the cross-section indices given a variable name and mode (in 'xy','xz','yz') """
    if mode not in ['xy','xz','yz']:
        raise ValueError('\"mode\" should be in {\"xy\", \"xz\", \"yz\"}')

    # Get a list of all the cross-section files
    files = glob.glob('{}.{}.*.*'.format(variable, mode))
    if len(files) == 0:
        raise Exception('Cannot find any cross-section')

    # Get a list with all the cross-section files for one time
    time  = files[0].split('.')[-1]
    files = glob.glob('{}.{}.*.{}'.format(variable, mode, time))

    # Get the indices
    indices = [int(f.split('.')[-2]) for f in files]
    indices.sort()
    return indices

_opts = {
   'blue'   : '\033[94m',
   'green'  : '\033[92m',
   'purple' : '\033[95m',
   'red'    : '\033[91m',
   'yellow' : '\033[93m',
   'bf'     : '\033[1m',
   'ul'     : '\033[4m',
   'end'    : '\033[0m'
}

def print_header(message, time=True):
    """
    Format of print statements indicating new main routine
    """
    if time:
        now = datetime.datetime.now()
        print('[{}] {}{}{}'.format(now.strftime('%d-%m: %H:%M'), _opts['green'], message, _opts['end']))
    else:
        print('{}{}{}{}'.format(_opts['green'], _opts['bf'], message, _opts['end']))

def print_message(message):
    """
    Format of print statements
    """
    print(' - {}'.format(message))

def print_warning(message):
    """
    Format of print warnings
    """
    print('{}{}WARNING:{} {}'.format(_opts['yellow'], _opts['bf'], _opts['end'], message))

def print_error(message):
    """
    Format of print errors
    """
    print('{}{}ERROR:{} {}'.format(_opts['red'], _opts['bf'], _opts['end'], message))



def run_scripts(scripts):
    def exec_function(lib, function, *args):
        rc = getattr(lib, function)(*args)

        if rc != 0:
            raise Exception('{}: {}() returned {}'.format(script, function, rc))


    if scripts is not None:
        # Loop over, and execute all functions
        for script, functions in scripts.items():
            if (script == __file__):
                lib = sys.modules[__name__]
            else:
                # Module name = script name minus the `.py`
                module = script.replace('.py', '')
                # The full module name is relative to the source file, with dots instead of slashes
                full_module = os.path.relpath(os.getcwd(),sys.path[0]).replace('/','.')+'.'+module
                # Import module; this executes all code that is not in classes/functions
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
    fnames = glob.glob('../'+origin+'/*_input.nc')
    fnames += glob.glob('../'+origin+'/grid.0000000')
    fnames += glob.glob('../'+origin+'/fftwplan.0000000')
    fnames += glob.glob('../'+origin+'/*.'+timestr)
    for file in fnames:
        shutil.copy(file, '.')


def restart_post(origin, timestr):
    #Write a real function that compares relevant files between dir1 and dir2

    fnames = glob.glob('*.'+timestr)
    for file in fnames:
        if not filecmp.cmp('../'+origin+'/'+file, file):
            raise Warning(file + ' is not identical')

def compare(origin, file, starttime=-1, vars={}):
    nc_new = nc.Dataset(file, mode="r")
    nc_old = nc.Dataset('../'+origin+'/'+file, mode="r")

    blacklist=['iter']
    rtol=1e-3
    atol=1e-8
    if len(vars) == 0:
        for key in nc_new.variables.keys():
            if key not in blacklist:
                vars.update({key: [rtol, atol]})

    for key, opts in vars.items():
        var_new = np.mean(nc_new.variables[key][starttime:,...],axis=0)
        var_old = np.mean(nc_old.variables[key][starttime:,...],axis=0)
        if not np.allclose(var_new, var_old, rtol=opts[0], atol=opts[1], equal_nan=True):
            with np.errstate(all='ignore'):
                raise Warning('{0} in {1} has a relative error of up to {2:.2%}'.format(key, file, np.max(np.abs((var_new - var_old)/var_old))))


def execute(command):
    sp = subprocess.Popen(command, executable='/bin/bash', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = sp.communicate()
    sp.wait()

    # Write the standard output and errors to logflies
    with open('stdout.log', 'w') as  f:
        f.write(out.decode('utf-8'))
    with open('stderr.log', 'w') as  f:
        f.write(err.decode('utf-8'))

    if sp.returncode != 0:
        raise Exception('\'{}\' returned \'{}\'.'.format(command, sp.returncode))

def test_cases(cases, executable, outputfile=''):

    if not os.path.exists(executable):
        raise Exception('ERROR: Executable {} does not exists'.format(executable))

    # Get the absolute path to the executable
    executable_rel = executable
    executable = os.path.abspath(executable)
    rootdir = os.getcwd()

    for case in cases:
        print_header('Testing case \'{}\' for executable \'{}\''.format(case.name, executable_rel))
        # Determine whether to run serial or parallel

        # Move to working directory
        rootdir = os.getcwd()
        rundir  = rootdir + '/' + case.casedir +  '/' + case.rundir + '/'
        casedir = rootdir + '/' + case.casedir +  '/'
        if case.rundir is not '':
            try:
                shutil.rmtree(rundir)
            except Exception:
                pass
            os.mkdir(rundir)
            os.chdir(rundir)

            try:
                for fname in case.files:
                    shutil.copy(casedir+fname, rundir)
            except:
                print_warning(case.name +  ': Cannot find {} for copying,  skipping case!'.format(casedir+fname))
                os.chdir(rootdir)
                continue
        else:
            case.keep = True

        try:
            # Update .ini file for testing
            for variable, value in case.options.items():
                replace_namelist_value(variable, value, '{0}.ini'.format(case.name))
            mode, ntasks = determine_mode()

            # Create input data, and do other pre-processing
            run_scripts(case.pre)

            for phase in case.phases:
                case.time = timeit.default_timer()
                if mode == 'serial':
                    execute('{} {} {}'.format(executable, phase, case.name))
                elif mode == 'parallel':
                    execute('mpirun -n {} {} {} {}'.format(ntasks, executable, phase, case.name))
                case.time = timeit.default_timer() - case.time

            # Run the post-processing steps
            run_scripts(case.post)
            case.success = True
        except Exception as e:
            print(str(e))
            print_warning('Case Failed!')
            case.success = False
        finally:
            # Go back to root of all cases
            os.chdir(rootdir)

    #Write the output file and remove unnecssary dirs
    if outputfile is not '':
        with open(outputfile, 'w') as csvFile:
            write = csv.writer(csvFile)
            write.writerow(['Name', 'Run Dir', 'Success', 'Time', 'Options'])
            for case in cases:
                write.writerow([case.name, case.rundir, case.success, case.time, case.options])
        csvFile.close()

    for case in cases:
        if case.success and not case.keep:
            rundir  = rootdir + '/' + case.name +  '/' + case.rundir + '/'
            shutil.rmtree(rundir)

def generator_restart(cases):
    cases_out = []
    for case in cases:
        nl = Read_namelist('{0}/{0}.ini'.format(case.name))
        #Everything relevant is in the time group, so merge that with the overriding options
        options = {'iotimeprec' : 0}
        options.update(nl['time'])
        if case.options is not None:
            options.update(case.options)

        iotimeprec  = options['iotimeprec']
        endtime     = options['endtime']
        savetime    = int(endtime/2)
        endtimestr  = '{0:07d}'.format(endtime*10**(-iotimeprec))
        savetimestr = '{0:07d}'.format(savetime*10**(-iotimeprec))

        case_init = case
        case_init.rundir = 'init'
        case_init.options.update({'savetime' : savetime, 'endtime': endtime })

        case_restart = copy.deepcopy(case)
        case_restart.rundir = 'restart'
        case_restart.phases = ['run']
        case_restart.options.update({'starttime' : savetime, 'endtime': endtime })
        case_restart.pre     = {__file__ : [['restart_pre',  case_init.rundir, savetimestr]]}
        case_restart.post    = {__file__ : [['restart_post', case_init.rundir, endtimestr]]}

        cases_out.append(case_init)
        cases_out.append(case_restart)

    return cases_out

def primeFactors(n):
    import math

    result = []
    for i in range(2,int(math.sqrt(n))+1):
        # while i divides n , print i ad divide n
        while n % i== 0:
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
            itot  = nl['grid']['itot']
            jtot  = nl['grid']['jtot']
            xsize = nl['grid']['xsize']
            ysize = nl['grid']['ysize']

        for proc in procs:
            if dir == 'x':
                option = {'npx' : proc}
            elif dir == 'y':
                option = {'npy' : proc}
            elif dir == 'xy':
                primes = primeFactors(proc)
                npy = 1
                npx = 1
                for i in range(0,len(primes),2):
                    npy *= primes[i]
                    if i+1 < len(primes):
                        npx *= primes[i+1]
                option = {'npy' : npy, 'npx' : npx}
            if type == 'weak':
                option.update({'itot' : itot*npx, 'jtot' : jtot*npy, 'xsize' : xsize*npx, 'ysize' : ysize * npy})
            new_case = copy.deepcopy(case)
            new_case.options.update(option)
            new_case.rundir = '{0:03d}'.format(proc)
            cases_out.append(new_case)

    return cases_out

def generator_parameter_change(cases, **kwargs):
    cases_out = []
    if len(kwargs) > 0:
        for case in cases:
            key, value = list(kwargs.items())[0]
            for val in value:
                new_case = copy.deepcopy(case)
                new_case.options.update({key : val})
                new_case.rundir += (key+str(val)).replace('.','')

                cases_out.append(new_case)
        del kwargs[key]
        if len(kwargs) > 0:
            cases_out = generator_parameter_change(cases_out, **kwargs)

    return cases_out

class Case:
    def __init__(self, name, options={}, pre={}, post={}, phases = ['init','run'], casedir = '', rundir='', files=[], keep=False):

        self.name     = name        # Case name
        self.options  = options     # Override existing namelist options
        self.pre      = pre         # List of pre-processing python scripts
        self.post     = post        # List of post-processing python scripts
        self.phases   = phases      # List of the run phases we have to go through
        self.casedir  = casedir     # Directory of the case; self.name by default
        self.rundir   = rundir      # Relative run directory
        self.files    = files       # List of files necessary to run the case
        self.success  = None        # Whether the entire case was run succesfully or not
        self.time     = None        # Duration of the last phase (usually run)
        self.keep     = keep        # Whether to keep the results of succefull simulations afterwards

        # By default; run {name}_input.py in preprocessing phase
        self.pre = pre if pre else {'{}_input.py'.format(name): None}
        print(self.name, files)
        self.files = files if files else ['{0}.ini'.format(name), '{}_input.py'.format(name)]
        print(self.files)
        self.casedir = casedir if casedir else name
       
