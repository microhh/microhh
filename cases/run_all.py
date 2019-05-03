import subprocess
import importlib
import shutil
import re
import os
import glob
import sys
import filecmp
import timeit
import csv
import copy
sys.path.append('../python/')
import microhh_tools as mht

from messages import *


def determine_mode(namelist):
    npx = namelist['npx'] if 'npx' in namelist.keys() else 1
    npy = namelist['npy'] if 'npy' in namelist.keys() else 1
    mode = 'serial' if npx*npy == 1 else 'parallel'
    return mode, npx*npy


def run_scripts(case_name, run_dir, scripts):
    def exec_function(lib, function, *args):
        rc = getattr(lib, function)(*args)

        if rc != 0:
            print_error('{}: {}() returned {}'.format(script, function, rc))
            return 1
        return 0

    nerror = 0
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
                lib = importlib.import_module(full_module)
            # If any specific routines are specified, run them
            if functions is not None:
                for function in functions:
                    args = function[1:]
                    nerror += exec_function(lib, function[0], *args)

    return nerror

def restart_pre(origin, timestr):
    fnames = glob.glob('../'+origin+'/*_input.nc')
    fnames += glob.glob('../'+origin+'/grid.0000000')
    fnames += glob.glob('../'+origin+'/fftwplan.0000000')
    fnames += glob.glob('../'+origin+'/*.'+timestr)
    for file in fnames:
        shutil.copy(file, '.')
    return 0


def restart_post(origin, timestr):
    #Write a real function that compares relevant files between dir1 and dir2

    fnames = glob.glob('*.'+timestr)
    for file in fnames:
        if not filecmp.cmp('../'+origin+'/'+file, file):
            return file + ' is not identical'

    return 0

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
        print_error('\'{}\' returned \'{}\'.'.format(command, sp.returncode))
        return 1

    return 0


def test_cases(cases, executable, outputfile='test_results.csv'):
    if not os.path.exists(executable):
        raise Exception('ERROR: Executable {} does not exists'.format(executable))

    # Get the absolute path to the executable
    executable_rel = executable
    executable = os.path.abspath(executable)
    rootdir = os.getcwd()

    for case in cases:
        # Determine whether to run serial or parallel
        mode, ntasks = determine_mode(case.options)

        print_header('Testing case \'{}\' for executable \'{}\' ({})'.format(case.name, executable_rel, mode))
        nerror = 0
        # Move to working directory
        rootdir = os.getcwd()
        rundir  = rootdir + '/' + case.name +  '/' + case.rundir + '/'
        casedir = rootdir + '/' + case.name +  '/'

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
            print_warning('Cannot find required files, skipping!')

        # Update .ini file for testing
        for variable, value in case.options.items():
            mht.replace_namelist_value(variable, value, '{0}.ini'.format(case.name))

        # Create input data, and do other pre-processing
        nerror += run_scripts(case.name, case.rundir, case.pre)
        for phase in case.phases:
            case.time = timeit.default_timer()
            if mode == 'serial':
                nerror += execute('{} {} {}'.format(executable, phase, case.name))
            elif mode == 'parallel':
                nerror += execute('mpirun -n {} {} {} {}'.format(ntasks, executable, phase, case.name))
            case.time = timeit.default_timer() - case.time

        # Run the post-processing steps
        nerror += run_scripts(case.name, case.rundir, case.post)

        case.success = (nerror == 0)

        # Go back to root of all cases
        os.chdir(rootdir)

    #Write the output file and remove unnecssary dirs
    with open(outputfile, 'w') as csvFile:
        write = csv.writer(csvFile)
        write.writerow(['Name', 'Run Dir', 'Success', 'Time', 'Options'])
        for case in cases:
            write.writerow([case.name, case.rundir, case.success, case.time, case.options])
    csvFile.close()

    for case in cases:
        rundir  = rootdir + '/' + case.name +  '/' + case.rundir + '/'
        if case.success and not case.keep:
            shutil.rmtree(rundir)

def generator_restart(cases):
    cases_out = []
    for case in cases:
        nl = mht.Read_namelist('{0}/{0}.ini'.format(case.name))
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

        cases_out.append([case_init, case_restart])

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

def generator_strongscaling(cases, procs, dir='y'):
    cases_out = []
    for case in cases:
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
            new_case = copy.deepcopy(case)
            new_case.options.update(option)
            new_case.rundir = '{0:03d}'.format(proc)
            cases_out.append(new_case)

    return cases_out

class Case:
    def __init__(self, name, options={}, pre={}, post={}, phases = ['init','run'], rundir='', files=[], keep=False):

        self.name     = name        # Case / directory name
        self.options  = options     # Override existing namelist options
        self.pre      = pre         # List of pre-processing python scripts
        self.post     = post        # List of post-processing python scripts
        self.phases   = phases      # List of the run phases we have to go through
        self.rundir   = rundir      # Relative run directory
        self.files    = files       # List of files necessary to run the case
        self.success  = None        # Whether the entire case was run succesfully or not
        self.time     = None        # Duration of the last phase (usually run)
        self.keep     = keep        # Whether to keep the results of succefull simulations afterwards

        # By default; run {name}_input.py in preprocessing phase
        self.pre.update({'{}_input.py'.format(name): None})
        if self.files == []:
            self.files = ['{0}.ini'.format(name)]
        self.files += list(self.pre.keys()) + list(self.post.keys())


if __name__ == '__main__':

    if (True):
        # Serial
        #cases = generator_restart([Case('bomex')])
        cases = generator_strongscaling([Case('bomex')],procs=[1,2,4,8], dir='xy')
        print(cases)
        test_cases(cases, '../build_mpi/microhh')

    if (False):
        # Serial
        cases = [
            Case('bomex',     { 'itot' : 16, 'jtot' : 16,  'ktot' : 32, 'endtime': 3600 }),
            Case('drycblles', { 'itot' : 16, 'jtot' : 16, 'ktot' : 32 }, post={'validate.py': ['test1', 'test2']})
                ]

        test_cases(cases, '../build/microhh')


        # Parallel
        cases = [
            Case('bomex',     { 'npx': 1, 'npy': 2, 'itot' : 16, 'jtot' : 16,  'ktot' : 32, 'endtime': 3600 }),
            Case('drycblles', { 'npx': 1, 'npy': 2, 'itot' : 16, 'jtot' : 16, 'ktot' : 32 }, post={'validate.py': ['test1', 'test2']})
                ]

        test_cases(cases, '../build_parallel/microhh')


    if (False):
        blacklist = ['prandtlslope','moser600']

        # Settings for all test cases:
        settings_serial   = { 'itot' : 16, 'jtot' : 8,  'ktot' : 16, 'endtime': 10 }
        settings_parallel = { 'itot' : 16, 'jtot' : 8,  'ktot' : 16, 'endtime': 10, 'npx': 1, 'npy': 2 }

        dirs  = glob.glob('*')
        cases_serial   = []
        cases_parallel = []
        for dir in dirs:
            if os.path.isdir(dir) and dir not in blacklist:
                print(dir)
                cases_serial  .append(Case(dir, settings_serial  ))
                cases_parallel.append(Case(dir, settings_parallel))

        test_cases(cases_serial, '../build/microhh')

        #test_cases(cases_parallel, '../build_parallel/microhh')


#    if (False):
#        # DANGER: checkout all ini files
#        dirs  = glob.glob('*')
#        for dir in dirs:
#            if os.path.isdir(dir):
#                execute('git checkout {0}/{0}.ini'.format(dir))
