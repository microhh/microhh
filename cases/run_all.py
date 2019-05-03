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

from messages import *

def replace_namelist_value(namelist_file, variable, new_value):
    with open(namelist_file, 'r') as source:
        lines = source.readlines()
    with open(namelist_file, 'w') as source:
        for line in lines:
            source.write(re.sub(r'({}).*'.format(variable), r'\1={}'.format(new_value), line))


def determine_mode(namelist):
    npx = namelist['npx'] if 'npx' in namelist.keys() else 1
    npy = namelist['npy'] if 'npy' in namelist.keys() else 1
    mode = 'serial' if npx*npy == 1 else 'parallel'
    return mode, npx*npy


def run_scripts(case_name, scripts):
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

                # Import module; this executes all code that is not in classes/functions
                lib = importlib.import_module('{0}.{1}'.format(case_name, module))

            # If any specific routines are specified, run them
            if functions is not None:
                for function in functions:
                    args = function[1:]
                    nerror += exec_function(lib, function[0], *args)
    return nerror

def restart_pre(origin, time):
    fnames = glob.glob('../'+origin+'/*_input.nc')
    fnames += glob.glob('../'+origin+'/grid.0000000')
    fnames += glob.glob('../'+origin+'/fftwplan.0000000')
    fnames += glob.glob('../'+origin+'/*'+time)
    for file in fnames:
        shutil.copy(file, '.')
    return 0        
    
       
def restart_post(origin, time):
    #Write a real function that compares relevant files between dir1 and dir2

    fnames = glob.glob('*'+time)
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


def test_cases(cases, executable):
    if not os.path.exists(executable):
        raise Exception('ERROR: Executable {} does not exists'.format(executable))

    # Get the absolute path to the executable
    executable_rel = executable
    executable = os.path.abspath(executable)

    with open('test_results.csv', 'w') as csvFile:
        write = csv.writer(csvFile)
        write.writerow(["Name", "Run Dir", "Success", "Time", "Options"])
        for case in cases:
            # Determine whether to run serial or parallel
            mode, ntasks = determine_mode(case.options)

            print_header('Testing case \'{}\' for executable \'{}\' ({})'.format(case.name, executable_rel, mode))
            nerror = 0
            # Move to working directory
            rootdir = os.getcwd()
            os.chdir(case.name)
            casedir = os.getcwd()
            rundir = casedir +  '/' + case.rundir + '/'
            try:
                shutil.rmtree(rundir)
            except Exception:
                pass
                
            os.mkdir(rundir)    
            
            try:
                for fname in case.files:
                    shutil.copy(fname, rundir)
            except:
                print_warning('Cannot find required files, skipping!')

            os.chdir(rundir)
            
            # Update .ini file for testing
            for variable, value in case.options.items():
                replace_namelist_value('{0}.ini'.format(case.name), variable, value)

            # Create input data, and do other pre-processing
            nerror += run_scripts(case.name, case.pre)

            for phase in case.phases:
                case.time = timeit.default_timer()
                if mode == 'serial':
                    nerror += execute('{} {} {}'.format(executable, phase, case.name))
                elif mode == 'parallel':
                    nerror += execute('mpirun -n {} {} {} {}'.format(ntasks, executable, phase, case.name))
                case.time = timeit.default_timer() - case.time
                
            # Run the post-processing steps
            nerror += run_scripts(case.name, case.post)

            case.success = (nerror == 0)
            
            write.writerow([case.name, case.rundir, case.success, case.time, case.options])
            # Go back to root of all cases
            os.chdir(rootdir)
            
            
    csvFile.close()

class Case:
    def __init__(self, name, options, pre=None, post=None, phases = ['init','run'], rundir='.', files=None):

        self.name     = name        # Case / directory name
        self.options  = options     # Override existing namelist options
        self.pre      = pre         # List of pre-processing python scripts
        self.post     = post        # List of post-processing python scripts
        self.phases   = phases      # List of the run phases we have to go through
        self.rundir   = rundir      # Relative run directory
        self.files    = files       # List of files necessary to run the case
        self.success  = None        # Whether the entire case was run succesfully or not
        self.time     = None        # Duration of the last phase (usually run)
        
        # By default; run {name}_input.py in preprocessing phase
        if self.pre is None:
            self.pre = {'{}_input.py'.format(name): None}
        if self.files is None:
            self.files = ['{}_input.py'.format(name), '{0}.ini'.format(name)]


if __name__ == '__main__':

    if (True):
        # Serial
        cases = [
            Case('bomex',     {'savetime'  : 1800 , 'endtime': 3600 }, rundir='run1'),
            Case('bomex',     {'starttime' : 1800, 'endtime': 3600 }, pre={__file__ : [['restart_pre', 'run1', '1800']]}, post={__file__ : [['restart_post', 'run1', '3600']]}, phases=['run'], rundir='run2'),

                ]

        test_cases(cases, '../build/microhh')

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

