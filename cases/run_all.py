import subprocess
import importlib
import shutil
import re
import os
import glob

from messages import *

def replace_namelist_value(namelist_file, variable, new_value):
    with open(namelist_file, "r") as source:
        lines = source.readlines()
    with open(namelist_file, "w") as source:
        for line in lines:
            source.write(re.sub(r'({}).*'.format(variable), r'\1={}'.format(new_value), line))


def determine_mode(namelist):
    npx = namelist['npx'] if 'npx' in namelist.keys() else 1
    npy = namelist['npy'] if 'npy' in namelist.keys() else 1
    mode = 'serial' if npx*npy == 1 else 'parallel'
    return mode, npx*npy


def run_scripts(case_name, scripts):
    def exec_function(lib, function):
        rc = getattr(lib, function)()

        if rc != 0:
            print_error('{}: {}() returned {}'.format(script, function, rc))

    if scripts is not None:
        # Loop over, and execute all functions
        for script, functions in scripts.items():

            # Module name = script name minus the `.py`
            module = script.replace('.py', '')

            # Import module; this executes all code that is not in classes/functions
            lib = importlib.import_module('{0}.{1}'.format(case_name, module))

            # If any specific routines are specified, run them
            if functions is not None:
                if isinstance(functions, list):
                    for function in functions:
                        exec_function(lib, function)
                else:
                    exec_function(lib, functions)


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
        print_error('\"{}\" returned \"{}\".'.format(command, sp.returncode))


def test_cases(cases, executable):
    if not os.path.exists(executable):
        raise Exception('ERROR: Executable {} does not exists'.format(executable))

    # Get the absolute path to the executable
    executable_rel = executable
    executable = os.path.abspath(executable)

    for case in cases:

        # Determine whether to run serial or parallel
        mode, ntasks = determine_mode(case.options)

        print_header('Testing case \"{}\" for executable \"{}\" ({})'.format(case.name, executable_rel, mode))

        # Move to working directory
        os.chdir(case.name)

        # Script only works if there is a `[name].ini` and `[name]prof.py` present
        if (os.path.exists('{0}.ini'.format(case.name)) and os.path.exists('{0}prof.py'.format(case.name))):

            # Backup original .ini file
            shutil.copy('{0}.ini'.format(case.name), '{0}.ini.original'.format(case.name))

            # Update .ini file for testing
            for variable, value in case.options.items():
                replace_namelist_value('{0}.ini'.format(case.name), variable, value)

            # Create input data, and do other pre-processing
            run_scripts(case.name, case.pre)

            # Run init and run phases
            execute('rm -f *000*')      # BvS: do we want this?
            for phase in ['init', 'run']:
                if mode == 'serial':
                    execute('{} {} {}'.format(executable, phase, case.name))
                elif mode == 'parallel':
                    execute('mpirun -n {} {} {} {}'.format(ntasks, executable, phase, case.name))

            # Run the post-processing steps
            run_scripts(case.name, case.post)

            # Restore original .ini file and remove backup
            shutil.copy('{0}.ini.original'.format(case.name), '{0}.ini'.format(case.name))
            os.remove('{0}.ini.original'.format(case.name))

        else:
            print_warning('cannot find {0}.ini and/or {0}prof.py, skipping...!'.format(case.name))

        # Go back to root of all cases
        os.chdir('..')


class Case:
    def __init__(self, name, options, pre=None, post=None):

        self.name     = name        # Case / directory name
        self.options  = options     # Override existing namelist options
        self.pre      = pre         # List of pre-processing python scripts
        self.post     = post        # List of post-processing python scripts

        # By default; run {name}prof.py in preprocessing phase
        if self.pre is None:
            self.pre = {'{}prof.py'.format(name): None}


if __name__ == "__main__":

    if (False):
        # Serial
        cases = [
            Case("bomex",     { "itot" : 16, "jtot" : 16,  "ktot" : 32, "endtime": 3600 }),
            Case("drycblles", { "itot" : 16, "jtot" : 16, "ktot" : 32 }, post={'validate.py': ['test1', 'test2']})
                ]

        test_cases(cases, '../build/microhh')
        test_cases(cases, '../build/microhh_single')

        # Parallel
        cases = [
            Case("bomex",     { "npx": 1, "npy": 2, "itot" : 16, "jtot" : 16,  "ktot" : 32, "endtime": 3600 }),
            Case("drycblles", { "npx": 1, "npy": 2, "itot" : 16, "jtot" : 16, "ktot" : 32 }, post={'validate.py': ['test1', 'test2']})
                ]

        test_cases(cases, '../build_parallel/microhh')
        test_cases(cases, '../build_parallel/microhh_single')

    if (False):
        blacklist = ['prandtlslope']

        # Settings for all test cases:
        settings_serial   = { "itot" : 16, "jtot" : 8,  "ktot" : 16, "endtime": 10 }
        settings_parallel = { "itot" : 16, "jtot" : 8,  "ktot" : 16, "endtime": 10, "npx": 1, "npy": 2 }

        dirs  = glob.glob('*')
        cases_serial   = []
        cases_parallel = []
        for dir in dirs:
            if os.path.isdir(dir) and dir not in blacklist:
                cases_serial  .append(Case(dir, settings_serial  ))
                cases_parallel.append(Case(dir, settings_parallel))

        test_cases(cases_serial, '../build/microhh')
        test_cases(cases_serial, '../build/microhh_single')

        test_cases(cases_parallel, '../build_parallel/microhh')
        test_cases(cases_parallel, '../build_parallel/microhh_single')


    if (True):
        # DANGER: checkout all ini files
        dirs  = glob.glob('*')
        for dir in dirs:
            if os.path.isdir(dir):
                execute('git checkout {0}/{0}.ini'.format(dir))

