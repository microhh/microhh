import subprocess
import os
import sys
import glob
import re
import filecmp

from drycbllesprof import create_profiles

def execute(task, silent):
    if (silent):
        subprocess.call('{} > /dev/null'.format(task), shell=True, executable='/bin/bash')
    else:
        subprocess.call('{}'.format(task), shell=True, executable='/bin/bash')

def execute_with_return(task):
    sp = subprocess.Popen(task, shell=True, executable='/bin/bash', stdout=subprocess.PIPE)
    return sp.stdout.read().decode("utf-8")

def clean_workdir():
    for search in ['*.000*', '*.out']:
        files = glob.glob(search)
        for f in files:
            os.remove(f)

def replace_namelist_value(variable, new_value, namelist_file):
    """ Replace a variables value in an existing namelist """
    with open(namelist_file, "r") as source:
        lines = source.readlines()
    with open(namelist_file, "w") as source:
        for line in lines:
            source.write(re.sub(r'({}).*'.format(variable), r'\1={}'.format(new_value), line))

def update_namelists(mode, npx, npy):
    for namelist in ['drycblles.ini', 'drycblles_restart.ini']:
        replace_namelist_value('npx', npx, namelist)
        replace_namelist_value('npy', npy, namelist)

def link_executables(path):
    if path[-1] != '/':
        path += '/'

    if os.path.exists('microhh'):
        os.remove('microhh')
    if os.path.exists('microhh_single'):
        os.remove('microhh_single')

    os.symlink('{}microhh'.format(path),        'microhh')
    os.symlink('{}microhh_single'.format(path), 'microhh_single')

def files_are_identical(file1, file2):
    #return filecmp.cmp(file1, file2, False)    # Not the same results as command line `cmp`

    code = execute_with_return('cmp {} {}'.format(file1, file2))
    return True if code == "" else False


if __name__ == '__main__':

    serial_build   = '../../build/'
    parallel_build = '../../build_parallel/'

    # Number of tasks for parallel test
    npx = 2
    npy = 2

    # Mute MicroHH output
    silent = True

    # Keep track of whether all tasks are a success
    great_success = True

    # Test both the serial and parallel version
    # The serial version also runs with mpirun when ntasks=1
    for mode in ['serial', 'parallel']:
        print('Running {} tests'.format(mode))

        # Update the number of tasks in the namelists

        # Link the correct executables and update namelists
        if mode == 'serial':
            link_executables(serial_build)
            update_namelists(mode, 1, 1)
            ntasks = 1
        else:
            link_executables(parallel_build)
            update_namelists(mode, npx, npy)
            ntasks = npx*npy

        # Test both single and double precision code:
        for model in ['microhh', 'microhh_single']:
            print('Testing executable {}'.format(model))

            # Remove in/output previous runs and create input
            clean_workdir()
            create_profiles()

            # Run cold start
            execute('mpirun -n {} ./{} init drycblles'.format(ntasks, model), silent)
            execute('mpirun -n {} ./{} run drycblles' .format(ntasks, model), silent)

            # Rename output for comparison later
            variables = ['u', 'v', 'w', 'th', 'time']
            for var in variables:
                os.rename('{0:}.0003600'.format(var), '{0:}.0003600ref'.format(var))

            # Run warm start
            execute('mpirun -n {} ./{} run drycblles_restart'.format(ntasks, model), silent)

            # Compare restart files
            for var in variables:
                if not files_are_identical('{0:}.0003600'.format(var), '{0:}.0003600ref'.format(var)):
                    print('Executable {} in {} mode for variable {} not bitwise identical!'.format(model, mode, var))
                    great_success = False

    # Set namelists back to serial setup
    update_namelists('serial', 1, 1)

    # Return correct exit code to let Travis fail/succeed
    if great_success:
        print('ALL TESTS PASSED!')
        sys.exit(0)
    else:
        print('ONE OR MORE TESTS FAILED')
        sys.exit(1)
