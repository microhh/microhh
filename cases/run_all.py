import subprocess
import importlib
import shutil
import re
import os
import glob

def replace_namelist_value(namelist_file, variable, new_value):
    with open(namelist_file, "r") as source:
        lines = source.readlines()
    with open(namelist_file, "w") as source:
        for line in lines:
            source.write(re.sub(r'({}).*'.format(variable), r'\1={}'.format(new_value), line))


def determine_mode(namelist):
    npx = namelist['npx'] if 'npx' in namelist.keys() else 1
    npy = namelist['npx'] if 'npx' in namelist.keys() else 1
    mode = 'serial' if npx*npy == 1 else 'parallel'
    return mode, npx*npy


def run_scripts(case_name, scripts):
    if scripts is not None:
        for script, functions in scripts.items():

            # Import module as name minus the `.py` at the end
            module = script.replace('.py', '')

            # Import module; this executes all code that is not in classes/functions
            lib = importlib.import_module('{0}.{1}'.format(case_name, module))

            # Call specific routines in the module:
            if functions is not None:
                if isinstance(functions, list):
                    for function in functions:
                        rc = getattr(lib, function)()
                        if rc != 0:
                            print('ERROR: {}: {}() returned {}'.format(script, function, rc))
                else:
                    rc = getattr(lib, functions)()
                    if rc != 0:
                        print('ERROR: {}: {}() returned {}'.format(script, functions, rc))


def execute(command):
    sp = subprocess.Popen(command, executable='/bin/bash', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = sp.communicate()
    if sp.returncode != 0:
        with open('stdout.log', 'w') as  f:
            f.write(out.decode('utf-8'))
        with open('stderr.log', 'w') as  f:
            f.write(err.decode('utf-8'))

        print('ERROR: \"{}\" returned \"{}\". Output written to stdout.log and stderr.log'.format(command, sp.returncode))


def test_cases(cases, executable):
    if not os.path.exists(executable):
        raise Exception('ERROR: Executable {} does not exists'.format(executable))

    # Aaarghhh
    executable = '../{}'.format(executable)

    for case in cases:
        print('Testing case \"{}\" for executable \"{}\"'.format(case.name, executable))

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

            # Determine whether to run serial or parallel
            mode, ntasks = determine_mode(case.options)

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
            print('WARNING: cannot find {0}.ini and/or {0}prof.py, skipping...!'.format(case.name))

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

    if (True):
        # Serial
        cases = [
            Case("bomex",     { "itot" : 16, "jtot" : 16,  "ktot" : 32, "endtime": 3600 }),
            Case("drycblles", { "itot" : 16, "jtot" : 16, "ktot" : 32 }, post={'validate.py': ['test1', 'test2']})
                ]

        test_cases(cases, '../build/microhh')
        test_cases(cases, '../build/microhh_single')

    if (False):
        blacklist = ['prandtlslope']

        # Settings for all test cases:
        settings = { "itot" : 16, "jtot" : 16,  "ktot" : 16, "endtime": 10 }

        dirs  = glob.glob('*')
        cases = []
        for dir in dirs:
            if os.path.isdir(dir) and dir not in blacklist and 'restart' not in dir:
                cases.append(Case(dir, settings))

        test_cases(cases, '../build/microhh')
        test_cases(cases, '../build/microhh_single')

#    if (False):
#        #
#        # undo Barts mistakes..
#        #
#
#        dirs  = glob.glob('*')
#        for dir in dirs:
#            if os.path.isdir(dir):
#                execute('git checkout {0}/{0}.ini'.format(dir))
