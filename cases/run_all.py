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
import numpy as np
import netCDF4 as nc
sys.path.append('../python/')
import microhh_tools as mht

from messages import *


def determine_mode():
    namelist = mht.Read_namelist()['master']

    npx = namelist['npx'] if 'npx' in namelist.keys() else 1
    npy = namelist['npy'] if 'npy' in namelist.keys() else 1
    mode = 'serial' if npx*npy == 1 else 'parallel'
    return mode, npx*npy


def run_scripts(scripts):
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

def compare(origin, file, starttime=-1, vars={}):
    nerror = 0
    try:
        nc_new = nc.Dataset(file, mode="r")
        nc_old = nc.Dataset('../'+origin+'/'+file, mode="r")
    except Exception:
        pass
        return sys.exc_info()[1]

    blacklist=['iter']
    rtol=1e-3
    atol=1e-8
    if len(vars) == 0:
        for key in nc_new.variables.keys():
            if key not in blacklist:
                vars.update({key: [rtol, atol]})

    for key, opts in vars.items():
        tmp = nc_new.variables[key][:]
        var_new = np.mean(nc_new.variables[key][starttime:,...],axis=0)
        var_old = np.mean(nc_old.variables[key][starttime:,...],axis=0)
        if not np.allclose(var_new, var_old, rtol=opts[0], atol=opts[1], equal_nan=True):
            with np.errstate(all='ignore'):
                print('{0} in {1} has a relative error of up to {2:.2%}'.format(key, file, np.max(np.abs((var_new - var_old)/var_old))))
            nerror = 1
    return nerror


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


def test_cases(cases, executable, outputfile=''):
    nerror = 0

    if not os.path.exists(executable):
        raise Exception('ERROR: Executable {} does not exists'.format(executable))

    # Get the absolute path to the executable
    executable_rel = executable
    executable = os.path.abspath(executable)
    rootdir = os.getcwd()

    for case in cases:
        print_header('Testing case \'{}\' for executable \'{}\''.format(case.name, executable_rel))
        # Determine whether to run serial or parallel
        nerror = 0
        # Move to working directory
        rootdir = os.getcwd()
        if case.rundir is '':
            print_error(case.name + ' case.rundir is empty; not allowed!')
            continue
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
            print_warning(case.name +  ': Cannot find {} for copying,  skipping case!'.format(casedir+fname))
            os.chdir(rootdir)
            continue

        # Update .ini file for testing
        for variable, value in case.options.items():
            mht.replace_namelist_value(variable, value, '{0}.ini'.format(case.name))
        mode, ntasks = determine_mode()

        # Create input data, and do other pre-processing
        nerror += run_scripts(case.pre)
        for phase in case.phases:
            case.time = timeit.default_timer()
            if mode == 'serial':
                nerror += execute('{} {} {}'.format(executable, phase, case.name))
            elif mode == 'parallel':
                nerror += execute('mpirun -n {} {} {} {}'.format(ntasks, executable, phase, case.name))
            case.time = timeit.default_timer() - case.time

        # Run the post-processing steps
        nerror += run_scripts(case.post)

        case.success = (nerror == 0)
        if not case.success:
            nerror += 1
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
        rundir  = rootdir + '/' + case.name +  '/' + case.rundir + '/'
        if case.success and not case.keep:
            shutil.rmtree(rundir)

    return nerror

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

def generator_scaling(cases, procs, type='strong', dir='y'):
    cases_out = []
    for case in cases:
        if type == 'weak':
            nl = mht.Read_namelist('{0}/{0}.ini'.format(case.name))
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
        if self.pre == {}:
            self.pre = {'{}_input.py'.format(name): None}
        if self.files == []:
            self.files = ['{0}.ini'.format(name), '{}_input.py'.format(name)]

def taylorgreen(executable, float_type):
    kwargs = {'itot' : [16, 32, 64, 128, 256], 'swadvec' : ['2', '4', '4m']}
    cases = generator_parameter_change([Case('taylorgreen', keep=True)], **kwargs )
    for case in cases:
        options = copy.deepcopy(case.options)
        for key, value in options.items():
            if key == 'itot':
                case.options.update({'ktot' : int(value/2)})
            elif key == 'swadvec':
                if value == '2':
                    case.options.update({'swspatialorder' : 2})
                else:
                    case.options.update({'swspatialorder' : 4})

    nerror = test_cases(cases,executable,outputfile='taylorgreen.csv')
    os.chdir('taylorgreen')
    import taylorgreen.taylorgreenconv as conv
    conv.main(float_type)

    for case in cases:
        if case.success:
            shutil.rmtree(case.rundir)
    os.chdir('..')

    return nerror

def conservation(executable, float_type):
    kwargs = {'rkorder' : [3, 4], 'dtmax' : [10, 5, 2.5, 1.25]}
    cases = generator_parameter_change([Case('conservation', keep=True)], **kwargs )

    nerror = test_cases(cases,executable,outputfile='conservation.csv')
    os.chdir('conservation')
    import conservation.conservationplot as cons
    cons.main()

    #for case in cases:
        #if case.success:
            #shutil.rmtree(case.rundir)
    os.chdir('..')

    return nerror

def compare_execs(cases, execs, options=None):
    nerror = 0

    for i in range(len(execs)):
        runcases = copy.deepcopy(cases)
        rundir = 'run{0:03d}'.format(i)
        for case in runcases:
            if i == 0:
                case.keep = True
            else:
                case.keep = False
                compfile = '{}_default_0000000.nc'.format(case.name)
                case.post = {__file__ : [['compare', 'run000', '{}_default_0000000.nc'.format(case.name)]]}
            case.rundir = rundir

            if options is not None:
                case.options.update(options[i])

        nerror += test_cases(runcases, execs[i],outputfile='compare_execs_{0:03d}.csv'.format(i))


    if nerror == 0:
        for case in cases:
            if case.success:
                shutil.rmtree(case.name+'/run000')

    return nerror

if __name__ == '__main__':

    nerror = 0

    exec_gpu_single = '../microhh_gpu_single'
    exec_gpu_double = '../microhh_gpu_double'
    exec_cpu_single = '../microhh_cpu_single'
    exec_cpu_double = '../microhh_cpu_double'
    exec_mpi_single = '../microhh_mpi_single'
    exec_mpi_double = '../microhh_mpi_double'
    float_type = 'double'

    if (False):
        # Serial
        #cases = generator_restart([Case('bomex')])
        cases = generator_scaling([Case('bomex')],procs=[1,2,4,8], type='weak',dir='xy')
        nerror += test_cases(cases, '../build_mpi/microhh')

    if (False):
        nerror += taylorgreen(exec_cpu_single, 'float')
    if (True):
        nerror += conservation(exec_gpu_single, 'float')


    if (False):
        blacklist = ['prandtlslope','moser600']

        # Settings for all test cases:

        dirs  = glob.glob('*')
        cases = []
        for dir in dirs:
            if os.path.isdir(dir) and dir not in blacklist:
                cases.append(Case(dir))
        opt_mpi = {'npx' : 8}
        cases = [Case('bomex'), Case('drycbl'), Case('drycblles')]
        nerror += compare_execs(cases, [exec_mpi_double, exec_mpi_single, exec_gpu_double, exec_gpu_single], options=[opt_mpi, opt_mpi,{} , {}])


    sys.exit(nerror)

#    if (False):
#        # DANGER: checkout all ini files
#        dirs  = glob.glob('*')
#        for dir in dirs:
#            if os.path.isdir(dir):
#                execute('git checkout {0}/{0}.ini'.format(dir))
