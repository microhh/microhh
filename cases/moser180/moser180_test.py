import sys
import os
import copy
import shutil

sys.path.append('../../python/')
import microhh_tools as mht

opt_mpi = [('master', 'npx', 2), ('master', 'npy', 4)]
opt_small = [('grid', 'itot', 64), ('grid', 'jtot', 48), ('grid', 'ktot', 32)]


def run(executable='microhh', mode='cpu', casedir='.', experiment='local'):
    options = []
    if mode == 'cpumpi':
        options.extend(opt_mpi)
    cases = [
        mht.Case(
            'moser180',
            casedir=casedir,
            rundir=experiment,
            options=options)]

    mht.run_cases(
        cases,
        executable,
        mode,
        outputfile='{}/moser180_{}.csv'.format(casedir, experiment))


def run_small(executable='microhh', mode='cpu',
              casedir='.', experiment='local'):
    options = opt_small.copy()
    if mode == 'cpumpi':
        options.extend(opt_mpi)
    cases = [
        mht.Case(
            'moser180',
            casedir=casedir,
            rundir='{}_small'.format(experiment),
            options=options)]

    mht.run_cases(
        cases,
        executable,
        mode,
        outputfile='{}/moser180_small_{}.csv'.format(casedir, experiment))


def run_restart(executable='microhh', mode='cpu',
                casedir='.', experiment='local'):
    options = opt_small.copy()
    if mode == 'cpumpi':
        options.extend(opt_mpi)
    base_case = mht.Case(
        'moser180',
        casedir=casedir,
        rundir=experiment,
        options=options)
    cases = mht.generator_restart(base_case, 60.)

    mht.run_cases(
        cases,
        executable,
        mode,
        outputfile='{}/moser180_restart_{}.csv'.format(casedir, experiment))


if __name__ == '__main__':
    kwargs = dict([arg.split('=') for arg in sys.argv[2:]])
    if len(sys.argv) > 1:
        globals()[sys.argv[1]](**kwargs)
    else:
        run()
