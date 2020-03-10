import sys
import os
import copy
import shutil

sys.path.append('../../python/')
import microhh_tools as mht

# Case configuration dicts
opt_mpi = {
        'master': {'npx': 2, 'npy': 2}}

opt_small = {
        'grid': {'itot': 64, 'jtot': 48, 'ktot': 32},
        'time': {'endtime': 10}}

opt_small_restart = {
        'grid': {'itot': 32, 'jtot': 16, 'ktot': 32},
        'time': {'endtime': 200, 'savetime': 100}}

# Case configuration dicts with name label for permutations.
dict_opts = {
        'all_enabled': {},
        'advec': {'advec': {'swadvec': 0}},
        'diff': {'diff': {'swdiff': 0}}}

list_permutations = [ dict_opts ]


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

    # Deep copy the small version of the reference case and disable stats.
    options = opt_small.copy()

    if mode == 'cpumpi':
        mht.merge_options(options, opt_mpi)

    base_case = mht.Case(
        'moser180',
        casedir=casedir,
        rundir=experiment,
        options=options)

    base_cases = mht.generator_parameter_permutations(base_case, list_permutations)
    cases = []
    for case in base_cases:
        cases.extend(mht.generator_restart(case, 300.))

    mht.run_cases(
        cases,
        executable,
        mode,
        outputfile='{}/bomex_restart_{}.csv'.format(casedir, experiment))

if __name__ == '__main__':
    kwargs = dict([arg.split('=') for arg in sys.argv[2:]])
    if len(sys.argv) > 1:
        globals()[sys.argv[1]](**kwargs)
    else:
        run()
