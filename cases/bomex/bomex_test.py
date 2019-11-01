import sys
import os
import copy
import shutil

sys.path.append('../../python/')
import microhh_tools as mht

# Case configuration dicts
opt_mpi = {
        'master': {'npx': 2, 'npy': 4}}

opt_small = {
        'grid': {'itot': 8, 'jtot': 8},
        'time': {'endtime': 7200}}

# Case configuration dicts with name label for permutations.
dict_opts = {
        'all_enabled': {},
        'advec': {'advec': {'swadvec': 0}},
        'diff': {'diff': {'swdiff': 0}},
        'thermo': {'thermo': {'swthermo': '0'}},
        'basestate': {'thermo': {'swupdatebasestate': 0}},
        'buffer': {'buffer': {'swbuffer': 0}}}

list_permutations = [ dict_opts ]


def run(executable='microhh', mode='cpu', casedir='.', experiment='local'):

    options = []
    if mode == 'cpumpi':
        options.extend(opt_mpi)

    cases = [
        mht.Case(
            'bomex',
            casedir=casedir,
            rundir=experiment,
            options=options)]

    mht.run_cases(
        cases,
        executable,
        mode,
        outputfile='{}/bomex_{}.csv'.format(casedir, experiment))


def run_small(executable='microhh', mode='cpu',
              casedir='.', experiment='local'):

    # Deep copy the reference case.
    options = opt_small.copy()
    if mode == 'cpumpi':
        mht.merge_options(options, opt_mpi)

    cases = [
        mht.Case(
            'bomex',
            casedir=casedir,
            rundir='{}_small'.format(experiment),
            options=options)]

    mht.run_cases(
        cases,
        executable,
        mode,
        outputfile='{}/bomex_small_{}.csv'.format(casedir, experiment))


def run_restart(executable='microhh', mode='cpu',
                casedir='.', experiment='local'):

    # Deep copy the small version of the reference case.
    options = opt_small.copy()

    if mode == 'cpumpi':
        mht.merge_options(options, opt_mpi)

    base_case = mht.Case(
        'bomex',
        casedir=casedir,
        rundir=experiment,
        options=options)

    base_cases = mht.generator_parameter_permutations(base_case, list_permutations)
    cases = []
    for case in base_cases:
        cases.extend(mht.generator_restart(case, 1800.))

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
