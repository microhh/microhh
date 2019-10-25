import sys
import os
import copy
import shutil

sys.path.append('../../python/')
import microhh_tools as mht

dict_mpi = { 'default_run': { 'master': { 'npx': 2, 'npy': 4 } } }
dict_mpi_small = { 'default_run': { 'master': { 'npx': 2, 'npy': 4 }, 'grid': { 'itot': 64, 'jtot': 48, 'ktot': 32 } } }
dict_small = { 'default_run': { 'grid': { 'itot': 64, 'jtot': 48, 'ktot': 32 } } }

def run_test(executable='microhh', float_type='dp', mode='cpu', casedir='.', experiment='local'):
    base_case = mht.Case('moser180', casedir=casedir, rundir='default_run_{}'.format(experiment))
    if mode == 'cpumpi':
        cases = mht.generator_parameter_permutations(base_case, experiment, [dict_mpi])
    else:
        cases = [base_case]

    mht.run_cases(
            cases,
            executable,
            mode,
            outputfile='{}/moser180_{}.csv'.format(casedir, experiment))


def run_restart_test(executable='microhh', float_type='dp', mode='cpu', casedir='.', experiment='local'):
    base_case = mht.Case('moser180', casedir=casedir)
    if mode == 'cpumpi':
        base_case = mht.generator_parameter_permutations(base_case, experiment, [dict_mpi_small])[0]
    else:
        base_case = mht.generator_parameter_permutations(base_case, experiment, [dict_small])[0]

    cases = mht.generator_restart(base_case, experiment, 600.)

    mht.run_cases(
            cases,
            executable,
            mode,
            outputfile='{}/moser180_restart_{}.csv'.format(casedir, experiment))


if __name__ == '__main__':
    if len(sys.argv) > 1:
        run_test(sys.argv[1:])
    else:
        run_test()
