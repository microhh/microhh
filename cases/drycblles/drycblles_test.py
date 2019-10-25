import sys
import os
import copy
import shutil

sys.path.append('../../python/')
import microhh_tools as mht

dict_mpi = { 'default_run': { 'master': { 'npx': 2, 'npy': 4 } } }

def run_test(executable='microhh', float_type='dp', mode='cpu', casedir='.', experiment='local'):
    base_case = mht.Case('drycblles', casedir=casedir, keep=True)
    if mode == 'cpumpi':
        cases = mht.generator_parameter_permutations(base_case, [dict_mpi])
    else:
        cases = [base_case]

    mht.run_cases(
            cases,
            executable,
            mode,
            experiment,
            outputfile='{}/drycblles_{}.csv'.format(casedir, experiment))


if __name__ == '__main__':
    if len(sys.argv) > 1:
        run_test(sys.argv[1:])
    else:
        run_test()
