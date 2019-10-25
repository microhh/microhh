import sys
import os
import copy
import shutil

sys.path.append('../../python/')
import microhh_tools as mht


def run_test(executable='microhh', float_type='dp', casedir='.', experiment='local'):
    base_case = mht.Case('drycblles', casedir=casedir, keep=True)
    mht.test_cases(
            [ base_case ],
            executable,
            experiment,
            outputfile='{}/drycblles_{}.csv'.format(casedir, experiment))

if __name__ == '__main__':
    if len(sys.argv) > 1:
        run_test(sys.argv[1:])
    else:
        run_test()
