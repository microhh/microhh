import sys
import os
import copy
import shutil

sys.path.append('../../python/')
import microhh_tools as mht

list_resolution = [
        { 'grid' : { 'itot' :  16, 'ktot' :   8 } },
        { 'grid' : { 'itot' :  32, 'ktot' :  16 } },
        { 'grid' : { 'itot' :  64, 'ktot' :  32 } },
        { 'grid' : { 'itot' : 128, 'ktot' :  64 } },
        { 'grid' : { 'itot' : 256, 'ktot' : 128 } } ]

list_order = [
        { 'grid' : { 'swspatialorder', 2 }, 'advec' : { 'swadvec', '2'  } },
        { 'grid' : { 'swspatialorder', 4 }, 'advec' : { 'swadvec', '4'  } },
        { 'grid' : { 'swspatialorder', 4 }, 'advec' : { 'swadvec', '4m' } } ]

def run_test(executable='microhh', float_type='dp', casedir='.'):
    base_case = mht.Case('taylorgreen', casedir=casedir, keep=True)
    cases = mht.generator_parameter_change(base_case, [ list_resolution, list_order ])
    mht.test_cases(cases, executable, outputfile='taylorgreen.csv')

if __name__ == "__main__":
    if len(sys.argv)>1:
        run_test(sys.argv[1:])
    else:
        run_test()

