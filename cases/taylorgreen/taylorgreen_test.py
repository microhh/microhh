import sys
import os
import copy
import shutil

sys.path.append('../../python/')
import microhh_tools as mht

dict_resolution = {
        'itot016': { 'grid': { 'itot':  16, 'ktot':   8 } },
        'itot032': { 'grid': { 'itot':  32, 'ktot':  16 } },
        'itot064': { 'grid': { 'itot':  64, 'ktot':  32 } },
        'itot128': { 'grid': { 'itot': 128, 'ktot':  64 } },
        'itot256': { 'grid': { 'itot': 256, 'ktot': 128 } } }

dict_order = {
        'swadvec2' : { 'grid': { 'swspatialorder': 2 }, 'advec': { 'swadvec': '2'  } },
        'swadvec4' : { 'grid': { 'swspatialorder': 4 }, 'advec': { 'swadvec': '4'  } },
        'swadvec4m': { 'grid': { 'swspatialorder': 4 }, 'advec': { 'swadvec': '4m' } } }

def run_test(executable='microhh', float_type='dp', casedir='.'):
    base_case = mht.Case('taylorgreen', casedir=casedir, keep=True)
    cases = mht.generator_parameter_permutations(base_case, [ dict_resolution, dict_order ])
    mht.test_cases(cases, executable, outputfile='taylorgreen.csv')

if __name__ == '__main__':
    if len(sys.argv) > 1:
        run_test(sys.argv[1:])
    else:
        run_test()
