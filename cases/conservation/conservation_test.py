import sys
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt

sys.path.append('../../python/')
import microhh_tools as mht
from . import conservation_funcs

no_opts = {}

opt_mpi = {
        'master': {'npx': 2, 'npy': 2}}

dict_rk = {
    'rk3': {'time': {'rkorder': 3}},
    'rk4': {'time': {'rkorder': 4}}}

dict_dt = {
    'dt1000': {'time': {'dtmax': 10.  }, 'stats': {'sampletime': 10.}},
    'dt0500': {'time': {'dtmax':  5.  }, 'stats': {'sampletime': 10.}},
    'dt0250': {'time': {'dtmax':  2.5 }, 'stats': {'sampletime': 10.}},
    'dt0125': {'time': {'dtmax':  1.25}, 'stats': {'sampletime': 10.}}}


def check_convergence(case_name, case_dir, experiment_name):
    """
    Check conservations and convergence of time integration schemes.
    """
    dts = np.array([10., 5., 2.5, 1.25])

    rk3 = conservation_funcs.Parse_conservation(case_dir, experiment_name, 'rk3')
    rk4 = conservation_funcs.Parse_conservation(case_dir, experiment_name, 'rk4')

    def convergence(errors):
        return (np.log(errors[-1]) - np.log(errors[0])) / (np.log(dts[-1]) - np.log(dts[0]))

    print(f'Convergence TKE RK3 = {convergence(np.abs(rk3.tke_loss))}')
    print(f'Convergence TKE RK4 = {convergence(np.abs(rk4.tke_loss))}')


def run_test(executable='microhh', prec='dp', mode='cpu', case_dir='.', experiment='local'):

    case_name = 'conservation'

    mht.run_permutations(
            case_name, no_opts, opt_mpi, [dict_rk, dict_dt],
            executable=executable, mode=mode, case_dir=case_dir, experiment=experiment)

    check_convergence(case_name, case_dir, experiment)


if __name__ == '__main__':

    run_test()