import sys
import numpy as np

sys.path.append('../../python/')
import microhh_tools as mht
from . import taylorgreen_funcs

no_opts = {}

opt_mpi = {}

dict_resolution = {
    'itot016': {'grid': {'itot':  16, 'ktot':   8}},
    'itot032': {'grid': {'itot':  32, 'ktot':  16}},
    'itot064': {'grid': {'itot':  64, 'ktot':  32}},
    'itot128': {'grid': {'itot': 128, 'ktot':  64}},
    'itot256': {'grid': {'itot': 256, 'ktot': 128}}}

dict_order = {
    'swadvec2' : {'grid': {'swspatialorder': 2}, 'advec': {'swadvec': '2' }},
    'swadvec4' : {'grid': {'swspatialorder': 4}, 'advec': {'swadvec': '4' }},
    'swadvec4m': {'grid': {'swspatialorder': 4}, 'advec': {'swadvec': '4m'}}}


def check_convergence(case_name, case_dir, experiment):

    # Log the stream output to a file.
    sys.stdout = open('{}/{}_{}.log'.format(case_dir, case_name, experiment), 'w')

    ns = np.array([16, 32, 64, 128, 256])
    dxs = 1. / ns

    time = 1
    visc = (8. * np.pi**2. * 100.)**(-1.)

    def parse(name):
        cases = []
        err_u = []
        err_w = []
        err_p = []
        for n in ns:
            c = taylorgreen_funcs.Parse_TaylorGreen(time, visc, f'{case_dir}/{experiment}_itot{n:03d}_{name}')
            err_u.append(c.u_err)
            err_w.append(c.w_err)
            err_p.append(c.p_err)
            cases.append(c)

        return cases, err_u, err_w, err_p

    def convergence(errors):
        return (np.log(errors[-1]) - np.log(errors[0])) / (np.log(dxs[-1]) - np.log(dxs[0]))

    # 2nd order.
    advec2,  err_u_2,  err_w_2,  err_p_2  = parse('swadvec2')

    # 4th order with `advec4m`.
    advec4m, err_u_4m, err_w_4m, err_p_4m = parse('swadvec4m')

    # 4th order with `advec4`.
    advec4,  err_u_4,  err_w_4,  err_p_4  = parse('swadvec4')

    print(f'Errors p 2nd = {err_p_2}')
    print(f'Convergence u 2nd = {convergence(err_u_2)}')
    print(f'Convergence w 2nd = {convergence(err_w_2)}')
    print(f'Convergence p 2nd = {convergence(err_p_2)}')

    print(f'Errors p 4thm = {err_p_4m}')
    print(f'Convergence u 4thm = {convergence(err_u_4m)}')
    print(f'Convergence w 4thm = {convergence(err_w_4m)}')
    print(f'Convergence p 4thm = {convergence(err_p_4m)}')

    print(f'Errors p 4th = {err_p_4}')
    print(f'Convergence u 4th = {convergence(err_u_4)}')
    print(f'Convergence w 4th = {convergence(err_w_4)}')
    print(f'Convergence p 4th = {convergence(err_p_4)}')


def run_test(executable='microhh', prec='dp', mode='cpu', case_dir='.', experiment='local'):

    case_name = 'taylorgreen'

    mht.run_permutations(
            case_name, no_opts, opt_mpi, [dict_resolution, dict_order],
            executable=executable, mode=mode, case_dir=case_dir, experiment=experiment)

    check_convergence(case_name, case_dir, experiment)


if __name__ == '__main__':
    run_test()
