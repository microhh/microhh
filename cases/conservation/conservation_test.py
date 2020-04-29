import sys
import shutil
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sys.path.append('../../python/')
import microhh_tools as mht

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


def get_data_from_nc(case_dir, experiment_name):
    data = nc.Dataset('{}/{}/conservation.default.0000000.nc'.format(case_dir, experiment_name), 'r')
    time = data.variables['time'][:]
    dz = data.variables['zh'][1] - data.variables['zh'][0]

    mom = dz * (data.groups['default'].variables['u'][:,:] + data.groups['default'].variables['v'][:,:]).sum(axis=1)
    tke = 0.5 * dz * ( data.groups['default'].variables['u_2'][:,:] \
                     + data.groups['default'].variables['v_2'][:,:] \
                     + data.groups['default'].variables['w_2'][:,:-1]).sum(axis=1)
    mass = dz * (data.groups['default'].variables['s'][:,:]).sum(axis=1)

    mom  /= mom [1]
    tke  /= tke [1]
    mass /= mass[1]

    data_out = np.loadtxt('{}/{}/conservation.out'.format(case_dir, experiment_name), skiprows=1)
    time = data_out[:,1]
    mom = data_out[:,7] / data_out[1,7]
    tke = data_out[:,8] / data_out[1,8]
    mass = np.zeros_like(time)

    return time, mom, tke, mass


def plot(case_name, case_dir, experiment_name):
    time100_3rd, mom100_3rd, tke100_3rd, mass100_3rd = get_data_from_nc(case_dir, '{}_rk3_dt1000'.format(experiment_name))
    time200_3rd, mom200_3rd, tke200_3rd, mass200_3rd = get_data_from_nc(case_dir, '{}_rk3_dt0500'.format(experiment_name))
    time400_3rd, mom400_3rd, tke400_3rd, mass400_3rd = get_data_from_nc(case_dir, '{}_rk3_dt0250'.format(experiment_name))
    time800_3rd, mom800_3rd, tke800_3rd, mass800_3rd = get_data_from_nc(case_dir, '{}_rk3_dt0125'.format(experiment_name))
    time100_4th, mom100_4th, tke100_4th, mass100_4th = get_data_from_nc(case_dir, '{}_rk4_dt1000'.format(experiment_name))
    time200_4th, mom200_4th, tke200_4th, mass200_4th = get_data_from_nc(case_dir, '{}_rk4_dt0500'.format(experiment_name))
    time400_4th, mom400_4th, tke400_4th, mass400_4th = get_data_from_nc(case_dir, '{}_rk4_dt0250'.format(experiment_name))
    time800_4th, mom800_4th, tke800_4th, mass800_4th = get_data_from_nc(case_dir, '{}_rk4_dt0125'.format(experiment_name))

    file_name = '{0}/{1}_{2}.pdf'.format(case_dir, case_name, experiment_name)

    with PdfPages(file_name) as pdf:
        plt.figure()
        plt.plot(time100_3rd[1:], mom100_3rd [1:] - mom100_3rd [1], 'b-' , label='rk3: dt=10.0')
        plt.plot(time100_4th[1:], mom100_4th [1:] - mom100_4th [1], 'r-' , label='rk4: dt=10.0')
        plt.plot(time200_3rd[1:], mom200_3rd [1:] - mom200_3rd [1], 'b--', label='rk3: dt=5.0')
        plt.plot(time200_4th[1:], mom200_4th [1:] - mom200_4th [1], 'r--', label='rk4: dt=5.0')
        plt.plot(time400_3rd[1:], mom400_3rd [1:] - mom400_3rd [1], 'b:' , label='rk3: dt=2.5')
        plt.plot(time400_4th[1:], mom400_4th [1:] - mom400_4th [1], 'r:' , label='rk4: dt=2.5')
        plt.plot(time800_3rd[1:], mom800_3rd [1:] - mom800_3rd [1], 'b:' , label='rk3: dt=1.25')
        plt.plot(time800_4th[1:], mom800_4th [1:] - mom800_4th [1], 'r:' , label='rk4: dt=1.25')
        plt.legend(loc=0, frameon=False)
        plt.title('momentum')
        pdf.savefig()

        plt.figure()
        plt.plot(time100_3rd[1:], tke100_3rd [1:] - tke100_3rd [1], 'b-' , label='rk3: dt=10.0')
        plt.plot(time100_4th[1:], tke100_4th [1:] - tke100_4th [1], 'r-' , label='rk4: dt=10.0')
        plt.plot(time200_3rd[1:], tke200_3rd [1:] - tke200_3rd [1], 'b--', label='rk3: dt=5.0')
        plt.plot(time200_4th[1:], tke200_4th [1:] - tke200_4th [1], 'r--', label='rk4: dt=5.0')
        plt.plot(time400_3rd[1:], tke400_3rd [1:] - tke400_3rd [1], 'b:' , label='rk3: dt=2.5')
        plt.plot(time400_4th[1:], tke400_4th [1:] - tke400_4th [1], 'r:' , label='rk4: dt=2.5')
        plt.plot(time800_3rd[1:], tke800_3rd [1:] - tke800_3rd [1], 'b:' , label='rk3: dt=1.25')
        plt.plot(time800_4th[1:], tke800_4th [1:] - tke800_4th [1], 'r:' , label='rk4: dt=1.25')
        plt.legend(loc=0, frameon=False)
        plt.title('energy')
        pdf.savefig()

        plt.figure()
        plt.plot(time100_3rd[1:], mass100_3rd[1:] - mass100_3rd[1], 'b-' , label='rk3: dt=10.0')
        plt.plot(time100_4th[1:], mass100_4th[1:] - mass100_4th[1], 'r-' , label='rk4: dt=10.0')
        plt.plot(time200_3rd[1:], mass200_3rd[1:] - mass200_3rd[1], 'b--', label='rk3: dt=5.0')
        plt.plot(time200_4th[1:], mass200_4th[1:] - mass200_4th[1], 'r--', label='rk4: dt=5.0')
        plt.plot(time400_3rd[1:], mass400_3rd[1:] - mass400_3rd[1], 'b:' , label='rk3: dt=2.5')
        plt.plot(time400_4th[1:], mass400_4th[1:] - mass400_4th[1], 'r:' , label='rk4: dt=2.5')
        plt.plot(time800_3rd[1:], mass800_3rd[1:] - mass800_3rd[1], 'b:' , label='rk3: dt=1.25')
        plt.plot(time800_4th[1:], mass800_4th[1:] - mass800_4th[1], 'r:' , label='rk4: dt=1.25')
        plt.legend(loc=0, frameon=False)
        plt.title('mass')
        pdf.savefig()

        tkeloss100_3rd = tke100_3rd[-1] - tke100_3rd[1]
        tkeloss200_3rd = tke200_3rd[-1] - tke100_3rd[1]
        tkeloss400_3rd = tke400_3rd[-1] - tke100_3rd[1]
        tkeloss800_3rd = tke800_3rd[-1] - tke100_3rd[1]
        tkeloss100_4th = tke100_4th[-1] - tke100_4th[1]
        tkeloss200_4th = tke200_4th[-1] - tke100_4th[1]
        tkeloss400_4th = tke400_4th[-1] - tke100_4th[1]
        tkeloss800_4th = tke800_4th[-1] - tke100_4th[1]

        dts = np.array([10., 5., 2.5, 1.25])
        tkelosses_3rd = np.array([tkeloss100_3rd, tkeloss200_3rd, tkeloss400_3rd, tkeloss800_3rd])
        tkelosses_4th = np.array([tkeloss100_4th, tkeloss200_4th, tkeloss400_4th, tkeloss800_4th])

        off2 = 0.001
        off3 = 0.0001
        off4 = 0.00002
        slope2 = off2*(dts[:] / dts[0])**2.
        slope3 = off3*(dts[:] / dts[0])**3.
        slope4 = off4*(dts[:] / dts[0])**4.

        plt.figure()
        plt.loglog(dts, abs(tkelosses_3rd), 'bo-', label='rk3')
        plt.loglog(dts, abs(tkelosses_4th), 'go-', label='rk4')
        plt.loglog(dts, slope2, 'k--', label='2nd')
        plt.loglog(dts, slope3, 'k:' , label='3rd')
        plt.loglog(dts, slope4, 'k-.', label='4th')
        plt.xlabel('dt')
        plt.ylabel('tke loss')
        plt.legend(loc=0, frameon=False)
        pdf.savefig()


def run_test(executable='microhh', prec='dp', mode='cpu', case_dir='.', experiment='local'):

    case_name = 'conservation'

    mht.run_permutations(
            case_name, no_opts, opt_mpi, [dict_rk, dict_dt],
            executable=executable, mode=mode, case_dir=case_dir, experiment=experiment)

    plot(case_name, case_dir, experiment)


if __name__ == '__main__':

    run_test()
 
