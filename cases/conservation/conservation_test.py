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
    'dt1000': {'time': {'dt': 10.  }},
    'dt0500': {'time': {'dt':  5.  }},
    'dt0250': {'time': {'dt':  2.5 }},
    'dt0125': {'time': {'dt':  1.25}}}


def get_data_from_nc(experiment_name):
    data = nc.Dataset('{}/conservation.default.0000000.nc'.format(experiment_name), 'r')
    time = data.variables['time'][:]
    mom = data.group['default'].variable['u'][:,:].sum(axis=1)
    tke = data.group['budget'].variable['tke'][:,:].sum(axis=1)
    mass = data.group['default'].variable['s'][:,:].sum(axis=1)
    return time, mom, tke, mass


def plot(experiment_name):
    time100_3rd, mom100_3rd, tke100_3rd, mass100_3rd = get_data_from_nc('{}_rk3_dt1000'.format(experiment_name))
    time200_3rd, mom200_3rd, tke200_3rd, mass200_3rd = get_data_from_nc('{}_rk3_dt0500'.format(experiment_name))
    time400_3rd, mom400_3rd, tke400_3rd, mass400_3rd = get_data_from_nc('{}_rk3_dt0250'.format(experiment_name))
    time800_3rd, mom800_3rd, tke800_3rd, mass800_3rd = get_data_from_nc('{}_rk3_dt0125'.format(experiment_name))

    time100_4th, mom100_4th, tke100_4th, mass100_4th = get_data_from_nc('{}_rk4_dt1000'.format(experiment_name))
    time200_4th, mom200_4th, tke200_4th, mass200_4th = get_data_from_nc('{}_rk4_dt0500'.format(experiment_name))
    time400_4th, mom400_4th, tke400_4th, mass400_4th = get_data_from_nc('{}_rk4_dt0250'.format(experiment_name))
    time800_4th, mom800_4th, tke800_4th, mass800_4th = get_data_from_nc('{}_rk4_dt0125'.format(experiment_name))

    file_name = '{}.pdf'.format(case_name)

    with PdfPages(file_name) as pdf:
        plt.figure()
        plt.subplot(131)
        plt.plot(time100_3rd[1:100], mom100_3rd [1:100] - mom100_3rd [1], 'b-' , label='rk3: dt=10.0')
        plt.plot(time100_4th[1:100], mom100_4th [1:100] - mom100_4th [1], 'r-' , label='rk4: dt=10.0')
        plt.plot(time200_3rd[1:200], mom200_3rd [1:200] - mom200_3rd [1], 'b--', label='rk3: dt=5.0')
        plt.plot(time200_4th[1:200], mom200_4th [1:200] - mom200_4th [1], 'r--', label='rk4: dt=5.0')
        plt.plot(time400_3rd[1:400], mom400_3rd [1:400] - mom400_3rd [1], 'b:' , label='rk3: dt=2.5')
        plt.plot(time400_4th[1:400], mom400_4th [1:400] - mom400_4th [1], 'r:' , label='rk4: dt=2.5')
        plt.plot(time800_3rd[1:800], mom800_3rd [1:800] - mom800_3rd [1], 'b:' , label='rk3: dt=1.25')
        plt.plot(time800_4th[1:800], mom800_4th [1:800] - mom800_4th [1], 'r:' , label='rk4: dt=1.25')
        plt.legend(loc=0, frameon=False)
        plt.title('momentum conservation')
        plt.subplot(132)
        plt.plot(time100_3rd[1:100], tke100_3rd [1:100] - tke100_3rd [1], 'b-' , label='rk3: dt=10.0')
        plt.plot(time100_4th[1:100], tke100_4th [1:100] - tke100_4th [1], 'r-' , label='rk4: dt=10.0')
        plt.plot(time200_3rd[1:200], tke200_3rd [1:200] - tke200_3rd [1], 'b--', label='rk3: dt=5.0')
        plt.plot(time200_4th[1:200], tke200_4th [1:200] - tke200_4th [1], 'r--', label='rk4: dt=5.0')
        plt.plot(time400_3rd[1:400], tke400_3rd [1:400] - tke400_3rd [1], 'b:' , label='rk3: dt=2.5')
        plt.plot(time400_4th[1:400], tke400_4th [1:400] - tke400_4th [1], 'r:' , label='rk4: dt=2.5')
        plt.plot(time800_3rd[1:800], tke800_3rd [1:800] - tke800_3rd [1], 'b:' , label='rk3: dt=1.25')
        plt.plot(time800_4th[1:800], tke800_4th [1:800] - tke800_4th [1], 'r:' , label='rk4: dt=1.25')
        plt.legend(loc=0, frameon=False)
        plt.title('energy conservation')
        plt.subplot(133)
        plt.plot(time100_3rd[1:100], mass100_3rd[1:100] - mass100_3rd[1], 'b-' , label='rk3: dt=10.0')
        plt.plot(time100_4th[1:100], mass100_4th[1:100] - mass100_4th[1], 'r-' , label='rk4: dt=10.0')
        plt.plot(time200_3rd[1:200], mass200_3rd[1:200] - mass200_3rd[1], 'b--', label='rk3: dt=5.0')
        plt.plot(time200_4th[1:200], mass200_4th[1:200] - mass200_4th[1], 'r--', label='rk4: dt=5.0')
        plt.plot(time400_3rd[1:400], mass400_3rd[1:400] - mass400_3rd[1], 'b:' , label='rk3: dt=2.5')
        plt.plot(time400_4th[1:400], mass400_4th[1:400] - mass400_4th[1], 'r:' , label='rk4: dt=2.5')
        plt.plot(time800_3rd[1:800], mass800_3rd[1:800] - mass800_3rd[1], 'b:' , label='rk3: dt=1.25')
        plt.plot(time800_4th[1:800], mass800_4th[1:800] - mass800_4th[1], 'r:' , label='rk4: dt=1.25')
        plt.legend(loc=0, frameon=False)
        plt.title('mass conservation')
        pdf.savefig()

        tkeloss100_3rd = tke100_3rd[100] - tke100_3rd[1]
        tkeloss200_3rd = tke200_3rd[200] - tke100_3rd[1]
        tkeloss400_3rd = tke400_3rd[400] - tke100_3rd[1]
        tkeloss800_3rd = tke800_3rd[800] - tke100_3rd[1]
        tkeloss100_4th = tke100_4th[100] - tke100_4th[1]
        tkeloss200_4th = tke200_4th[200] - tke100_4th[1]
        tkeloss400_4th = tke400_4th[400] - tke100_4th[1]
        tkeloss800_4th = tke800_4th[800] - tke100_4th[1]

        dts           = np.array([10., 5., 2.5, 1.25])
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


if __name__ == '__main__':

    case_name = 'conservation'

    kwargs = dict([arg.split('=') for arg in sys.argv[2:]])

    if len(sys.argv) > 1:
        function_name = sys.argv[1]

        if function_name == 'run_case':
            mht.run_permutations(case_name, no_opts, opt_mpi, [dict_rk, dict_dt], **kwargs)
        else:
            raise Exception('\"{}\" is an invalid option'.format(function_name))

    else:
        mht.run_permutations(case_name, no_opts, opt_mpi, [dict_rk, dict_dt])

    experiment_name = 'local'

    plot(experiment_name)
