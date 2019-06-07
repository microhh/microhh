import sys
import shutil
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
sys.path.append('../../python/')
import microhh_tools as mht

def plot(filename='results.pdf'):
    data100_3rd = np.loadtxt('rkorder3dtmax10/conservation.out', skiprows=1)
    data200_3rd = np.loadtxt('rkorder3dtmax5/conservation.out', skiprows=1)
    data400_3rd = np.loadtxt('rkorder3dtmax25/conservation.out', skiprows=1)
    data800_3rd = np.loadtxt('rkorder3dtmax125/conservation.out', skiprows=1)

    data100_4th = np.loadtxt('rkorder4dtmax10/conservation.out', skiprows=1)
    data200_4th = np.loadtxt('rkorder4dtmax5/conservation.out', skiprows=1)
    data400_4th = np.loadtxt('rkorder4dtmax25/conservation.out', skiprows=1)
    data800_4th = np.loadtxt('rkorder4dtmax125/conservation.out', skiprows=1)

    time100_3rd = data100_3rd[:,1]
    mom100_3rd  = data100_3rd[:,7] / data100_3rd[1,7]
    tke100_3rd  = data100_3rd[:,8] / data100_3rd[1,8]
    mass100_3rd = data100_3rd[:,9] / data100_3rd[1,9]

    time200_3rd = data200_3rd[:,1]
    mom200_3rd  = data200_3rd[:,7] / data200_3rd[1,7]
    tke200_3rd  = data200_3rd[:,8] / data200_3rd[1,8]
    mass200_3rd = data200_3rd[:,9] / data200_3rd[1,9]

    time400_3rd = data400_3rd[:,1]
    mom400_3rd  = data400_3rd[:,7] / data400_3rd[1,7]
    tke400_3rd  = data400_3rd[:,8] / data400_3rd[1,8]
    mass400_3rd = data400_3rd[:,9] / data400_3rd[1,9]

    time800_3rd = data800_3rd[:,1]
    mom800_3rd  = data800_3rd[:,7] / data800_3rd[1,7]
    tke800_3rd  = data800_3rd[:,8] / data800_3rd[1,8]
    mass800_3rd = data800_3rd[:,9] / data800_3rd[1,9]

    time100_4th = data100_4th[:,1]
    mom100_4th  = data100_4th[:,7] / data100_4th[1,7]
    tke100_4th  = data100_4th[:,8] / data100_4th[1,8]
    mass100_4th = data100_4th[:,9] / data100_4th[1,9]

    time200_4th = data200_4th[:,1]
    mom200_4th  = data200_4th[:,7] / data200_4th[1,7]
    tke200_4th  = data200_4th[:,8] / data200_4th[1,8]
    mass200_4th = data200_4th[:,9] / data200_4th[1,9]

    time400_4th = data400_4th[:,1]
    mom400_4th  = data400_4th[:,7] / data400_4th[1,7]
    tke400_4th  = data400_4th[:,8] / data400_4th[1,8]
    mass400_4th = data400_4th[:,9] / data400_4th[1,9]

    time800_4th = data800_4th[:,1]
    mom800_4th  = data800_4th[:,7] / data800_4th[1,7]
    tke800_4th  = data800_4th[:,8] / data800_4th[1,8]
    mass800_4th = data800_4th[:,9] / data800_4th[1,9]

    with PdfPages(filename) as pdf:
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

executable = 'microhh'

kwargs = {'rkorder' : [3, 4], 'dtmax' : [10, 5, 2.5, 1.25]}
cases = mht.generator_parameter_change([mht.Case('conservation', casedir='.', keep=True)], **kwargs )

mht.test_cases(cases, executable, outputfile='conservation.csv')

plot('conservation.pdf')

for case in cases:
    if case.success:
        shutil.rmtree(case.rundir)

