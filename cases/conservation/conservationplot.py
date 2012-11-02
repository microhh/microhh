from pylab import *

data100_2nd = loadtxt('conservation100_2nd/conservation.out', skiprows=1)
data200_2nd = loadtxt('conservation200_2nd/conservation.out', skiprows=1)
data400_2nd = loadtxt('conservation400_2nd/conservation.out', skiprows=1)

data100_4th = loadtxt('conservation100_4th/conservation.out', skiprows=1)
data200_4th = loadtxt('conservation200_4th/conservation.out', skiprows=1)
data400_4th = loadtxt('conservation400_4th/conservation.out', skiprows=1)

time100_2nd = data100_2nd[:,1]
mom100_2nd  = data100_2nd[:,7] / data100_2nd[1,7]
tke100_2nd  = data100_2nd[:,8] / data100_2nd[1,8]
mass100_2nd = data100_2nd[:,9] / data100_2nd[1,9]

time200_2nd = data200_2nd[:,1]
mom200_2nd  = data200_2nd[:,7] / data200_2nd[1,7]
tke200_2nd  = data200_2nd[:,8] / data200_2nd[1,8]
mass200_2nd = data200_2nd[:,9] / data200_2nd[1,9]

time400_2nd = data400_2nd[:,1]
mom400_2nd  = data400_2nd[:,7] / data400_2nd[1,7]
tke400_2nd  = data400_2nd[:,8] / data400_2nd[1,8]
mass400_2nd = data400_2nd[:,9] / data400_2nd[1,9]

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

figure()
subplot(131)
plot(time100_2nd[1: 500], mom100_2nd [1: 500] - mom100_2nd [1], 'b-' , label='cfl08_2nd')
plot(time100_4th[1: 500], mom100_4th [1: 500] - mom100_4th [1], 'r-' , label='cfl08_4th')
plot(time200_2nd[1:1000], mom200_2nd [1:1000] - mom200_2nd [1], 'b--', label='cfl04_2nd')
plot(time200_4th[1:1000], mom200_4th [1:1000] - mom200_4th [1], 'r--', label='cfl04_4th')
plot(time400_2nd[1:2000], mom400_2nd [1:2000] - mom400_2nd [1], 'b:' , label='cfl02_2nd')
plot(time400_4th[1:2000], mom400_4th [1:2000] - mom400_4th [1], 'r:' , label='cfl02_4th')
legend(loc=0, frameon=False)
title('momentum conservation')
subplot(132)
plot(time100_2nd[1: 500], tke100_2nd [1: 500] - tke100_2nd [1], 'b-' , label='cfl08_2nd')
plot(time100_4th[1: 500], tke100_4th [1: 500] - tke100_4th [1], 'r-' , label='cfl08_4th')
plot(time200_2nd[1:1000], tke200_2nd [1:1000] - tke200_2nd [1], 'b--', label='cfl04_2nd')
plot(time200_4th[1:1000], tke200_4th [1:1000] - tke200_4th [1], 'r--', label='cfl04_4th')
plot(time400_2nd[1:2000], tke400_2nd [1:2000] - tke400_2nd [1], 'b:' , label='cfl02_2nd')
plot(time400_4th[1:2000], tke400_4th [1:2000] - tke400_4th [1], 'r:' , label='cfl02_4th')
legend(loc=0, frameon=False)
title('energy conservation')
subplot(133)
plot(time100_2nd[1: 500], mass100_2nd[1: 500] - mass100_2nd[1], 'b-' , label='cfl08_2nd')
plot(time100_4th[1: 500], mass100_4th[1: 500] - mass100_4th[1], 'r-' , label='cfl08_4th')
plot(time200_2nd[1:1000], mass200_2nd[1:1000] - mass200_2nd[1], 'b--', label='cfl04_2nd')
plot(time200_4th[1:1000], mass200_4th[1:1000] - mass200_4th[1], 'r--', label='cfl04_4th')
plot(time400_2nd[1:2000], mass400_2nd[1:2000] - mass400_2nd[1], 'b:' , label='cfl02_2nd')
plot(time400_4th[1:2000], mass400_4th[1:2000] - mass400_4th[1], 'r:' , label='cfl02_4th')
legend(loc=0, frameon=False)
title('mass conservation')

tkeloss100_2nd = tke100_2nd[100] - tke100_2nd[1]
tkeloss200_2nd = tke200_2nd[200] - tke100_2nd[1]
tkeloss400_2nd = tke400_2nd[400] - tke100_2nd[1]
tkeloss100_4th = tke100_4th[100] - tke100_4th[1]
tkeloss200_4th = tke200_4th[200] - tke100_4th[1]
tkeloss400_4th = tke400_4th[400] - tke100_4th[1]

dts           = array([8., 4., 2.])
tkelosses_2nd = array([tkeloss100_2nd, tkeloss200_2nd, tkeloss400_2nd])
tkelosses_4th = array([tkeloss100_4th, tkeloss200_4th, tkeloss400_4th])

off2 = 0.0001
off3 = 0.00005
off4 = 0.00002
slope2 = off2*(dts[:] / dts[0])**2.
slope3 = off3*(dts[:] / dts[0])**3.
slope4 = off4*(dts[:] / dts[0])**4.


figure()
loglog(dts, abs(tkelosses_2nd), 'bo-', label='2ndx')
loglog(dts, abs(tkelosses_4th), 'go-', label='4thx')
loglog(dts, slope2, 'k--', label='2nd')
loglog(dts, slope3, 'k:' , label='3rd')
loglog(dts, slope4, 'k-.', label='4th')
xlabel('dt')
ylabel('tke loss')
legend(loc=0, frameon=False)
