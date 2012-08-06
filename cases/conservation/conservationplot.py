from pylab import *

data = loadtxt('conservation.out', skiprows=1)

nt   = 100
time = data[:,1]
mom  = data[:,7]
tke  = data[:,8]
mass = data[:,9]

figure()
plot(time[1:nt], mom [1:nt] - mom [1], label = 'momdiff ')
plot(time[1:nt], tke [1:nt] - tke [1], label = 'tkediff ')
plot(time[1:nt], mass[1:nt] - mass[1], label = 'massdiff')
legend(loc=0, frameon=False)


