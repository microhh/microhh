from pylab import *

data100 = loadtxt('conservation100/conservation.out', skiprows=1)
data200 = loadtxt('conservation200/conservation.out', skiprows=1)
data400 = loadtxt('conservation400/conservation.out', skiprows=1)

time100 = data100[:,1]
mom100  = data100[:,7]
tke100  = data100[:,8]
mass100 = data100[:,9]

time200 = data200[:,1]
mom200  = data200[:,7]
tke200  = data200[:,8]
mass200 = data200[:,9]

time400 = data400[:,1]
mom400  = data400[:,7]
tke400  = data400[:,8]
mass400 = data400[:,9]

figure()
subplot(131)
plot(time100[1:100], mom100 [1:100] - mom100 [1], 'b-' , label = 'momdiff100 ')
plot(time200[1:200], mom200 [1:200] - mom200 [1], 'b--', label = 'momdiff200 ')
plot(time400[1:400], mom400 [1:400] - mom400 [1], 'b:' , label = 'momdiff400 ')
legend(loc=0, frameon=False)
subplot(132)
plot(time100[1:100], tke100 [1:100] - tke100 [1], 'g-' , label = 'tkediff100 ')
plot(time200[1:200], tke200 [1:200] - tke200 [1], 'g--', label = 'tkediff200 ')
plot(time400[1:400], tke400 [1:400] - tke400 [1], 'g:' , label = 'tkediff400 ')
legend(loc=0, frameon=False)
subplot(133)
plot(time100[1:100], mass100[1:100] - mass100[1], 'r-' , label = 'massdiff100')
plot(time200[1:200], mass200[1:200] - mass200[1], 'r--', label = 'massdiff200')
plot(time400[1:400], mass400[1:400] - mass400[1], 'r:' , label = 'massdiff400')
legend(loc=0, frameon=False)

timesteps = array([100, 200, 400])
momerror  = array([mom100[-1]-mom100[1], mom200[-1]-mom200[1], mom400[-1]-mom400[1]])
tkeerror  = array([tke100[-1]-tke100[1], tke200[-1]-tke200[1], tke400[-1]-tke400[1]])
masserror = array([mass100[-1]-mass100[1], mass200[-1]-mass200[1], mass400[-1]-mass400[1]])

"""
figure()
subplot(131)
loglog(timesteps, abs(momerror))
subplot(132)
loglog(timesteps, abs(tkeerror))
subplot(133)
loglog(timesteps, abs(masserror))
"""
