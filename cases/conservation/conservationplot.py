from pylab import *

data100_2nd = loadtxt('conservation100_2nd/conservation.out', skiprows=1)
data200_2nd = loadtxt('conservation200_2nd/conservation.out', skiprows=1)
data400_2nd = loadtxt('conservation400_2nd/conservation.out', skiprows=1)

data100_4th = loadtxt('conservation100_4th/conservation.out', skiprows=1)
data200_4th = loadtxt('conservation200_4th/conservation.out', skiprows=1)
data400_4th = loadtxt('conservation400_4th/conservation.out', skiprows=1)

time100_2nd = data100_2nd[:,1]
mom100_2nd  = data100_2nd[:,7]
tke100_2nd  = data100_2nd[:,8]
mass100_2nd = data100_2nd[:,9]

time200_2nd = data200_2nd[:,1]
mom200_2nd  = data200_2nd[:,7]
tke200_2nd  = data200_2nd[:,8]
mass200_2nd = data200_2nd[:,9]

time400_2nd = data400_2nd[:,1]
mom400_2nd  = data400_2nd[:,7]
tke400_2nd  = data400_2nd[:,8]
mass400_2nd = data400_2nd[:,9]

time100_4th = data100_4th[:,1]
mom100_4th  = data100_4th[:,7]
tke100_4th  = data100_4th[:,8]
mass100_4th = data100_4th[:,9]

time200_4th = data200_4th[:,1]
mom200_4th  = data200_4th[:,7]
tke200_4th  = data200_4th[:,8]
mass200_4th = data200_4th[:,9]

time400_4th = data400_4th[:,1]
mom400_4th  = data400_4th[:,7]
tke400_4th  = data400_4th[:,8]
mass400_4th = data400_4th[:,9]

figure()
subplot(131)
plot(time100_2nd[1:200], mom100_2nd [1:200] - mom100_2nd [1], 'b-' , label = 'momdiff100 ')
plot(time200_2nd[1:400], mom200_2nd [1:400] - mom200_2nd [1], 'b--', label = 'momdiff200 ')
plot(time400_2nd[1:800], mom400_2nd [1:800] - mom400_2nd [1], 'b:' , label = 'momdiff400 ')
legend(loc=0, frameon=False)
subplot(132)
plot(time100_2nd[1:200], tke100_2nd [1:200] - tke100_2nd [1], 'b-' , label = 'tkediff100 ')
plot(time200_2nd[1:400], tke200_2nd [1:400] - tke200_2nd [1], 'b--', label = 'tkediff200 ')
plot(time400_2nd[1:800], tke400_2nd [1:800] - tke400_2nd [1], 'b:' , label = 'tkediff400 ')
legend(loc=0, frameon=False)
subplot(133)
plot(time100_2nd[1:200], mass100_2nd[1:200] - mass100_2nd[1], 'b-' , label = 'massdiff100')
plot(time200_2nd[1:400], mass200_2nd[1:400] - mass200_2nd[1], 'b--', label = 'massdiff200')
plot(time400_2nd[1:800], mass400_2nd[1:800] - mass400_2nd[1], 'b:' , label = 'massdiff400')
legend(loc=0, frameon=False)

figure()
subplot(131)
plot(time100_4th[1:200], mom100_4th [1:200] - mom100_4th [1], 'r-' , label = 'momdiff100 ')
plot(time200_4th[1:400], mom200_4th [1:400] - mom200_4th [1], 'r--', label = 'momdiff200 ')
plot(time400_4th[1:800], mom400_4th [1:800] - mom400_4th [1], 'r:' , label = 'momdiff400 ')
legend(loc=0, frameon=False)
subplot(132)
plot(time100_4th[1:200], tke100_4th [1:200] - tke100_4th [1], 'r-' , label = 'tkediff100 ')
plot(time200_4th[1:400], tke200_4th [1:400] - tke200_4th [1], 'r--', label = 'tkediff200 ')
plot(time400_4th[1:800], tke400_4th [1:800] - tke400_4th [1], 'r:' , label = 'tkediff400 ')
legend(loc=0, frameon=False)
subplot(133)
plot(time100_4th[1:200], mass100_4th[1:200] - mass100_4th[1], 'r-' , label = 'massdiff100')
plot(time200_4th[1:400], mass200_4th[1:400] - mass200_4th[1], 'r--', label = 'massdiff200')
plot(time400_4th[1:800], mass400_4th[1:800] - mass400_4th[1], 'r:' , label = 'massdiff400')
legend(loc=0, frameon=False)

figure()
subplot(131)
plot(time100_2nd[1:200], mom100_2nd [1:200] - mom100_2nd [1], 'b-' , label = 'momdiff100_2nd' )
plot(time100_4th[1:200], mom100_4th [1:200] - mom100_4th [1], 'r-' , label = 'momdiff100_4th' )
plot(time200_2nd[1:400], mom200_2nd [1:400] - mom200_2nd [1], 'b--', label = 'momdiff200_2nd' )
plot(time200_4th[1:400], mom200_4th [1:400] - mom200_4th [1], 'r--', label = 'momdiff200_4th' )
plot(time400_2nd[1:800], mom400_2nd [1:800] - mom400_2nd [1], 'b:' , label = 'momdiff400_2nd' )
plot(time400_4th[1:800], mom400_4th [1:800] - mom400_4th [1], 'r:' , label = 'momdiff400_4th' )
legend(loc=0, frameon=False)
subplot(132)
plot(time100_2nd[1:200], tke100_2nd [1:200] - tke100_2nd [1], 'b-' , label = 'tkediff100_2nd' )
plot(time100_4th[1:200], tke100_4th [1:200] - tke100_4th [1], 'r-' , label = 'tkediff100_4th' )
plot(time200_2nd[1:400], tke200_2nd [1:400] - tke200_2nd [1], 'b--', label = 'tkediff200_2nd' )
plot(time200_4th[1:400], tke200_4th [1:400] - tke200_4th [1], 'r--', label = 'tkediff200_4th' )
plot(time400_2nd[1:800], tke400_2nd [1:800] - tke400_2nd [1], 'b:' , label = 'tkediff400_2nd' )
plot(time400_4th[1:800], tke400_4th [1:800] - tke400_4th [1], 'r:' , label = 'tkediff400_4th' )
legend(loc=0, frameon=False)
subplot(133)
plot(time100_2nd[1:200], mass100_2nd[1:200] - mass100_2nd[1], 'b-' , label = 'massdiff100_2nd')
plot(time100_4th[1:200], mass100_4th[1:200] - mass100_4th[1], 'r-' , label = 'massdiff100_4th')
plot(time200_2nd[1:400], mass200_2nd[1:400] - mass200_2nd[1], 'b--', label = 'massdiff200_2nd')
plot(time200_4th[1:400], mass200_4th[1:400] - mass200_4th[1], 'r--', label = 'massdiff200_4th')
plot(time400_2nd[1:800], mass400_2nd[1:800] - mass400_2nd[1], 'b:' , label = 'massdiff400_2nd')
plot(time400_4th[1:800], mass400_4th[1:800] - mass400_4th[1], 'r:' , label = 'massdiff400_4th')
legend(loc=0, frameon=False)

"""
timesteps = array([100, 200, 400])
momerror  = array([mom100[-1]-mom100[1], mom200[-1]-mom200[1], mom400[-1]-mom400[1]])
tkeerror  = array([tke100[-1]-tke100[1], tke200[-1]-tke200[1], tke400[-1]-tke400[1]])
masserror = array([mass100[-1]-mass100[1], mass200[-1]-mass200[1], mass400[-1]-mass400[1]])

figure()
subplot(131)
loglog(timesteps, abs(momerror))
subplot(132)
loglog(timesteps, abs(tkeerror))
subplot(133)
loglog(timesteps, abs(masserror))
"""
