from pylab import *

ss1024sm = loadtxt("strongscaling.1024.supermuc", skiprows=1);
ss1024sm_procs      = ss1024sm[:,0]
ss1024sm_times      = ss1024sm[:,1]
ss1024sm_speedup    = ss1024sm_times[0] / ss1024sm_times[:]
ss1024sm_linspeedup = ss1024sm_procs[:] / ss1024sm_procs[0]
ss1024sm_eff        = ss1024sm_speedup / ss1024sm_linspeedup

ss2048sm = loadtxt("strongscaling.2048.supermuc", skiprows=1);
ss2048sm_procs      = ss2048sm[:,0]
ss2048sm_times      = ss2048sm[:,1]
ss2048sm_speedup    = ss2048sm_times[0] / ss2048sm_times[:]
ss2048sm_linspeedup = ss2048sm_procs[:] / ss2048sm_procs[0]
ss2048sm_eff        = ss2048sm_speedup / ss2048sm_linspeedup

procsx      = array([256, 512, 1024, 2048, 4096, 8192, 16384])
linspeedupy = array([0.5, 1, 2, 4, 8, 16, 32])
effy        = ones(7)

close('all')
figure()
loglog(ss1024sm_procs, ss1024sm_speedup, 'bo-', label='1024x1024x1024')
loglog(ss2048sm_procs, ss2048sm_speedup, 'ro-', label='2048x2048x1024')
loglog(procsx, linspeedupy, 'k-' )
xlim(procsx     [0], procsx     [-1])
ylim(linspeedupy[0], linspeedupy[-1])
xlabel('#processes')
ylabel('speedup')
legend(loc=0, frameon=False)
grid()
title('strong scaling SuperMUC')

figure()
semilogx(ss1024sm_procs, ss1024sm_eff, 'bo-', label='1024x1024x1024')
semilogx(ss2048sm_procs, ss2048sm_eff, 'ro-', label='2048x2048x1024')
semilogx(procsx, effy, 'k-' )
xlim(procsx[0], procsx[-1])
ylim(0., 1.5)
xlabel('#processes')
ylabel('efficiency')
legend(loc=0, frameon=False)
grid()
title('strong scaling SuperMUC')

wssmallsm = loadtxt("weakscaling.small.supermuc", skiprows=1);
wssmallsm_procs      = wssmallsm[:,0]
wssmallsm_times      = wssmallsm[:,1]
wssmallsm_eff        = wssmallsm_times[0] / wssmallsm_times[:]

wslargesm = loadtxt("weakscaling.large.supermuc", skiprows=1);
wslargesm_procs      = wslargesm[:,0]
wslargesm_times      = wslargesm[:,1]
wslargesm_eff        = wslargesm_times[0] / wslargesm_times[:]

figure()
semilogx(wssmallsm_procs, wssmallsm_eff, 'bo-', label='32x16x1024 per proc')
semilogx(wslargesm_procs, wslargesm_eff, 'ro-', label='32x32x1024 per proc')
semilogx(procsx, effy, 'k-' )
xlim(procsx[0], procsx[-1])
ylim(0., 1.5)
xlabel('#processes')
ylabel('efficiency')
legend(loc=0, frameon=False)
grid()
title('weak scaling SuperMUC')

