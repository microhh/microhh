import numpy as np
from pylab import *

# Read table A1 paper
A1 = np.loadtxt('andren1994_tableA1')
z  = A1[:,0]
u  = A1[:,1]
v  = A1[:,2]
q2 = A1[:,3]

ug = np.zeros(z.size)
vg = np.zeros(z.size)

ug[:] = 10.

proffile = open('andren1994.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s} \n'.format('z','u','v','ug','vg'))
for k in range(z.size):
  proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} \n'.format(z[k], u[k], v[k], ug[k], vg[k]))
proffile.close()
