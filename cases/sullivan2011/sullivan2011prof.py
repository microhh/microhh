import numpy

# set the height
kmax  = 64
zsize = 2048.
dz = zsize / kmax

# set the height
z    = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
s    = numpy.zeros(numpy.size(z))

for k in range(kmax):
  # temperature
  if(z[k] <= 974.):
    s[k] = 300.
  elif(z[k] <= 1074):
    s[k] = 300. + (z[k]-974.)*0.08
  else:
    s[k] = 308. + (z[k]-1074.)*0.003

# write the data to a file
proffile = open('sullivan2011.prof','w')
proffile.write('{0:^20s} {1:^20s}\n'.format('z','s'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E}\n'.format(z[k], s[k]))
proffile.close()

