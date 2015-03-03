import numpy

# set the height
# Get number of vertical levels and size from .ini file
with open('sullivan2011.ini') as f:
  for line in f:
    if(line.split('=')[0]=='ktot'):
      kmax = int(line.split('=')[1])
    if(line.split('=')[0]=='zsize'):
      zsize = float(line.split('=')[1])

dz = zsize / kmax

# set the height
z    = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
th   = numpy.zeros(numpy.size(z))

for k in range(kmax):
  # temperature
  if(z[k] <= 974.):
    th[k] = 300.
  elif(z[k] <= 1074):
    th[k] = 300. + (z[k]-974.)*0.08
  else:
    th[k] = 308. + (z[k]-1074.)*0.003

# write the data to a file
proffile = open('sullivan2011.prof','w')
proffile.write('{0:^20s} {1:^20s}\n'.format('z','th'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E}\n'.format(z[k], th[k]))
proffile.close()

