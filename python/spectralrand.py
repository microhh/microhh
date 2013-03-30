from pylab import *

itot = 512
jtot = 512

# the mean wavenumber of the perturbation
k = 20.

# the variance in wave number of the perturbation
ks2 = 2.**2.

# create random field
#sraw = rand(itot, jtot)

# filter the data
#stemp = zeros((itot+6, jtot+6))
#stemp[3:itot+3, 3:jtot+3] = sraw[:,:]
#stemp[0:3,:] = stemp[itot:itot+3,:]
#stemp[itot+3:itot+6,:] = stemp[3:6,:]
#stemp[:,0:3] = stemp[:,jtot:jtot+3]
#stemp[:,jtot+3:jtot+6] = stemp[:,3:6]
#
#kernel = zeros((7,7))
#for ii in range(7):
#  for jj in range(7):
#    kernel[ii,jj] = exp(-(ii-3)**2./2.)*exp(-(jj-3)**2./2.)
#
#s = zeros((itot,jtot))
#for i in range(itot):
#  for j in range(jtot):
#    s[i,j] = sum(stemp[i:i+7,j:j+7]*kernel[:,:])

# create random field
s = rand(itot, jtot)
sfft = rfft2(s)

# set mean to zero
sfft[0,0] = 0.

# calculate the radial wave numbers
l = zeros(sfft.shape)
for i in range(0,itot/2+1):
  for j in range(0,jtot/2+1):
    l[i,j] = (i**2. + j**2.)**.5
for i in range(itot/2+1,itot):
  for j in range(0,jtot/2+1):
    l[i,j] = ((itot-i)**2. + j**2.)**.5

# filter on radial wave number using a gaussian function
sfftfilter = zeros((sfft.shape[0], sfft.shape[1]), dtype=np.complex)
factor = exp(-(l-k)**2. / (2.*ks2))
sfftfilter = factor*sfft

sfilter = irfft2(sfftfilter)

# normalize the variance to 1.
sfilter /= var(sfilter)**0.5

close('all')
figure()
pcolormesh(sfilter.transpose())
colorbar();
xlim(0,itot)
ylim(0,jtot)
xlabel('x')
ylabel('y')
title('random field with spectral properties and var = 1')

figure()
pcolormesh(abs(sfftfilter.transpose())[0:itot/2+1,:])
colorbar();
xlim(0,itot/2+1)
ylim(0,jtot/2+1)
xlabel('k')
ylabel('l')
title('shape of the spectral filter for positive wave numbers')

