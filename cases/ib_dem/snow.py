import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
import numpy as np

# fraction of landscape not snow covered
cover = 0.8

# resolution
itot = 64
jtot = 64

# the mean wavenumber of perturbation 1
k1 = 2.
fac1 = 3.

# the mean wavenumber of perturbation 2
k2 = 3*k1
fac2 = 1.

# the variance in wave number of perturbation
ks2 = 1.5**2.

x = np.linspace(0.+ 0.5//itot, 1.-0.5//itot, itot)
y = np.linspace(0.+ 0.5//jtot, 1.-0.5//jtot, jtot)

xfft = np.linspace(0, itot//2, itot//2+1)
yfft = np.linspace(0, jtot//2, jtot//2+1)

np.random.seed(103)
arnd = 2.*np.pi*np.random.rand(jtot, itot//2+1)
afftrnd = np.cos(arnd) + 1j*np.sin(arnd)
# calculate the radial wave numbers
l = np.zeros(afftrnd.shape)
for i in range(0,itot//2+1):
    for j in range(0,jtot//2+1):
        l[i,j] = (i**2. + j**2.)**.5
for i in range(itot//2+1,itot):
    for j in range(0,jtot//2+1):
        l[i,j] = ((itot-i)**2. + j**2.)**.5

# filter on radial wave number using a gaussian function
#fftfilter = zeros(sfft.shape, dtype=np.complex)
factor  = fac1*np.exp(-(l-k1)**2. / (2.*ks2))
factor += fac2*np.exp(-(l-k2)**2. / (2.*ks2))

# create a line for plotting with the spectral bands
factork = np.linspace(0., 25., 1000)
factory1 = fac1*np.exp(-(factork-k1)**2. / (2.*ks2))
factory2 = fac2*np.exp(-(factork-k2)**2. / (2.*ks2))

# create the filtered field
afft = factor*afftrnd

# make sure the mean is exactly 0
afft[0,0] = 0.

a = np.fft.irfft2(afft)

# normalize the variance to 1
a /= np.std(a)

asort = np.sort(a, axis=None)
index = asort.size-1 - int(cover*asort.size)
if(index == -1):
  athres = asort[0]-1
else:
  athres = asort[index]

amask = a.copy()

a /= 6
a -= 1.1*a.min()

crange = np.linspace(0.,1.,2)
fig = pl.figure(figsize=(6,6.6))
fig.subplots_adjust(0.11,0.07,0.97,0.97,0.1,0.4)
gs = gridspec.GridSpec(4, 1)
sub1 = pl.subplot(gs[0:3,0])
pl.imshow(amask, extent=(0,500,0,500), origin='lower', vmin=athres, vmax=a.max())
pl.colorbar()
pl.contour(500*x, 500*y, amask, colors='k')
pl.xlabel(r'$x$ (m)')
pl.ylabel(r'$y$ (m)')
pl.title('a) Snow cover and height lines', loc='left')
sub2 = pl.subplot(gs[3,0])
pl.plot(factork, (1./3.)*(factory1+factory2), 'r-', alpha=0.4, linewidth=4)
pl.plot(factork, (1./3.)*factory1, 'k--')
pl.plot(factork, (1./3.)*factory2, 'k:')
pl.xlabel(r'$k$')
pl.ylabel(r'$E(k)$')
pl.ylim(0., 1.1)
pl.title('b) Normalized energy spectrum of orography', loc='left')
pl.tight_layout()
