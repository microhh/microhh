import numpy as np
import matplotlib.pyplot as plt

# Set grid cells.
itot = 128
jtot = 128
ktot = 64

# Set standard deviation random noise
th_std = 1e-2

# Set properties of gaussian filter in wave number space
K_mean = 2.
K_std = 1.

# DO NOT EDIT BELOW.
phase = 2.*np.pi*np.random.rand(ktot, jtot, itot//2+1)
th_fft = np.cos(phase) + 1j*np.sin(phase)

# Remove the mean.
th_fft[0,0,0] = 0;

# Calculate the radial wave numbers.
K = np.zeros(th_fft.shape)

kk = np.arange(itot//2+1)
ll = np.zeros(jtot)
ll[:jtot//2+1] = np.arange(jtot//2+1)
for j in range(jtot//2+1,jtot):
    ll[j] = jtot-j
mm = np.zeros(ktot)
mm[:ktot//2+1] = np.arange(ktot//2+1)
for k in range(ktot//2+1,ktot):
    mm[k] = ktot-k

K[:,:,:] = (kk[None,None,:]**2 + ll[None,:,None]**2 + mm[:,None,None]**2)**.5

# Filter on radial wave number using a gaussian function.
factor = np.exp(-(K-K_mean)**2. / (2.*K_std**2))
th_fft *= factor

th = np.fft.irfftn(th_fft)
th *= th_std/th.std()

th_file = open("th.0000000", "rb")
th_ini = np.fromfile(th_file, dtype=np.float64)
th_file.close()

th_ini += th.flatten()
th_ini.tofile("th.0000000")

# Create some test plots.
plt.figure()
plt.subplot(221)
plt.pcolormesh(th[:,:,20])
plt.colorbar()
plt.subplot(222)
plt.pcolormesh(th[:,20,:])
plt.colorbar()
plt.subplot(223)
plt.pcolormesh(th[20,:,:])
plt.colorbar()
plt.subplot(224)
plt.hist(th.flatten(), 30, density=True)
plt.tight_layout()
plt.show()
