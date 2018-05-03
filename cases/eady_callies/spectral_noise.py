import numpy as np
import matplotlib.pyplot as plt

# Set grid cells.
itot = 128
jtot = 128
ktot = 64

# Set standard deviation random noise
b_std = 1e-7

# Set properties of gaussian filter in wave number space
K_mean = 2.
K_std = 1.

# DO NOT EDIT BELOW.
phase = 2.*np.pi*np.random.rand(ktot, jtot, itot//2+1)
b_fft = np.cos(phase) + 1j*np.sin(phase)

# Remove the mean.
b_fft[0,0,0] = 0;

# Calculate the radial wave numbers.
K = np.zeros(b_fft.shape)

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
b_fft *= factor

b = np.fft.irfftn(b_fft)
b *= b_std/b.std()

b_file = open("b.0000000", "rb")
b_ini = np.fromfile(b_file, dtype=np.float64)
b_file.close()

b_ini += b.flatten()
b_ini.tofile("b.0000000")

# Create some test plots.
plt.figure()
plt.subplot(221)
plt.pcolormesh(b[:,:,20])
plt.colorbar()
plt.subplot(222)
plt.pcolormesh(b[:,20,:])
plt.colorbar()
plt.subplot(223)
plt.pcolormesh(b[20,:,:])
plt.colorbar()
plt.subplot(224)
plt.hist(b.flatten(), 30, density=True)
plt.tight_layout()
plt.show()
