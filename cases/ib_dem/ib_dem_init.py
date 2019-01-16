import numpy as np
import matplotlib.pylab as pl

from microhh_tools import *    # available in MICROHH_DIR/python; copy to this directory

# Class to create a smooth stretched vertical grid
class Grid:
    def __init__(self, kmax, nloc1, nbuf1, dz1, dz2):
        dn         = 1./kmax
        n          = np.linspace(dn, 1.-dn, kmax)
        nloc1     *= dn
        nbuf1     *= dn
        dzdn1      = dz1/dn
        dzdn2      = dz2/dn

        dzdn       = dzdn1 + 0.5*(dzdn2-dzdn1)*(1. + np.tanh((n-nloc1)/nbuf1))
        self.dz    = dzdn*dn

        self.kmax  = kmax
        self.z     = np.zeros(self.dz.size)
        stretch    = np.zeros(self.dz.size)

        self.z[0]  = 0.5*self.dz[0]
        stretch[0] = 1.

        for k in range(1, self.kmax):
              self.z [k] = self.z[k-1] + 0.5*(self.dz[k-1]+self.dz[k])
              stretch[k] = self.dz[k]/self.dz[k-1]

        self.zsize = self.z[kmax-1] + 0.5*self.dz[kmax-1]
        print('kmax={0:}, zsize={1:.4f}'.format(kmax,self.zsize))

    def plot(self):
        pl.figure()
        pl.plot(self.z, self.dz, '-x')
        pl.xlabel('z (m)')
        pl.ylabel('dz (m)')

def z_spectral(itot, jtot, k1, fac1, k2, fac2, ks2, dz, z_offset):
    x = np.linspace(0.+ 0.5//itot, 1.-0.5//itot, itot)
    y = np.linspace(0.+ 0.5//jtot, 1.-0.5//jtot, jtot)

    xfft = np.linspace(0, itot//2, itot//2+1)
    yfft = np.linspace(0, jtot//2, jtot//2+1)

    np.random.seed(103)
    arnd = 2.*np.pi*np.random.rand(jtot, itot//2+1)
    afftrnd = np.cos(arnd) + 1j*np.sin(arnd)

    # Calculate the radial wave numbers
    l = np.zeros(afftrnd.shape)
    for i in range(0,itot//2+1):
        for j in range(0,jtot//2+1):
            l[i,j] = (i**2. + j**2.)**.5
    for i in range(itot//2+1,itot):
        for j in range(0,jtot//2+1):
            l[i,j] = ((itot-i)**2. + j**2.)**.5

    # Filter on radial wave number using a gaussian function
    factor  = fac1*np.exp(-(l-k1)**2. / (2.*ks2))
    factor += fac2*np.exp(-(l-k2)**2. / (2.*ks2))

    # Create the filtered field
    afft = factor*afftrnd

    # Make sure the mean is exactly 0
    afft[0,0] = 0.

    a = np.fft.irfft2(afft)

    # Add ofset to set lower boundary at z=0
    a -= a.min()

    # Scale to requested max height
    a *= dz / a.max()

    # Add offset
    a += z_offset

    return a


if __name__ == "__main__":
    # Read the .ini file
    ini = Read_namelist()

    # Create stretched grid
    grid = Grid(96, 40, 5, 0.0004, 0.0007)
    #grid.plot()

    # Write `zsize` and `ktot` back to .ini file
    replace_namelist_value('zsize', grid.zsize)
    replace_namelist_value('ktot',  grid.kmax)

    # Create initial profiles:
    z = grid.z
    u = ini['force']['uflux'] * np.ones(z.size)
    s = z

    # Write the data to a .prof file for MicroHH
    proffile = open('ib_dem.prof','w')
    proffile.write('{0:^20s} {1:^20s} {2:^20s}\n'.format('z', 'u', 's'))
    for k in range(z.size):
        proffile.write('{0:1.14E} {1:1.14E} {2:1.14E}\n'.format(z[k], u[k], s[k]))
    proffile.close()

    # Create 2D IB height
    zIB = z_spectral(ini['grid']['itot'], ini['grid']['jtot'], k1=2, fac1=4, k2=6, fac2=1, ks2=1.5**2., dz=0.01, z_offset=grid.z[2])

    dx = ini['grid']['xsize'] / ini['grid']['itot']
    dy = ini['grid']['ysize'] / ini['grid']['jtot']
    x  = np.arange(0.5*dx, ini['grid']['xsize'], dx)
    y  = np.arange(0.5*dy, ini['grid']['ysize'], dy)

    #pl.figure()
    #pl.pcolormesh(x,y,zIB)
    #pl.colorbar()

    # Write to input file for MicroHH
    write_restart_file('dem.0000000', zIB[np.newaxis,:,:], ini['grid']['itot'], ini['grid']['jtot'], 1)

    #sbot_ib = zIB - zIB.min()
    #sbot_ib /= sbot_ib.max()
    #sbot_ib = 1. - sbot_ib

    tanh_x = np.tanh( (x - ini['grid']['xsize']/2.) / 0.01 )
    sbot_ib = np.ones(zIB.shape) * tanh_x[None,:]

    # pl.figure()
    # pl.pcolormesh(x,y,sbot_ib)
    # pl.colorbar()
    # pl.show()

    write_restart_file('sbot_ib.0000000', sbot_ib[np.newaxis,:,:], ini['grid']['itot'], ini['grid']['jtot'], 1)

