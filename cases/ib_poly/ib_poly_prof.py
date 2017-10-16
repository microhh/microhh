import numpy as np
import matplotlib.pylab as pl
from scipy.special import erf

# Create stretched grid
class Stretched_grid:
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
        print('Grid: kmax=%i, zsize=%f'%(kmax,self.zsize))
        
    def plot(self):
        pl.figure()
        pl.plot(self.dz, self.z, '-x')
        pl.xlabel('dz (m)')
        pl.ylabel('z (m)')      


if (__name__ == '__main__'):
    import microhh_tools as mht

    # Create stretched grid
    gr = Stretched_grid(64, 40, 10, 2.0, 10.0)
    #gr.plot()
    
    # Write back vertical extent domain
    mht.replace_namelist_value('zsize', gr.zsize)
    mht.replace_namelist_value('ktot', gr.z.size)
    
    # Create initial vertical profiles
    u = 0.09/0.4 * np.log(gr.z)
    v = np.zeros_like(gr.z)
    
    # Write to .prof file as input for MicroHH
    variables = {'z':gr.z, 'u':u, 'v':v}
    mht.write_output('ib_poly.prof', variables, gr.z.size)
