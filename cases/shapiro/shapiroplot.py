import numpy as np
import struct as st
import pylab as pl

class Data:
    def __init__(self, path, nx, nz, N2):
        ny = 1

        n   = nx*nz
        fin = open("{0}/grid.{1:07d}".format(path, 0),"rb")
        raw = fin.read(nx*8)
        self.x = np.array(st.unpack('{0}d'.format(nx), raw))
        raw = fin.read(nx*8)
        self.xh = np.array(st.unpack('{0}d'.format(nx), raw))
        raw = fin.read(ny*8)
        self.y = np.array(st.unpack('{0}d'.format(ny), raw))
        raw = fin.read(ny*8)
        self.yh = np.array(st.unpack('{0}d'.format(ny), raw))
        raw = fin.read(nz*8)
        self.z = np.array(st.unpack('{0}d'.format(nz), raw))
        raw = fin.read(nz*8)
        self.zh = np.array(st.unpack('{0}d'.format(nz), raw))
        fin.close()
        
        fin = open("{0}/u.0010000".format(path), "rb")
        raw = fin.read(nx*nz*8)
        tmp = np.array(st.unpack('{0}d'.format(nx*nz), raw))
        del(raw)
        tmp = tmp.reshape((nz, nx))
        self.u = np.roll(tmp, nx/4, 1)
        del(tmp)
         
        fin = open("{0}/w.0010000".format(path), "rb")
        raw = fin.read(nx*nz*8)
        tmp = np.array(st.unpack('{0}d'.format(nx*nz), raw))
        del(raw)
        tmp = tmp.reshape((nz, nx))
        self.w = np.roll(tmp, nx/4, 1)
        del(tmp)
         
        fin = open("{0}/b.0010000".format(path), "rb")
        raw = fin.read(nx*nz*8)
        tmp = np.array(st.unpack('{0}d'.format(nx*nz), raw))
        del(raw)
        tmp = tmp.reshape((nz, nx))
        self.b = np.roll(tmp, nx/4, 1)
        del(tmp)
        for k in range(nz):
            self.b[k,:] -= N2 * self.z[k]

        # Load reference data
        nx_hires = 513
        nz_hires = 1025
        fin = open("u.out", "rb")
        raw = fin.read(nx_hires*nz_hires*8)
        tmp = np.array(st.unpack('{0}d'.format(nx_hires*nz_hires), raw))
        del(raw)
        u_hires = tmp.reshape((nz_hires, nx_hires))
        del(tmp)
         
        fin = open("w.out", "rb")
        raw = fin.read(nx_hires*nz_hires*8)
        tmp = np.array(st.unpack('{0}d'.format(nx_hires*nz_hires), raw))
        del(raw)
        w_hires = tmp.reshape((nz_hires, nx_hires))
        del(tmp)
         
        fin = open("b.out", "rb")
        raw = fin.read(nx_hires*nz_hires*8)
        tmp = np.array(st.unpack('{0}d'.format(nx_hires*nz_hires), raw))
        del(raw)
        b_hires = tmp.reshape((nz_hires, nx_hires))
        del(tmp)
        
        # Get data for 64 x 128 case
        stepx   = nx_hires / nx
        startxh = 0
        startx  = stepx / 2
        stepz   = nz_hires / nz
        startz  = stepz / 2
        startzh = 0
        self.u_ref = u_hires[startz :nz_hires-1:stepz, startxh:nx_hires-1:stepx]
        self.w_ref = w_hires[startzh:nz_hires-1:stepz, startx :nx_hires-1:stepx]
        self.b_ref = b_hires[startz :nz_hires-1:stepz, startx :nx_hires-1:stepx]

        self.u_error = 0.
        self.w_error = 0.
        self.b_error = 0.

        self.dx = self.x[1] - self.x[0]
        for k in range(self.z.size):
            self.u_error += np.sqrt(self.dx*sum((self.u[k,:] - self.u_ref[k,:])**2))
            self.w_error += np.sqrt(self.dx*sum((self.w[k,:] - self.w_ref[k,:])**2))
            self.b_error += np.sqrt(self.dx*sum((self.b[k,:] - self.b_ref[k,:])**2))

        self.u_error = self.u_error / self.z.size
        self.w_error = self.w_error / self.z.size
        self.b_error = self.b_error / self.z.size

N2 = 0.0004
data32 = Data("32x64" , 32,  64, N2)
data64 = Data("64x128", 64, 128, N2)

sl = [data32, data64]

print('errors 32 x  64: ', data32.u_error, data32.w_error, data32.b_error)
print('errors 64 x 128: ', data64.u_error, data64.w_error, data64.b_error)
print('convergence u', (np.log(data64.u_error)-np.log(data32.u_error)) / (np.log(data64.dx)-np.log(data32.dx)) )
print('convergence w', (np.log(data64.w_error)-np.log(data32.w_error)) / (np.log(data64.dx)-np.log(data32.dx)) )
print('convergence p', (np.log(data64.b_error)-np.log(data32.b_error)) / (np.log(data64.dx)-np.log(data32.dx)) )

# Visualisation of fields.
pl.close('all')
for s in sl:
    pl.figure(figsize=(10,6))
    pl.subplot(131)
    pl.pcolormesh(s.xh, s.z, s.u - s.u_ref)
    pl.xlim(0, 5.12)
    pl.ylim(0, 5.12)
    pl.title('u_error')
    pl.colorbar()
    pl.subplot(132)
    pl.pcolormesh(s.x, s.zh, s.w - s.w_ref)
    pl.xlim(0, 5.12)
    pl.ylim(0, 5.12)
    pl.title('w_error')
    pl.colorbar()
    pl.subplot(133)
    pl.pcolormesh(s.x, s.z, s.b - s.b_ref)
    pl.xlim(0, 5.12)
    pl.ylim(0, 5.12)
    pl.title('b_error')
    pl.colorbar()
