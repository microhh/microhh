import struct
from numpy import *

class microhh:
  def __init__(self, iter, itot, ktot, path):
    nx = itot
    ny = 1
    nz = ktot
    
    n  = nx*ny*nz
    
    fin = open("{0:s}/grid.{1:07d}".format(path, 0),"rb")
    raw = fin.read(nx*8)
    self.x = array(struct.unpack('<{}d'.format(nx), raw))
    raw = fin.read(nx*8)
    self.xh = array(struct.unpack('<{}d'.format(nx), raw))
    raw = fin.read(ny*8)
    self.y = array(struct.unpack('<{}d'.format(ny), raw))
    raw = fin.read(ny*8)
    self.yh = array(struct.unpack('<{}d'.format(ny), raw))
    raw = fin.read(nz*8)
    self.z = array(struct.unpack('<{}d'.format(nz), raw))
    raw = fin.read(nz*8)
    self.zh = array(struct.unpack('<{}d'.format(nz), raw))
    fin.close()
    
    fin = open("{0:s}/u.xz.00000.{1:07d}".format(path, iter),"rb")
    raw = fin.read(n*8)
    tmp = array(struct.unpack('<{}d'.format(n), raw))
    del(raw)
    self.u = tmp.reshape((nz, ny, nx))
    del(tmp)
    fin.close()
    
    fin = open("{0:s}/w.xz.00000.{1:07d}".format(path, iter),"rb")
    raw = fin.read(n*8)
    tmp = array(struct.unpack('<{}d'.format(n), raw))
    del(raw)
    self.w = tmp.reshape((nz, ny, nx))
    del(tmp)
    fin.close()
    
    fin = open("{0:s}/p.xz.00000.{1:07d}".format(path, iter),"rb")
    raw = fin.read(n*8)
    tmp = array(struct.unpack('<{}d'.format(n), raw))
    del(raw)
    self.p = tmp.reshape((nz, ny, nx))
    del(tmp)
    fin.close()

class getref:
  def __init__(self, x, xh, z, zh, visc, time):
    self.u = zeros((zh.size, 1, x .size))
    self.w = zeros((z .size, 1, xh.size))
    self.p = zeros((z .size, 1, x .size))
    
    for k in range(z.size):
      self.u[k,0,:] =  sin(2.*pi*xh)*cos(2.*pi*z [k])*exp(-8.*pi**2.*visc*time)
      self.w[k,0,:] = -cos(2.*pi*x )*sin(2.*pi*zh[k])*exp(-8.*pi**2.*visc*time)
      self.p[k,0,:] = (0.25*(cos(4.*pi*x) + cos(4.*pi*z[k]))-0.25)*(exp(-8.*pi**2.*visc*time)**2.)

class geterror:
  def __init__(self, data, ref):
    dx = 1.  / data.x.size
    dz = 0.5 / data.z.size

    zh = append(data.zh, 0.5)
    dz = zh[1::] - zh[0:-1]

    z = append(-data.z[0], data.z)
    dzh = z[1::] - z[0:-1]

    self.u = 0.
    self.w = 0.
    self.p = 0.

    for k in range(data.z.size):
      self.u = self.u + sum( (dx*dz [k]*abs(data.u[k,:] - ref.u[k,:])) )
      self.w = self.w + sum( (dx*dzh[k]*abs(data.w[k,:] - ref.w[k,:])) )
      self.p = self.p + sum( (dx*dz [k]*abs(data.p[k,:] - ref.p[k,:])) )
