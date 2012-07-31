import numpy
import struct

class microhh:
  def __init__(self, iter, itot, jtot, ktot):
    nx = itot
    ny = jtot
    nz = ktot
    
    n  = nx*ny*nz
    
    fin = open("grid.{:07d}".format(0),"rb")
    raw = fin.read(nx*8)
    self.x = numpy.array(struct.unpack('<{}d'.format(nx), raw))
    raw = fin.read(nx*8)
    self.xh = numpy.array(struct.unpack('<{}d'.format(nx), raw))
    raw = fin.read(ny*8)
    self.y = numpy.array(struct.unpack('<{}d'.format(ny), raw))
    raw = fin.read(ny*8)
    self.yh = numpy.array(struct.unpack('<{}d'.format(ny), raw))
    raw = fin.read(nz*8)
    self.z = numpy.array(struct.unpack('<{}d'.format(nz), raw))
    raw = fin.read(nz*8)
    self.zh = numpy.array(struct.unpack('<{}d'.format(nz), raw))
    fin.close()
    
    fin = open("u.{:07d}".format(iter),"rb")
    raw = fin.read(n*8)
    tmp = numpy.array(struct.unpack('<{}d'.format(n), raw))
    del(raw)
    self.u = tmp.reshape((nz, ny, nx))
    del(tmp)
    fin.close()
    
    fin = open("v.{:07d}".format(iter),"rb")
    raw = fin.read(n*8)
    tmp = numpy.array(struct.unpack('<{}d'.format(n), raw))
    del(raw)
    self.v = tmp.reshape((nz, ny, nx))
    del(tmp)
    fin.close()
    
    fin = open("w.{:07d}".format(iter),"rb")
    raw = fin.read(n*8)
    tmp = numpy.array(struct.unpack('<{}d'.format(n), raw))
    del(raw)
    self.w = tmp.reshape((nz, ny, nx))
    del(tmp)
    fin.close()
    
    fin = open("p.{:07d}".format(iter),"rb")
    raw = fin.read(n*8)
    tmp = numpy.array(struct.unpack('<{}d'.format(n), raw))
    del(raw)
    self.p = tmp.reshape((nz, ny, nx))
    del(tmp)
    fin.close()

    fin = open("s.{:07d}".format(iter),"rb")
    raw = fin.read(n*8)
    tmp = numpy.array(struct.unpack('<{}d'.format(n), raw))
    del(raw)
    self.s = tmp.reshape((nz, ny, nx))
    del(tmp)
    fin.close()

