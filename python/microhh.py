import numpy
import struct

class microhh:
  def __init__(self, iter, itot, jtot, ktot):
    nx = itot
    ny = jtot
    nz = ktot
    
    n  = nx*ny*nz
    
    fin = open("grid.{:06d}".format(0),"rb")
    raw = fin.read(nx*8)
    self.x = numpy.array(struct.unpack('{}d'.format(nx), raw))
    raw = fin.read(nx*8)
    self.xh = numpy.array(struct.unpack('{}d'.format(nx), raw))
    raw = fin.read(ny*8)
    self.y = numpy.array(struct.unpack('{}d'.format(ny), raw))
    raw = fin.read(ny*8)
    self.yh = numpy.array(struct.unpack('{}d'.format(ny), raw))
    raw = fin.read(nz*8)
    self.z = numpy.array(struct.unpack('{}d'.format(nz), raw))
    raw = fin.read(nz*8)
    self.zh = numpy.array(struct.unpack('{}d'.format(nz), raw))
    fin.close()
    
    fin = open("u.{:06d}".format(iter),"rb")
    raw = fin.read(n*8)
    tmp = numpy.array(struct.unpack('{}d'.format(n), raw))
    self.u = tmp.reshape((nz, ny, nx))
    fin.close()
    
    fin = open("v.{:06d}".format(iter),"rb")
    raw = fin.read(n*8)
    tmp = numpy.array(struct.unpack('{}d'.format(n), raw))
    self.v = tmp.reshape((nz, ny, nx))
    fin.close()
    
    fin = open("w.{:06d}".format(iter),"rb")
    raw = fin.read(n*8)
    tmp = numpy.array(struct.unpack('{}d'.format(n), raw))
    self.w = tmp.reshape((nz, ny, nx))
    fin.close()
    
    fin = open("p.{:06d}".format(iter),"rb")
    raw = fin.read(n*8)
    tmp = numpy.array(struct.unpack('{}d'.format(n), raw))
    self.p = tmp.reshape((nz, ny, nx))
    fin.close()

    fin = open("s.{:06d}".format(iter),"rb")
    raw = fin.read(n*8)
    tmp = numpy.array(struct.unpack('{}d'.format(n), raw))
    self.s = tmp.reshape((nz, ny, nx))
    fin.close()

