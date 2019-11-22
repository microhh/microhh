import numpy as np

class Grid:
    def __init__(self, xsize, ysize, zsize, itot, jtot, ktot, ijgc=1, kgc=1, TF=np.float64):
        """
        Simple (equidistant) grid
        """

        np.random.seed(666)

        self.xsize = TF(xsize)
        self.ysize = TF(ysize)
        self.zsize = TF(zsize)

        self.itot = np.int32(itot)
        self.jtot = np.int32(jtot)
        self.ktot = np.int32(ktot)

        self.igc = np.int32(ijgc)
        self.jgc = np.int32(ijgc)
        self.kgc = np.int32(kgc)

        self.icells = np.int32(itot+2*ijgc)
        self.jcells = np.int32(jtot+2*ijgc)
        self.kcells = np.int32(ktot+2*kgc)

        self.ijcells = np.int32(self.icells*self.jcells)
        self.ncells  = np.int32(self.icells*self.jcells*self.kcells)

        self.istart = np.int32(self.igc)
        self.jstart = np.int32(self.jgc)
        self.kstart = np.int32(self.kgc)

        self.iend = np.int32(itot+self.igc)
        self.jend = np.int32(jtot+self.jgc)
        self.kend = np.int32(ktot+self.kgc)

        self.dx = TF(self.xsize / self.itot)
        self.dy = TF(self.ysize / self.jtot)
        self.dz = np.random.random(self.kcells).astype(TF)

        self.dxi = TF(1/self.dx)
        self.dyi = TF(1/self.dy)

        self.dzi   = np.random.random(self.kcells).astype(TF)
        self.dzi4  = np.random.random(self.kcells).astype(TF)
        self.dzhi4 = np.random.random(self.kcells).astype(TF)



class Field3d:
    def __init__(self, ncells, ijcells, TF=np.float64):
        """
        Simple 3D field incl. some surface fields
        """

        self.fld  = np.random.random(ncells).astype(TF)
        self.tend = np.zeros(ncells, dtype=TF)

        self.fld_bot = np.random.random(ijcells).astype(TF)
        self.fld_top = np.random.random(ijcells).astype(TF)

        self.flux_bot = np.random.random(ijcells).astype(TF)
        self.flux_top = np.random.random(ijcells).astype(TF)


class Fields:
    def __init__(self, fields, ncells, ijcells, kcells, TF=np.float64):

        np.random.seed(666)

        for field in fields:
            setattr(self, field, Field3d(ncells, ijcells, TF))

        self.rhoref  = np.random.random(kcells).astype(TF)
        self.rhorefh = np.random.random(kcells).astype(TF)


if __name__ == '__main__':
    TF = np.float32
    grid   = Grid(3200, 3200, 3200, 32, 32, 32, 2, 1, TF)
    fields = Fields(['u','v','w','s'], grid.ncells, grid.ijcells, grid.kcells, TF)
