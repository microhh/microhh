import numpy as np
import struct as st

class Grid_cell:
    def __init__(self, i, j, k, x=None, y=None, z=None):
        self.i = i
        self.j = j
        self.k = k

        self.x = x
        self.y = y
        self.z = z

class Neighbours:
    def __init__(self):
        self.i = []
        self.j = []
        self.k = []
        self.x = []
        self.y = []
        self.z = []
        self.d = []

    def add(self, i, j, k, x, y, z, d):
        self.i.append(i)
        self.j.append(j)
        self.k.append(k)
        self.x.append(x)
        self.y.append(y)
        self.z.append(z)
        self.d.append(d)

    def sort(self, clip=None):
        index = np.array(self.d).argsort()
        n = index.size if clip is None else clip

        assert n < index.size

        self.i = np.array(self.i)[index][:n]
        self.j = np.array(self.j)[index][:n]
        self.k = np.array(self.k)[index][:n]
        self.x = np.array(self.x)[index][:n]
        self.y = np.array(self.y)[index][:n]
        self.z = np.array(self.z)[index][:n]
        self.d = np.array(self.d)[index][:n]

class Ghost_cell(Grid_cell):
    def __init__(self, i, j, k, x=None, y=None, z=None):
        Grid_cell.__init__(self, i, j, k, x, y, z)
        self.neighbours = Neighbours()

    def set_boundary_location(self, x, z):
        self.xb = x
        self.zb = z

class Immersed_boundary:
    """
    Base class for defining the immersed boundary
    """
    def __init__(self, grid):
        self.grid = grid
        self.ghost_cells_u = []
        self.ghost_cells_v = []
        self.ghost_cells_w = []
        self.ghost_cells_s = []

    def calculate(self):
        # u-component:
        self.find_ghost_cells(grid.xh, grid.z, self.ghost_cells_u)
        self.find_fluid_cells(grid.xh, grid.z, self.ghost_cells_u, n=2)
        self.write_output(self.ghost_cells_u, 'ib.boundary_u')

        # w-component:
        self.find_ghost_cells(grid.x, grid.zh, self.ghost_cells_w)
        self.find_fluid_cells(grid.x, grid.zh, self.ghost_cells_w, n=2)
        self.write_output(self.ghost_cells_w, 'ib.boundary_w')

        # scalar
        self.find_ghost_cells(grid.x, grid.z, self.ghost_cells_s)
        self.find_fluid_cells(grid.x, grid.z, self.ghost_cells_s, n=2)
        self.write_output(self.ghost_cells_s, 'ib.boundary_s')

    def find_ghost_cells(self, x, z, ghost_cells):
        dx = x[1] - x[0]

        for i in range(grid.istart,grid.iend):
            zb  = self.boundary_function(x[i  ])
            zbm = self.boundary_function(x[i-1])
            zbp = self.boundary_function(x[i+1])

            for k in range(grid.kstart, grid.kend):
                if (z[k] < zb and z[k+1] > zb  or 
                    z[k] < zb and z[k  ] > zbm or
                    z[k] < zb and z[k  ] > zbp):

                    g = Ghost_cell(i, 1, k, x[i], None, z[k])
                    d, xb, zb = self.nearest_boundary(x[i], z[k], x[i]-dx, x[i]+dx)
                    g.set_boundary_location(xb, zb)
                    ghost_cells.append(g)

    # For each IB ghost cell, find the nearest grid points outside of the boundary
    def find_fluid_cells(self, x, z, ghost_cells, n):
        dx = x[1] - x[0]

        for ghost_cell in ghost_cells:
            i = ghost_cell.i
            k = ghost_cell.k
            
            # BvS fix this for cases with more than 1 ghost cell
            min_di = 0 if i == 0        else -1
            max_di = 0 if i == x.size-1 else +1

            for di in range(min_di, max_di+1):
                zb = self.boundary_function(x[i+di])
                for dk in range(-1,3):
                    if (z[k+dk] > zb):
                        ii0 = max(0, i+di-1)
                        ii1 = min(i+di+1, x.size-1)
                        min_dist, x_tmp, z_tmp = self.nearest_boundary(x[i+di], z[k+dk], x[ii0], x[ii1])

                        # Estimate dz:
                        dz = z[k+dk+1] - z[k+dk]
                        if(min_dist > 0.25*min(dx,dz)):
                            distance = np.sqrt( (x[i]-x[i+di])**2 + (z[k]-z[k+dk])**2 )
                            ghost_cell.neighbours.add(i+di, 1, k+dk, x[i+di], 1, z[k+dk], distance)

            ghost_cell.neighbours.sort(clip=n)

    # Find the nearest location on the boundary
    def nearest_boundary(self, x, z, x0, x1):
        min_dist    = 1e12
        x_boundary  = 0
        z_boundary  = 0
      
        xx = np.linspace(x0, x1, 256)
        zb = self.boundary_function(xx)
        d  = np.sqrt((x-xx)**2+(z-zb)**2)

        return d.min(), xx[d.argmin()], zb[d.argmin()]

    def write_output(self, ghost_cells, file_name):
        f = open(file_name, 'w')
       
        f.write('{0:^5s} {1:^5s} {2:^18s} {3:^18s} '.format('i', 'k', 'xb', 'zb'))
        for i in range(ghost_cells[0].neighbours.i.size):
            f.write('{0:^5s} {1:^5s} '.format('ii','kk'))
        f.write('\n')

        for g in ghost_cells:
            f.write('{0:05d} {1:05d} {2:1.12E} {3:1.12E} '.format(g.i-1, g.k-1, g.xb, g.zb))
            for i in range(g.neighbours.i.size):
                f.write('{0:05d} {1:05d} '.format(g.neighbours.i[i]-1, g.neighbours.k[i]-1))
            f.write('\n')
        f.close()

    def plot(self, ghost_cells):
        pl.figure()
        ax = pl.gca()

        # Plot the immersed boundary
        x = np.linspace(0, self.grid.xsize, 256)
        y = self.boundary_function(x)
        ax.plot(x, y, linewidth=1.5, color='k')

        for g in ghost_cells:
            # Mark the ghost_cells
            pl.scatter(g.x, g.z, marker='o', s=20)
            
            # Mark the nearest location on the boundary
            pl.plot([g.x, g.xb],[g.z, g.zb], 'k-')

            # Mark the fluid cells:
            for x,z in zip(g.neighbours.x, g.neighbours.z):
                pl.scatter(x, z, marker='x')
                pl.plot([g.x, x],[g.z, z], 'k:')

        # Add half level grid lines
        ax.set_xticks(grid.xh)
        ax.set_yticks(grid.zh)
        ax.grid()
        ax.set_xlim(0, self.grid.xsize)
        ax.set_ylim(0, self.grid.zsize) 

class Immersed_gaussian(Immersed_boundary):
    """
    Gaussian shaped immersed boundary
    x0 = center of hill (m)
    A = amplitude of hill (m)
    sigma_x = width (std.dev) of hill (m)  
    """
    def __init__(self, grid, x0, A, sigma_x, z_offset=0):
        Immersed_boundary.__init__(self, grid)

        # Settings Gaussian hill
        self.x0 = x0             # Center (x) hill (m)
        self.A = A               # Amplitude hill (m)
        self.sigma_x = sigma_x   # Standard deviation (width) hill (m)
        self.z_offset = z_offset # Vertical offset hill (m)

    def boundary_function(self, x):
        """
        Function describing the shape of the hill
        """
        return self.z_offset + self.A * np.exp(-((x-self.x0)/(2*self.sigma_x))**2)

class Immersed_sine(Immersed_boundary):
    """
    Sine shaped immersed boundary
    x0 = center of hill (m)
    A = amplitude of hill (m)
    sigma_x = width (std.dev) of hill (m)  
    """
    def __init__(self, grid, A, t, z_offset=0):
        Immersed_boundary.__init__(self, grid)
        self.A = A               # Amplitude hills (m)
        self.t = t               # Wavelength hills (m) 
        self.z_offset = z_offset # Vertical offset hill (m)

        print(self.z_offset)

    def boundary_function(self, x):
        """
        Function describing the shape of the hills
        """
        return self.z_offset + self.A + self.A * np.sin(2*np.pi*x/self.t)

if __name__ == "__main__":
    import matplotlib.pylab as pl
    from microhh_tools import *

    pl.close('all')

    ini  = Read_namelist()
    grid = Read_grid(ini.grid.itot, ini.grid.jtot, ini.grid.ktot, ini.grid.zsize, n_ghost=1)  
   
    ib = Immersed_sine(grid, A=0.0025, t=0.05, z_offset=0.00016)
    ib.calculate()
    #ib.plot(ib.ghost_cells_u)


