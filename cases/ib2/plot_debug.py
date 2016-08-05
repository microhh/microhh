import numpy as np
import matplotlib.pylab as pl
from mpl_toolkits.mplot3d import Axes3D
import itertools

pl.close('all')

from microhh_tools import *

class Read_debug:
    def __init__(self, filename='debug.txt'):
        f = np.loadtxt(filename, delimiter=',') 
        self.i  = f[:,0].astype(np.int)
        self.j  = f[:,1].astype(np.int)
        self.k  = f[:,2].astype(np.int)
        self.xb = f[:,3]
        self.yb = f[:,4]
        self.zb = f[:,5]

        self.ni = np.zeros((self.i.size, 3))
        self.nj = np.zeros((self.j.size, 3))
        self.nk = np.zeros((self.k.size, 3))

        self.ni[:,0] = f[:,6]
        self.nj[:,0] = f[:,7]
        self.nk[:,0] = f[:,8]

        self.ni[:,1] = f[:,9]
        self.nj[:,1] = f[:,10]
        self.nk[:,1] = f[:,11]

        self.ni[:,2] = f[:,12]
        self.nj[:,2] = f[:,13]
        self.nk[:,2] = f[:,14]

# ----- Read debug.txt -----
db  = Read_debug("debug_s.txt")

# ----- Read MicroHH namelist -----
nl = Read_namelist()
ib = nl.immersed_boundary

# ----- Read MicroHH grid -----
grid = Read_grid(nl.grid.itot, nl.grid.jtot, nl.grid.ktot, nl.grid.zsize, n_ghost=2) 


# 2d 2d 2d 2d 2d 2d -----------------
if(True):
    pl.figure()
    ax = pl.subplot(111) 
    
    # ----- Plot boundary ----- 
    x = np.linspace(0, nl.grid.xsize, 256)
    if(ib.sw_ib == 'sine'):
        z = ib.z_offset + ib.amplitude + ib.amplitude * np.sin(2*np.pi*x/ib.wavelength_x)
    elif(ib.sw_ib == 'gaussian'):
        z = ib.z_offset + ib.amplitude * np.exp(-pow((x-ib.x0_hill)/(2*ib.sigma_x_hill), 2));
    
    ax.plot(x, z, color='k', alpha=0.3)
    
    # ----- Scatter ghost cells -----
    x = grid.x
    y = grid.y
    z = grid.z
   
    colors = itertools.cycle(cc) 
    for i,j,k,xb,yb,zb in zip(db.i, db.j, db.k, db.xb, db.yb, db.zb):
        color = next(colors)
        ax.scatter(x[i], z[k], facecolor=color, edgecolor='none')
        ax.plot([x[i],xb], [z[k],zb], ':', color=color)
    
    #colors = itertools.cycle(cc) 
    for i,j,k,ni,nj,nk in zip(db.i, db.j, db.k, db.ni, db.nj, db.nk):
        color = next(colors)
        for n in range(3):
            ax.scatter(x[ni[n]], z[nk[n]], marker='x', facecolor=color, edgecolor=color)
            pl.plot([x[i],x[ni[n]]], [z[k],z[nk[n]]], '-', color=color)

    ax.set_xticks(grid.xh)
    ax.set_yticks(grid.zh)
    ax.grid()

# 3d 3d 3d 3d 3d 3d -----------------
if(True):
    pl.figure()
    ax = pl.subplot(111, projection='3d') 
    
    # ----- Plot boundary ----- 
    x = np.linspace(0, nl.grid.xsize, 256)
    if(ib.sw_ib == 'sine'):
        z = ib.z_offset + ib.amplitude + ib.amplitude * np.sin(2*np.pi*x/ib.wavelength_x)
    elif(ib.sw_ib == 'gaussian'):
        z = ib.z_offset + ib.amplitude * np.exp(-pow((x-ib.x0_hill)/(2*ib.sigma_x_hill), 2));
    
    y = grid.y
    xx,yy = np.meshgrid(x,y)
    zz = np.repeat(z[np.newaxis,:], y.size, axis=0) 
    ax.plot_surface(xx,yy,zz,color='y',alpha=0.3)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    # ----- Scatter ghost cells -----
    x = grid.x
    y = grid.y
    z = grid.z
    
    for i,j,k,xb,yb,zb in zip(db.i, db.j, db.k, db.xb, db.yb, db.zb):
        ax.scatter(x[i], y[j], z[k], zorder=10)
        ax.plot([x[i],xb],[y[j],yb],[z[k],zb], 'k:', zorder=11)
    
    for i,j,k,ni,nj,nk in zip(db.i, db.j, db.k, db.ni, db.nj, db.nk):
        for n in range(3):
            ax.scatter(x[ni[n]], y[nj[n]], z[nk[n]], marker='x')
            pl.plot([x[i],x[ni[n]]], [y[j],y[nj[n]]], [z[k],z[nk[n]]], 'b-')
