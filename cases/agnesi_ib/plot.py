import numpy as np
import matplotlib.pylab as pl
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from microhh_tools import *

pl.ion()
pl.close('all')

nl = Read_namelist()
gr = Read_grid(nl.grid.itot, nl.grid.jtot, nl.grid.ktot, nl.grid.zsize)

class Ref:
    def __init__(self, file_name):
        f = np.loadtxt(file_name, delimiter=',')
        self.x = f[:,0]
        self.y = f[:,1]

class Debug:
    def __init__(self, file_name, x, y, z):
        f = np.loadtxt(file_name, delimiter=',')
        i = 0
        self.i  = np.array([x[j] for j in f[:,i]]); i+=1
        self.j  = np.array([y[j] for j in f[:,i]]); i+=1
        self.k  = np.array([z[j] for j in f[:,i]]); i+=1

        self.xb = f[:,i]; i+=1 
        self.yb = f[:,i]; i+=1
        self.zb = f[:,i]; i+=1

        self.i1 = np.array([x[j] for j in f[:,i]]); i+=1
        self.j1 = np.array([y[j] for j in f[:,i]]); i+=1
        self.k1 = np.array([z[j] for j in f[:,i]]); i+=1
        self.i2 = np.array([x[j] for j in f[:,i]]); i+=1
        self.j2 = np.array([y[j] for j in f[:,i]]); i+=1
        self.k2 = np.array([z[j] for j in f[:,i]]); i+=1
        self.i3 = np.array([x[j] for j in f[:,i]]); i+=1
        self.j3 = np.array([y[j] for j in f[:,i]]); i+=1
        self.k3 = np.array([z[j] for j in f[:,i]]); i+=1
        self.i4 = np.array([x[j] for j in f[:,i]]); i+=1
        self.j4 = np.array([y[j] for j in f[:,i]]); i+=1
        self.k4 = np.array([z[j] for j in f[:,i]])
        self.k4 = f[:,i]

if(True):
    # Immersed boundary shape
    ib = nl.immersed_boundary
    xb = np.linspace(0,6000,256)
    yb = ib.amplitude / (1. + pow((xb-ib.x0_hill)/ib.sigma_x_hill, 2))

    ds = Debug('debug_s.txt', gr.x, gr.y, gr.z)

if(False):
    # Read MicroHH data from restart files
    u = read_3d('u.0003600', nl.grid.itot, nl.grid.jtot, nl.grid.ktot)
    v = read_3d('v.0003600', nl.grid.itot, nl.grid.jtot, nl.grid.ktot)
    w = read_3d('w.0003600', nl.grid.itot, nl.grid.jtot, nl.grid.ktot)

    # Immersed boundary shape
    ib = nl.immersed_boundary
    xb = np.linspace(0,6000,256)
    yb = ib.amplitude / (1. + pow((xb-ib.x0_hill)/ib.sigma_x_hill, 2))

    for x in np.arange(0,6000.01,1000):
        ii = np.abs(gr.x-x).argmin()
        jj = 0

        xx = gr.x[ii]
        yy = gr.y[jj]

        z_ib = ib.amplitude / (1. + pow((xx-ib.x0_hill)/ib.sigma_x_hill, 2))
        
        pl.figure()
        pl.plot(u[:,jj,ii], gr.z-200, label=str(x))
        pl.plot([-1,3],[z_ib,z_ib])
        pl.grid()
        pl.ylim(z_ib-50,z_ib+50)
        pl.xlim(-0.2,0.2)



if(False):
    # Read MicroHH data from restart files
    u = read_3d('u.0259200', nl.grid.itot, nl.grid.jtot, nl.grid.ktot)
    v = read_3d('v.0259200', nl.grid.itot, nl.grid.jtot, nl.grid.ktot)
    w = read_3d('w.0259200', nl.grid.itot, nl.grid.jtot, nl.grid.ktot)
    
    # Read reference data from Lundquist (2012)
    u1000 = Ref('reference_data/ref_luq_u_1000.txt')
    u2000 = Ref('reference_data/ref_luq_u_2000.txt')
    u3000 = Ref('reference_data/ref_luq_u_3000.txt')
    u5000 = Ref('reference_data/ref_luq_u_5000.txt')
    
    v2000 = Ref('reference_data/ref_luq_v_2000.txt')
    v3000 = Ref('reference_data/ref_luq_v_3000.txt')
    
    w3000 = Ref('reference_data/ref_luq_w_3000.txt')
    
    
    
    #def fake_subplot(x, y, x0, scale=1, color='k', linetype='-'):
    #    pl.plot(x0+x*scale, y, color=color, ls=linetype) 
    
    
    def plot_case(x1, y1, x2, y2, x_grid, y_grid, x0, y0, scale, zoffs):
        # Zero-line
        pl.plot([x0,x0], [0,y1[-1]], 'k:')
    
        # Get location in MicroHH 3D field, and plot field
        i = np.abs(x_grid - x0).argmin()
        j = np.abs(y_grid - y0).argmin()
        pl.plot(x0+x1[:,j,i]*scale, y1+zoffs, 'k-')
        print(x_grid[i], y_grid[j])
    
        # plot reference case:
        pl.plot(x0+x2*scale, y2, 'x')
    
    # Immersed boundary shape
    ib = nl.immersed_boundary
    xb = np.linspace(0,6000,256)
    yb = 2953. 
    yb = ib.amplitude / (1. + pow((xb-ib.x0_hill)/ib.sigma_x_hill, 2) \
                            + pow((yb-ib.y0_hill)/ib.sigma_y_hill, 2))
    
    pl.figure()
    
    #----------------------------------------
    ax=pl.subplot(311)
    remove_top_right_ax()
    pl.plot(xb, yb, 'k-', linewidth=1.5)
    pl.ylim(0,3000)
    
    plot_case(u, gr.z, u1000.x, u1000.y, gr.xh, gr.y, 1000, 3000, 60, -200)
    plot_case(u, gr.z, u2000.x, u2000.y, gr.xh, gr.y, 2000, 3000, 60, -200)
    plot_case(u, gr.z, u3000.x, u3000.y, gr.xh, gr.y, 3000, 3000, 60, -200)
    plot_case(u, gr.z, u3000.x, u3000.y, gr.xh, gr.y, 4000, 3000, 60, -200)
    plot_case(u, gr.z, u5000.x, u5000.y, gr.xh, gr.y, 5000, 3000, 60, -200)
    
    ax=pl.subplot(312)
    remove_top_right_ax()
    pl.plot(xb, yb, 'k-', linewidth=1.5)
    pl.ylim(0,3000)
    
    plot_case(v, gr.z, v2000.x, v2000.y, gr.x, gr.yh, 1000, 3000, 100, -200)
    plot_case(v, gr.z, v2000.x, v2000.y, gr.x, gr.yh, 2000, 3000, 100, -200)
    plot_case(v, gr.z, v3000.x, v3000.y, gr.x, gr.yh, 3000, 3000, 100, -200)
    plot_case(v, gr.z, v2000.x, v2000.y, gr.x, gr.yh, 4000, 3000, 100, -200)
    plot_case(v, gr.z, v2000.x, v2000.y, gr.x, gr.yh, 5000, 3000, 100, -200)
    
    
    
    """
    x0 = 3000
    y0 = 3000
    scale = 60
    zoffs = -200
    
    # 0-line
    pl.plot([x0,x0],[0,3000], 'k:')
    
    # data
    pl.plot(1000 + u1000.x*scale, u1000.y, 'x')
    pl.plot(x0 + u3000.x*scale, u3000.y, 'x')
    pl.plot(5000 + u5000.x*scale, u5000.y, 'x')
    
    
    
    i = np.abs(gr.xh - x0).argmin()
    j = np.abs(gr.y  - y0).argmin()
    pl.plot(x0+u[:,j,i]*scale, gr.z+zoffs, 'k-')
    
    # fake ax
    x = np.linspace(0,10,3)
    pl.plot([x0+x[0]*scale, x0+x[-1]*scale], [3000,3000], 'k')
    for xx in x:
        pl.text(x0+xx*scale, 3000, str(xx), ha='center', va='bottom', size=8)
    
    #----------------------------------------
    pl.subplot(312)
    remove_top_right_ax()
    pl.plot(xb, yb, 'k-', linewidth=1.5)
    pl.ylim(0,3000)
    
    x0 = 3000
    y0 = 3000
    scale = 100
    zoffs = -200
    
    # 0-line
    pl.plot([x0,x0],[0,3000], 'k:')
    
    # data
    pl.plot(x0 + v3000.x*scale, v3000.y, 'x')
    i = np.abs(gr.x - x0).argmin()
    j = np.abs(gr.yh- y0).argmin()
    pl.plot(x0+v[:,j,i]*scale, gr.z+zoffs, 'k-')
    
    # fake ax
    x = np.linspace(-1,3,3)
    pl.plot([x0+x[0]*scale, x0+x[-1]*scale], [3000,3000], 'k')
    for xx in x:
        pl.text(x0+xx*scale, 3000, str(xx), ha='center', va='bottom', size=8)
    
    
    
    #----------------------------------------
    pl.subplot(313)
    remove_top_right_ax()
    pl.plot(xb, yb, 'k-', linewidth=1.5)
    pl.ylim(0,3000)
    
    x0 = 3000
    y0 = 3000
    scale = 1000
    zoffs = -200
    
    # 0-line
    pl.plot([x0,x0],[0,3000], 'k:')
    
    # data
    pl.plot(x0 + w3000.x*scale, w3000.y, 'x')
    i = np.abs(gr.x - x0).argmin()
    j = np.abs(gr.y - y0).argmin()
    pl.plot(x0+w[:,j,i]*scale, gr.z+zoffs, 'k-')
    
    # fake ax
    x = np.array([-0.1,0.3])
    pl.plot([x0+x[0]*scale, x0+x[-1]*scale], [3000,3000], 'k')
    for xx in x:
        pl.text(x0+xx*scale, 3000, str(xx), ha='center', va='bottom', size=8)
    
    
    pl.savefig('profiles.pdf')
    
    #pl.subplot(131)
    #pl.plot(u3000.x, u3000.y, 'x')
    #
    #pl.subplot(132)
    #pl.plot(v3000.x, v3000.y, 'x')
    #i = np.abs(gr.x  - 3000).argmin()
    #j = np.abs(gr.x  - 3000).argmin()
    #pl.plot(v[:,j,i], gr.z-200)
    """
