import numpy as np
import matplotlib.pylab as pl
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib import rc
from microhh_tools import *

from matplotlib import rc
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{fourier},\usepackage{amsmath}')
rc('font', size=14)
rc('legend', fontsize=11)
rc('xtick', direction='out')
rc('ytick', direction='out')
rc('xtick.major', size=5)
rc('xtick.minor', size=2)
rc('ytick.major', size=5)
rc('ytick.minor', size=2)

pl.ion()
pl.close('all')

igc = 2
jgc = 2
kgc = 1

nl = Read_namelist()
gr = Read_grid(nl.grid.itot, nl.grid.jtot, nl.grid.ktot, nl.grid.zsize)

class Ref:
    def __init__(self, file_name):
        f = np.loadtxt(file_name, delimiter=',')
        self.x = f[:,0]
        self.y = f[:,1]

# Read the IB debug output from MicroHH
class Debug:
    def __init__(self, file_name, x, y, z):
        f = np.loadtxt(file_name, delimiter=',')
        i = 0
        self.x  = np.array([x[j+igc] for j in f[:,i].astype(np.int)], dtype=np.int); i+=1
        self.y  = np.array([y[j+jgc] for j in f[:,i].astype(np.int)], dtype=np.int); i+=1
        self.z  = np.array([z[j+kgc] for j in f[:,i].astype(np.int)], dtype=np.int); i+=1

        self.xb = f[:,i]; i+=1
        self.yb = f[:,i]; i+=1
        self.zb = f[:,i]; i+=1

        self.x1 = np.array([x[j+igc] for j in f[:,i].astype(np.int)], dtype=np.int); i+=1
        self.y1 = np.array([y[j+jgc] for j in f[:,i].astype(np.int)], dtype=np.int); i+=1
        self.z1 = np.array([z[j+kgc] for j in f[:,i].astype(np.int)], dtype=np.int); i+=1
        self.x2 = np.array([x[j+igc] for j in f[:,i].astype(np.int)], dtype=np.int); i+=1
        self.y2 = np.array([y[j+jgc] for j in f[:,i].astype(np.int)], dtype=np.int); i+=1
        self.z2 = np.array([z[j+kgc] for j in f[:,i].astype(np.int)], dtype=np.int); i+=1
        self.x3 = np.array([x[j+igc] for j in f[:,i].astype(np.int)], dtype=np.int); i+=1
        self.y3 = np.array([y[j+jgc] for j in f[:,i].astype(np.int)], dtype=np.int); i+=1
        self.z3 = np.array([z[j+kgc] for j in f[:,i].astype(np.int)], dtype=np.int); i+=1
        self.x4 = np.array([x[j+igc] for j in f[:,i].astype(np.int)], dtype=np.int); i+=1
        self.y4 = np.array([y[j+jgc] for j in f[:,i].astype(np.int)], dtype=np.int); i+=1
        self.z4 = np.array([z[j+kgc] for j in f[:,i].astype(np.int)], dtype=np.int)

def add_sub_ax(fig, parent_ax, x0, y0, width, xlim=None, x00=None):
    # Transform from data coordinates to pixel coordinates:
    t = parent_ax.transData.transform((x0, y0))
    # Transform from pixel coordinates to figure (0..1) coordinates:
    t = fig.transFigure.inverted().transform(t)

    if xlim is not None and x00 is not None:
        xoffs = -(x00-xlim[0])/(xlim[1]-xlim[0])*width
    else:
        xoffs = 0

    # Add new axis
    new_ax = fig.add_axes([t[0]+xoffs, t[1], width, parent_ax.get_position().height])

    # Set transparent background color
    new_ax.set_axis_bgcolor('none')

    # Set y-extent equal to parent
    new_ax.set_ylim(parent_ax.get_ylim())

    # Remove y-ticks and labels
    new_ax.yaxis.set_ticks([])

    # Remove axis, except the top one
    new_ax.spines['right'].set_visible(False)
    new_ax.spines['left'].set_visible(False)
    new_ax.spines['bottom'].set_visible(False)
    new_ax.get_xaxis().tick_top()

    # Set label position to top
    new_ax.xaxis.set_label_position('top')

    # Remove top and right axis of parent ax
    parent_ax.spines['right'].set_visible(False)
    parent_ax.get_yaxis().tick_left()
    parent_ax.spines['top'].set_visible(False)
    parent_ax.get_xaxis().tick_bottom()

    if xlim is not None:
        new_ax.set_xlim(xlim)

    # Return new axis object for further manipulation
    return new_ax

def get_loc_3d(x, y, x0, y0):
    i = np.abs(x-x0).argmin()
    j = np.abs(y-y0).argmin()
    return i,j

if(True):
    # Read MicroHH data from restart files
    u = read_3d('u.0259200', nl.grid.itot, nl.grid.jtot, nl.grid.ktot)
    v = read_3d('v.0259200', nl.grid.itot, nl.grid.jtot, nl.grid.ktot)
    w = read_3d('w.0259200', nl.grid.itot, nl.grid.jtot, nl.grid.ktot)

    # Read reference data from Lundquist (2012)
    u1000_bf = Ref('reference_data/u_1000_bf.txt')
    u1000_id = Ref('reference_data/u_1000_id.txt')

    u3000_bf = Ref('reference_data/u_3000_bf.txt')
    u3000_id = Ref('reference_data/u_3000_id.txt')

    u5000_bf = Ref('reference_data/u_5000_bf.txt')
    u5000_id = Ref('reference_data/u_5000_id.txt')

    # Boundary shape
    ib = nl.immersed_boundary
    xb = np.linspace(0,6000,256)
    yb = 2953.
    yb = ib.amplitude / (1. + pow((xb-ib.x0_hill)/ib.sigma_x_hill, 2) \
                            + pow((yb-ib.y0_hill)/ib.sigma_y_hill, 2))

    fig = pl.figure(figsize=(12,5))
    fig.subplots_adjust(0.08, 0.11,0.96,0.89,0,0)

    ax=pl.subplot(111)
    ax.plot(xb, yb, 'k-', linewidth=1.5)
    ax.fill_between(xb, 0, yb, color='0.8')
    ax.set_ylim(0,3000)
    ax.set_xlim(0,6000)
    ax.set_ylabel('z (m)')
    ax.set_xlabel('x (m)')

    # Dummy legend
    ax.plot(-1,-1, 'k-', label='MicroHH')
    ax.plot(-1,-1, 'x', label='WRF terrain following')
    ax.plot(-1,-1, 'o', mfc='none', label='WRF IBM-IDW')
    pl.legend(frameon=False, loc=2)

    ax1000 = add_sub_ax(fig, ax, 1000., 0., 0.1, xlim=[-1,11], x00=0)
    ax1000.plot(u1000_bf.x, u1000_bf.y, 'x')
    ax1000.plot(u1000_id.x, u1000_id.y, 'o', mfc='none')
    i,j = get_loc_3d(gr.xh, gr.y, 1000., 3000.)
    ax1000.plot(u[:,j,i], gr.z-200, 'k-')
    ax1000.set_xlabel('u (m/s)')
    ax1000.grid()
    ax1000.set_xticks([0,5,10])

    ax2000 = add_sub_ax(fig, ax, 2000., 0., 0.1, xlim=[-1,11], x00=0)
    i,j = get_loc_3d(gr.xh, gr.y, 2000., 3000.)
    ax2000.plot(u[:,j,i], gr.z-200, 'k-')
    ax2000.set_xlabel('u (m/s)')
    ax2000.grid()
    ax2000.set_xticks([0,5,10])

    ax3000 = add_sub_ax(fig, ax, 3000., 0., 0.1, xlim=[-1,11], x00=0)
    ax3000.plot(u3000_bf.x, u3000_bf.y, 'x')
    ax3000.plot(u3000_id.x, u3000_id.y, 'o', mfc='none')
    i,j = get_loc_3d(gr.xh, gr.y, 3000., 3000.)
    ax3000.plot(u[:,j,i], gr.z-200, 'k-')
    ax3000.set_xlabel('u (m/s)')
    ax3000.grid()
    ax3000.set_xticks([0,5,10])

    ax4000 = add_sub_ax(fig, ax, 4000., 0., 0.1, xlim=[-1,11], x00=0)
    i,j = get_loc_3d(gr.xh, gr.y, 4000., 3000.)
    ax4000.plot(u[:,j,i], gr.z-200, 'k-')
    ax4000.set_xlabel('u (m/s)')
    ax4000.grid()
    ax4000.set_xticks([0,5,10])

    ax5000 = add_sub_ax(fig, ax, 5000., 0., 0.1, xlim=[-1,11], x00=0)
    ax5000.plot(u5000_bf.x, u5000_bf.y, 'x')
    ax5000.plot(u5000_id.x, u5000_id.y, 'o', mfc='none')
    i,j = get_loc_3d(gr.xh, gr.y, 5000., 3000.)
    ax5000.plot(u[:,j,i], gr.z-200, 'k-')
    ax5000.set_xlabel('u (m/s)')
    ax5000.grid()
    ax5000.set_xticks([0,5,10])
