# use the Agg backend to have fast processing without X11
import matplotlib
matplotlib.use('Agg')
import netCDF4
from pylab import *
from scipy.interpolate import RectBivariateSpline

# parameters
kmin = 0
kmax = 1024

smin = -6.5
smax = 11.25

iterstart = 0
iterstep  = 1
iterend   = 1000

noborder = True

niter = (iterend-iterstart) / iterstep + 1

# transition data
stats = netCDF4.Dataset("slngrad.xz.00000.nc","r")
z     = stats.variables["z"][kmin:kmax]
x     = stats.variables["x"][:]
s     = stats.variables["slngrad"][-1,kmin:kmax,:]

# Black/Red/Yellow/White color scheme
# create color map
cdict = {'red':   (( 0.00, 0.00, 0.00),
                   ( 0.70, 0.50, 0.50),
                   ( 0.94, 1.00, 1.00),
                   ( 1.00, 1.00, 1.00)),
         'green': (( 0.00, 0.00, 0.00),
                   ( 0.70, 0.00, 0.00),
                   ( 0.94, 1.00, 1.00),
                   ( 1.00, 1.00, 1.00)),
         'blue':  (( 0.00, 0.00, 0.00),
                   ( 0.70, 0.00, 0.00),
                   ( 0.94, 0.40, 0.40),
                   ( 1.00, 1.00, 1.00))}
"""
cdict = {'red':   (( 0.00, 0    , 0    ),
                   ( 0.36, 0    , 0    ),
                   ( 0.59, 0    , 0    ),
                   ( 0.75, 1    , 1    ),
                   ( 0.96, 0.667, 0.667),
                   ( 1.00, 0.392, 0.392)),
         'green': (( 0.00, 0    , 0    ),
                   ( 0.36, 0.026, 0.026),
                   ( 0.59, 0.718, 0.718),
                   ( 0.75, 1    , 1    ),
                   ( 0.96, 0    , 0    ),
                   ( 1.00, 0    , 0    )),
         'blue':  (( 0.00, 0    , 0    ),
                   ( 0.36, 0.026, 0.026),
                   ( 0.59, 1    , 1    ),
                   ( 0.75, 0    , 0    ),
                   ( 0.96, 0    , 0    ),
                   ( 1.00, 0    , 0    ))}
"""
my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

# enable LaTeX plotting
#rc('font', family='serif')
#rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('font',**{'family':'serif','serif':['Palatino']})
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
ioff()

ext = [x.min(), x.max(), z.min(), z.max()]

for i in range(niter):
  s = stats.variables["slngrad"][i,kmin:kmax,:]
  iter = iterstart + i*iterstep
  print("Processing iteration: ", iter)

  close('all')
  fig1 = figure(figsize=(20.48,10.24))
  if(noborder):
    fig1.subplots_adjust(0.000,0.000,1.0,1.0,0.,0.)
  else:
    fig1.subplots_adjust(0.05,0.05,0.97,0.97,0,0)
  plt = imshow(s, origin='lower', extent=ext, cmap=my_cmap, vmin=smin, vmax=smax)
  #colorbar(pad=0.02)
  xlabel(r'$x$')
  ylabel(r'$z$')
  axis(ext)
  if(noborder):
    plt.axes.set_frame_on(False)
    plt.axes.get_xaxis().set_visible(False)
    plt.axes.get_yaxis().set_visible(False)
    plt.axes.get_yaxis().set_visible(False)
  savefig('figs/slngrad.{:07d}.png'.format(iter), facecolor='black', edgecolor='none', dpi=100)
  del(s)

