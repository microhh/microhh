# use the Agg backend to have fast processing without X11
import matplotlib
matplotlib.use('Agg')
import netCDF4
from pylab import *
from scipy.interpolate import RectBivariateSpline

# parameters
kmin = 0
kmax = 1024

smin = 0.
smax = 1.

iterstart = 0
iterstep  = 1
iterend   = 1000

noborder = True

niter = (iterend-iterstart) / iterstep + 1

# transition data
stats = netCDF4.Dataset("s.xz.00000.nc","r")
z     = stats.variables["z"][kmin:kmax]
x     = stats.variables["x"][:]
s     = stats.variables["s"][-1,kmin:kmax,:]

# Blue to white color scheme
# create color map
cdict = {'red':   (( 0.00, 0.10, 0.10),
                   ( 1.00, 1.00, 1.00)),
         'green': (( 0.00, 0.20, 0.20),
                   ( 1.00, 1.00, 1.00)),
         'blue':  (( 0.00, 0.30, 0.30),
                   ( 1.00, 1.00, 1.00))}
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
  s = stats.variables["s"][i,kmin:kmax,:]
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
  savefig('figs/s.{:07d}.png'.format(iter), facecolor='black', edgecolor='none', dpi=100)
  del(s)

