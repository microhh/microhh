#
#  MicroHH
#  Copyright (c) 2011-2023 Chiel van Heerwaarden
#  Copyright (c) 2011-2023 Thijs Heus
#  Copyright (c) 2014-2023 Bart van Stratum
#
#  This file is part of MicroHH
#
#  MicroHH is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MicroHH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
#

import netCDF4
from pylab import *
from scipy.interpolate import RectBivariateSpline

# parameters
kmin = 0
kmax = 900

smin = -8.
smax = 13.

iterstart = 100
iterstep = 100
iterend = 60000

noborder = True

niter = (iterend - iterstart) / iterstep + 1
ioff()

# transition data
stats = netCDF4.Dataset("slngrad.xzcross.nc", "r")
z = stats.variables["z"][kmin:kmax]
x = stats.variables["x"][:]
s = stats.variables["slngrad"][-1, kmin:kmax, :]

# color scheme black/blue/yellow/red
cdict = {'red': ((0.00, 0, 0),
                 (0.36, 0, 0),
                 (0.59, 0, 0),
                 (0.75, 1, 1),
                 (0.96, 0.667, 0.667),
                 (1.00, 0.392, 0.392)),
         'green': ((0.00, 0, 0),
                   (0.36, 0.026, 0.026),
                   (0.59, 0.718, 0.718),
                   (0.75, 1, 1),
                   (0.96, 0, 0),
                   (1.00, 0, 0)),
         'blue': ((0.00, 0, 0),
                  (0.36, 0.026, 0.026),
                  (0.59, 1, 1),
                  (0.75, 0, 0),
                  (0.96, 0, 0),
                  (1.00, 0, 0))}
my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 256)

# enable LaTeX plotting
#rc('font', family='serif')
#rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('font', **{'family': 'serif', 'serif': ['Palatino']})
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

# create interpolation function
fint = RectBivariateSpline(z, x, s)
del(s)
# create interpolated value
zint = linspace(z.min(), z.max(), 1680)
sint = fint(zint, x)

ext = [x.min(), x.max(), zint.min(), zint.max()]

for i in range(niter):
    s = stats.variables["slngrad"][i, kmin:kmax, :]
    iter = iterstart + i * iterstep
    print("Processing iteration: ", iter)
    fint = RectBivariateSpline(z, x, s)
    sint = fint(zint, x)
    del(s)

    close('all')
    fig1 = figure(
        figsize=(19.2, 19.2 * (zint.max() - zint.min()) / (x.max() - x.min())))
    if(noborder):
        fig1.subplots_adjust(0.000, 0.000, 1.0, 1.0, 0., 0.)
    else:
        fig1.subplots_adjust(0.05, 0.05, 0.97, 0.97, 0, 0)
    plt = imshow(
        sint,
        origin='lower',
        extent=ext,
        cmap=my_cmap,
        vmin=smin,
        vmax=smax)
    # colorbar(pad=0.02)
    xlabel(r'$x$')
    ylabel(r'$z$')
    axis(ext)
    if(noborder):
        plt.axes.set_frame_on(False)
        plt.axes.get_xaxis().set_visible(False)
        plt.axes.get_yaxis().set_visible(False)
        plt.axes.get_yaxis().set_visible(False)
    savefig(
        'figs/slngrad.{:07d}.png'.format(iter),
        facecolor='black',
        edgecolor='none',
        dpi=100)
    del(fint, sint)
