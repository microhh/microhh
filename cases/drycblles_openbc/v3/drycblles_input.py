import matplotlib.pyplot as pl
import numpy as np
import netCDF4 as nc
import xarray as xr
from numba import jit, prange
import sys
import glob
from datetime import datetime
import asyncio

from lbc_input import lbc_input
import microhh_tools as mht

pl.close('all')

domain = sys.argv[1]

dtype = np.float64

swadvec = '2'   # needed for correct # gcs.

ini = mht.Read_namelist('drycblles.ini.base')

"""
Grid & nesting settings.
"""

# Grid settings outer domain.
#itot = 64
#jtot = 64
#ktot = 64

itot = 512
jtot = 512
ktot = 256

xsize = 6400
ysize = 6400
zsize = 3200

dx = xsize / itot
dy = ysize / jtot
dz = zsize / ktot

# Nest settings.
refinement_fac = 3

# Start index in parent domain.
i0_nest = 12
j0_nest = 16

# Size domain in parent coordinates!
# The nest itself has `refinement_fac` times as many grid points.
itot_nest = 48
jtot_nest = 32

# Number of lateral buffer points.
nbuffer = 5

# Number of ghost cells.
if swadvec == '2':
    nghost = 1
elif swadvec == '2i4':
    nghost = 2
elif swadvec == '2i5':
    nghost = 3

xstart_nest = i0_nest*dx
ystart_nest = j0_nest*dy

xend_nest = xstart_nest + itot_nest * dx
yend_nest = ystart_nest + jtot_nest * dy

xsize_nest = xend_nest - xstart_nest
ysize_nest = yend_nest - ystart_nest

dx_nest = xsize_nest / (itot_nest * refinement_fac)
dy_nest = ysize_nest / (jtot_nest * refinement_fac)

"""
Define initial fields/profiles.
"""
dthetadz = 0.003

z  = np.arange(0.5*dz, zsize, dz)
zh = np.arange(0, zsize, dz)

u  = np.zeros(np.size(z))
v  = np.zeros(np.size(z))
th = np.zeros(np.size(z))

for k in range(ktot):
    th[k] = 300. + dthetadz*z[k]

"""
Write case_input.nc
"""
nc_file = nc.Dataset("drycblles_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", ktot)
nc_z  = nc_file.createVariable("z" , dtype, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_u  = nc_group_init.createVariable("u" , dtype, ("z"))
nc_v  = nc_group_init.createVariable("v" , dtype, ("z"))
nc_th = nc_group_init.createVariable("th", dtype, ("z"))

nc_z [:] = z [:]
nc_u [:] = u [:] + 2
nc_v [:] = v [:] + 1
nc_th[:] = th[:]

nc_file.close()

"""
Update .ini file.
"""
ini['advec']['swadvec'] = swadvec

if domain == 'outer':

    x0 = np.floor((xstart_nest - nghost  * dx_nest) / dx) * dx
    x1 = np.ceil ((xstart_nest + nbuffer * dx_nest) / dx) * dx
    yzw = np.arange(x0, x1+1e-3, dx)

    x0 = np.floor((xend_nest - nbuffer * dx_nest) / dx) * dx
    x1 = np.ceil ((xend_nest + nghost  * dx_nest) / dx) * dx
    yze = np.arange(x0, x1+1e-3, dx)

    yz = np.concatenate((yzw, yze))

    y0 = np.floor((ystart_nest - nghost  * dy_nest) / dy) * dy
    y1 = np.ceil ((ystart_nest + nbuffer * dy_nest) / dy) * dy
    xzs = np.arange(y0, y1+1e-3, dy)

    y0 = np.floor((yend_nest - nbuffer * dy_nest) / dy) * dy
    y1 = np.ceil ((yend_nest + nghost  * dy_nest) / dy) * dy
    xzn = np.arange(y0, y1+1e-3, dy)

    xz = np.concatenate((xzs, xzn))

    ini['grid']['itot'] = itot
    ini['grid']['jtot'] = jtot
    ini['grid']['ktot'] = ktot

    ini['grid']['xsize'] = xsize
    ini['grid']['ysize'] = ysize
    ini['grid']['zsize'] = zsize

    ini['cross']['yz'] = list(yz)
    ini['cross']['xz'] = list(xz)

    ini['pres']['sw_openbc'] = False
    ini['boundary_lateral']['sw_openbc'] = False

elif domain == 'inner':

    ini['grid']['itot'] = itot_nest * refinement_fac
    ini['grid']['jtot'] = jtot_nest * refinement_fac
    ini['grid']['ktot'] = ktot

    ini['grid']['xsize'] = xsize_nest
    ini['grid']['ysize'] = ysize_nest
    ini['grid']['zsize'] = zsize

    ini['cross']['yz'] = (xend_nest+xstart_nest)/2
    ini['cross']['xz'] = (yend_nest+ystart_nest)/2

    ini['pres']['sw_openbc'] = True
    ini['boundary_lateral']['sw_openbc'] = True

ini.save('drycblles.ini', allow_overwrite=True)


"""
Create lateral boundaries
"""
class Timer:
    def __init__(self):
        self.tstart = datetime.now()
        self.elapsed = 0

    def pause(self):
        self.elapsed += (datetime.now() - self.tstart).total_seconds()

    def restart(self):
        self.tstart = datetime.now()

    def stop(self):
        self.elapsed += (datetime.now() - self.tstart).total_seconds()
        print(f'Elapsed: {self.elapsed}')


if domain == 'inner':

    fields = ['u', 'v', 'w', 'th', 's']

    t = Timer()

    # Read cross-sections.
    xz = {}
    yz = {}
    for fld in fields:
        xz[fld] = xr.open_dataset(f'outer/{fld}.xz.nc', decode_times=False)
        yz[fld] = xr.open_dataset(f'outer/{fld}.yz.nc', decode_times=False)
    time = xz[list(xz.keys())[0]].time.values

    t.pause()

    # Define nest grid.
    x = np.arange(dx_nest/2, xsize_nest, dx_nest)
    xh = np.arange(0, xsize_nest, dx_nest)
    y = np.arange(dy_nest/2, ysize_nest, dy_nest)
    yh = np.arange(0, ysize_nest, dy_nest)

    lbc = lbc_input(fields, x, y, z, xh, yh, zh, time, nghost, nbuffer, x_offset=xstart_nest, y_offset=ystart_nest)

    t.restart()

    for loc in ['west', 'east', 'north', 'south']:
        for fld in fields:
            # Short cuts.
            lbc_in = lbc[f'{fld}_{loc}']
            dims = lbc_in.dims

            # Dimensions in LBC file.
            xloc, yloc = dims[3], dims[2]

            # Dimensions in cross-section.
            xloc_in = 'xh' if 'xh' in xloc else 'x'
            yloc_in = 'yh' if 'yh' in yloc else 'y'

            # Switch between yz and xz crosses.
            cc = yz if loc in ['west','east'] else xz

            # Interpolate!
            ip = cc[fld].interp({yloc_in: lbc[yloc], xloc_in: lbc[xloc]}, method='nearest')

            # Check if interpolation was success.
            if np.any(np.isnan(ip[fld].values)):
                raise Exception('Interpolated BCs contain NaNs!')

            lbc_in[:] = ip[fld].values

    t.stop()

    #lbc.to_netcdf('lbc_old.nc')

    # Save as binary input files.
    for fld in fields:
        for loc in ['west', 'east', 'north', 'south']:
            lbc_in = lbc[f'{fld}_{loc}']
            lbc_in.values.astype(dtype).tofile('lbc_{}_{}.0000000'.format(fld, loc))

    """
    NEW NEW NEW!!!
    """
    lbc_new = lbc_input(fields, x, y, z, xh, yh, zh, time, nghost, nbuffer, x_offset=xstart_nest, y_offset=ystart_nest)

    grid = mht.Read_grid(
            itot, jtot, ktot, filename='outer/grid.0000000')

    def run_async(f):
        """
        Decorator to run processes asynchronous with `asyncio`.
        """
        def wrapped(*args, **kwargs):
            return asyncio.get_event_loop().run_in_executor(None, f, *args, **kwargs)
        return wrapped


    @run_async
    def read_and_interp(arr_out, var, edge, loc, t, time, index, nn_i, nn_j, itot, jtot, ktot, parent_dir, dtype):
        """
        Read cross-section directly from binary, and interpolate to LBC location.
        """

        if edge in ('west', 'east'):
            data = np.zeros((ktot, jtot, index.size), dtype)
            for i,ii in enumerate(index):
                raw = np.fromfile(f'{parent_dir}/{var}.yz.{loc}.{ii:05d}.{time:07d}', dtype=dtype)
                data[:,:,i] = raw.reshape((ktot, jtot))

        else:
            data = np.zeros((ktot, index.size, itot), dtype)
            for j,jj in enumerate(index):
                raw = np.fromfile(f'{parent_dir}/{var}.xz.{loc}.{jj:05d}.{time:07d}', dtype=dtype)
                data[:,j,:] = raw.reshape((ktot, itot))

        _interpolate(
                arr_out[t,:,:,:],
                data,
                nn_i,
                nn_j)


    @jit(nopython=True, fastmath=True, nogil=True, parallel=True)
    def _interpolate(arr_out, arr_in, nn_i, nn_j):
        """
        Fast NN interpolation using Numba.
        """

        itot = nn_i.size
        jtot = nn_j.size
        ktot = arr_out.shape[0]

        for k in prange(ktot):
            for j in prange(jtot):
                for i in prange(itot):
                    arr_out[k, j, i] = arr_in[k, nn_j[j], nn_i[i]]


    class LBC_interpolation:
        def __init__(self, lbc_ds, var, parent_dir, dims_x, dims_y, dims_z, t0, t1, dt, dtype):
            """
            Help class to interpolate LBCs directly from binary cross-section files.
            """

            self.lbc_ds = lbc_ds
            self.var = var
            self.parent_dir = parent_dir

            self.dims_x = dims_x
            self.dims_y = dims_y

            if var == 'w':
                self.dims_z = dims_z[:-1]
            else:
                self.dims_z = dims_z

            self.itot = self.dims_x.size
            self.jtot = self.dims_y.size
            self.ktot = self.dims_z.size

            self.t0 = t0
            self.t1 = t1
            self.dt = dt

            self.dtype = dtype

            # NN indexes:
            self.nn_i = {}
            self.nn_j = {}

            # Locations of cross-slices.
            self.index = {}

            # Data buffer for single time step.
            self.data = {}

            self.edges = ('west', 'east', 'south', 'north')

            # Location indicator in cross-section file name:
            if var == 'u':
                self.loc = '100'
            elif var == 'v':
                self.loc = '010'
            elif var == 'w':
                self.loc = '001'
            else:
                self.loc = '000'

            for edge in self.edges:
                var_name = f'{var}_{edge}'

                dim_name_x = lbc_ds[var_name].dims[3]
                dim_name_y = lbc_ds[var_name].dims[2]

                dim_size_x = lbc_ds[var_name].shape[3]
                dim_size_y = lbc_ds[var_name].shape[2]

                # Switch between west or east, or south or north.
                which = 'lower' if edge in ('west', 'south') else 'upper'

                if edge in ('west', 'east'):
                    cross_glob = f'{parent_dir}/{var}.yz.{self.loc}.*.{t0:07d}'
                    i_cross, x_cross = self._get_index_and_dims(cross_glob, dims_x, which)

                    self.nn_i[edge] = self._get_NN_indexes(lbc_ds[var_name][dim_name_x].values, x_cross)
                    self.nn_j[edge] = self._get_NN_indexes(lbc_ds[var_name][dim_name_y].values, dims_y)

                    self.index[edge] = i_cross
                    self.data[edge] = np.zeros((ktot, jtot, i_cross.size), dtype)

                else:
                    cross_glob = f'{parent_dir}/{var}.xz.{self.loc}.*.{t0:07d}'
                    j_cross, y_cross = self._get_index_and_dims(cross_glob, dims_y, which)

                    self.nn_i[edge] = self._get_NN_indexes(lbc_ds[var_name][dim_name_x].values, dims_x)
                    self.nn_j[edge] = self._get_NN_indexes(lbc_ds[var_name][dim_name_y].values, y_cross)

                    self.index[edge] = j_cross
                    self.data[edge] = np.zeros((ktot, j_cross.size, itot), dtype)

            # Interpolate!
            self.interpolate()


        def interpolate(self):
            """
            Interpolate all edges and time steps.
            """

            #for t,time in enumerate(range(self.t0, self.t1+1, self.dt)):
            #    for edge in self.edges:

            #        # Read new cross-slices directly from binaries:
            #        data = self.read_cross(edge, time)

            #        # Interpolate to LBC locations:
            #        self._interpolate(
            #                self.lbc_ds[f'{self.var}_{edge}'].values,
            #                self.data[edge],
            #                self.nn_i[edge],
            #                self.nn_j[edge],
            #                t)

            calls = []
            for t,time in enumerate(range(self.t0, self.t1+1, self.dt)):
                for edge in self.edges:

                    calls.append(
                        read_and_interp(
                                self.lbc_ds[f'{self.var}_{edge}'].values,
                                self.var, edge, self.loc,
                                t, time,
                                self.index[edge],
                                self.nn_i[edge], self.nn_j[edge],
                                self.itot, self.jtot, self.ktot,
                                self.parent_dir, self.dtype))

            loop = asyncio.get_event_loop()
            looper = asyncio.gather(*calls)
            results = loop.run_until_complete(looper)


        @staticmethod
        @jit(nopython=True, fastmath=True, nogil=True, parallel=True)
        def _interpolate(arr_out, arr_in, nn_i, nn_j, t):
            """
            Fast NN interpolation using Numba.
            """

            itot = nn_i.size
            jtot = nn_j.size
            ktot = arr_out.shape[1]

            for k in prange(ktot):
                for j in prange(jtot):
                    for i in prange(itot):
                        arr_out[t,k,j,i] = arr_in[k, nn_j[j], nn_i[i]]


        def read_cross(self, edge, time):
            """
            Read cross-sections from binaries for single time step.
            """

            if edge in ('west', 'east'):
                for i,ii in enumerate(self.index[edge]):
                    raw = np.fromfile(f'{self.parent_dir}/{self.var}.yz.{self.loc}.{ii:05d}.{time:07d}', dtype=self.dtype)
                    self.data[edge][:, :, i] = raw.reshape((self.ktot, self.jtot))

            else:
                for j,jj in enumerate(self.index[edge]):
                    raw = np.fromfile(f'{self.parent_dir}/{self.var}.xz.{self.loc}.{jj:05d}.{time:07d}', dtype=self.dtype)
                    self.data[edge][:, j, :] = raw.reshape((self.ktot, self.itot))





        def _get_index_and_dims(self, cross_search, dims_in, which):
            """
            Get indexes (-) and coordinates (m) of available cross-section planes.
            For a `xz` cross, get the y indexes/coordinates, and
            for a `yz` cross, get the x indexes/coordinates.
            """

            crosses = glob.glob(cross_search)
            crosses.sort()

            if len(crosses) == 0:
                raise Exception('Failed to find any crosses!')

            # Get indices (-) of cross planes.
            index = np.array([int(f.split('.')[-2]) for f in crosses])

            # Select west/south (`lower`) or east/north (`upper`) part.
            split = int(np.where((index[1:] - index[:-1]) > 1)[0]) + 1

            if which == 'lower':
                index = index[:split]
            else:
                index = index[split:]

            # Get coordinates (m) of cross planes.
            dims = dims_in[index]

            return index, dims


        def _get_NN_indexes(self, dim_out, dim_in):
            """
            Get nearest neighbour index of coordinates `dim_out` in `dim_in`
            """

            index_out = np.zeros(dim_out.size, np.uint16)

            # Allow NN interpolations to be half dx or dy out of bounds.
            margin = (dim_in[1] - dim_in[0]) / 2.

            # Find NN indexes.
            for i in range(index_out.size):
                if dim_out[i] < (dim_in[0] - margin) or dim_out[i] > (dim_out[-1] + margin):
                    raise Exception('NN index out of bounds!')

                index_out[i] = np.argmin(np.abs(dim_out[i] - dim_in))

            return index_out



    t0 = 0
    t1 = 1200
    dt = 60

    t = Timer()

    lbc_u = LBC_interpolation(lbc_new, 'u',  'outer/', grid.dim['xh'], grid.dim['y' ], grid.dim['z' ], t0, t1, dt, dtype)
    lbc_v = LBC_interpolation(lbc_new, 'v',  'outer/', grid.dim['x' ], grid.dim['yh'], grid.dim['z' ], t0, t1, dt, dtype)
    lbc_w = LBC_interpolation(lbc_new, 'w',  'outer/', grid.dim['x' ], grid.dim['y' ], grid.dim['zh'], t0, t1, dt, dtype)
    lbc_s = LBC_interpolation(lbc_new, 'th', 'outer/', grid.dim['x' ], grid.dim['y' ], grid.dim['z' ], t0, t1, dt, dtype)
    lbc_s = LBC_interpolation(lbc_new, 's',  'outer/', grid.dim['x' ], grid.dim['y' ], grid.dim['z' ], t0, t1, dt, dtype)

    t.stop()

    # Check old vs new method.
    print('Diffs old vs new:')
    for v in ['u', 'v', 'w', 'th', 's']:
        for edge in ['west', 'east', 'south', 'north']:
            print(v, edge, float(np.abs(lbc_new[f'{v}_{edge}'] - lbc[f'{v}_{edge}']).max()))

