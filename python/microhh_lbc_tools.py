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

import numpy as np
import xarray as xr
import pyproj


def get_cross_locations_for_lbcs(
        xstart_nest, ystart_nest,
        xend_nest, yend_nest,
        dx_parent, dy_parent,
        dx_child, dy_child,
        n_ghost, n_sponge):
    """
    Get `xz` and `yz` cross-section locations for nesting.
    NOTE: the resulting locations are a bit "spacious"
          (i.e. might contain more locations than strictly needed).
          This needs to be refined further...

    Parameters:
    -----------
    xstart_nest : float
        Start coordinate (m) child domain in x-direction.
    ystart_nest : float
        Start coordinate (m) child domain in y-direction.
    xend_nest : float
        End coordinate (m) child domain in x-direction.
    yend_nest : float
        End coordinate (m) child domain in y-direction.
    dx_parent : float
        Grid spacing parent domain in x-direction.
    dy_parent : float
        Grid spacing parent domain in y-direction.
    dx_child : float
        Grid spacing child domain in x-direction.
    dy_child : float
        Grid spacing child domain in y-direction.
    n_ghost : int
        Number of ghost cells.
    n_buffer : int
        Number of lateral sponge layer cells.
    """

    x0 = np.floor((xstart_nest - n_ghost  * dx_child) / dx_parent) * dx_parent - dx_parent/2
    x1 = np.ceil((xstart_nest + n_sponge * dx_child) / dx_parent) * dx_parent + dx_parent/2
    yzw = np.arange(x0, x1+1e-3, dx_parent)

    x0 = np.floor((xend_nest - n_sponge * dx_child) / dx_parent) * dx_parent - dx_parent/2
    x1 = np.ceil((xend_nest + n_ghost  * dx_child) / dx_parent) * dx_parent + dx_parent/2
    yze = np.arange(x0, x1+1e-3, dx_parent)

    yz = np.concatenate((yze, yzw))

    y0 = np.floor((ystart_nest - n_ghost  * dy_child) / dy_parent) * dy_parent - dx_parent/2
    y1 = np.ceil ((ystart_nest + n_sponge * dy_child) / dy_parent) * dy_parent + dx_parent/2
    xzs = np.arange(y0, y1+1e-3, dy_parent)

    y0 = np.floor((yend_nest - n_sponge * dy_child) / dy_parent) * dy_parent - dx_parent/2
    y1 = np.ceil ((yend_nest + n_ghost  * dy_child) / dy_parent) * dy_parent + dx_parent/2
    xzn = np.arange(y0, y1+1e-3, dy_parent)

    xz = np.concatenate((xzn, xzs))

    return xz, yz


def get_lbc_xr_dataset(
        fields,
        xsize, ysize,
        itot, jtot,
        z, zh,
        time,
        n_ghost, n_sponge,
        dtype=np.float64):
    """
    Create an Xarray Dataset that contains the correct
    variables / dimensions / coordinates needed to
    specify the lateral boundary conditions for an open-boundary
    simulation with MicroHH.

    Parameters:
    ----------
    fields : str
        List with prognostic field names.
    xsize : float
        Size of LES domain in x-direction.
    ysize : float
        Size of LES domain in y-direction.
    itot : int
        Number of grid points in x-direction.
    jtot : int
        Number of grid points in y-direction.
    z : np.ndarray
        Array with the full level grid point heights.
    zh : np.ndarray
        Array with the half level grid point heights.
    time : np.ndarray
        Array with LBC times.
    n_ghost : int
        Number of ghost cells.
    n_sponge : int
        Number of lateral sponge layer cells.
    dtype : np.float32 / np.float64
        Datatype used by MicroHH (SP=float32, DP=float64).
    """

    nt = time.size
    ktot = z.size

    dx = xsize/itot
    dy = ysize/jtot

    x = np.arange(dx/2, xsize, dx)
    xh = np.arange(0, xsize, dx)

    y = np.arange(dy/2, ysize, dy)
    yh = np.arange(0, ysize, dy)

    nlbc = n_ghost + n_sponge

    # Dimension sizes.
    dims = {
        'time': nt,
        'x': itot + 2*n_ghost,
        'xh': itot + 2*n_ghost,
        'xgw': nlbc,
        'xge': nlbc,
        'xhgw': nlbc + 1,
        'xhge': nlbc,
        'y': jtot + 2*n_ghost,
        'yh': jtot + 2*n_ghost,
        'ygn': nlbc,
        'ygs': nlbc,
        'yhgs': nlbc + 1,
        'yhgn': nlbc,
        'z': ktot}

    # Pad x,y dimensions with ghost cells.
    xp = np.zeros(x.size+2*n_ghost)
    xhp = np.zeros(x.size+2*n_ghost)

    yp = np.zeros(y.size+2*n_ghost)
    yhp = np.zeros(y.size+2*n_ghost)

    xp[n_ghost:-n_ghost] = x
    xhp[n_ghost:-n_ghost] = xh

    yp[n_ghost:-n_ghost] = y
    yhp[n_ghost:-n_ghost] = yh

    for i in range(n_ghost):
        xp[i] = x[0] - (n_ghost-i)*dx
        xhp[i] = xh[0] - (n_ghost-i)*dx

        yp[i] = y[0] - (n_ghost-i)*dy
        yhp[i] = yh[0] - (n_ghost-i)*dy
        
        xp[itot+n_ghost+i] = x[-1] + (i+1)*dx
        xhp[itot+n_ghost+i] = xh[-1] + (i+1)*dx

        yp[jtot+n_ghost+i] = y[-1] + (i+1)*dy
        yhp[jtot+n_ghost+i] = yh[-1] + (i+1)*dy
    
    # Define coordinates.
    coords = {
        'time': time,
        'x': xp,
        'xh': xhp,
        'y': yp,
        'yh': yhp,
        'z': z,
        'zh': zh,
        'xgw': xp[:nlbc],
        'xge': xp[itot+n_ghost-n_sponge:],
        'xhgw': xhp[:nlbc+1],
        'xhge': xhp[itot+n_ghost-n_sponge:],
        'ygs': yp[:nlbc],
        'ygn': yp[jtot+n_ghost-n_sponge:],
        'yhgs': yhp[:nlbc+1],
        'yhgn': yhp[jtot+n_ghost-n_sponge:]}

    # Create Xarray dataset.
    ds = xr.Dataset(coords=coords)

    def get_dim_size(dim_in):
        out = []
        for dim in dim_in:
            out.append(coords[dim].size)
        return out

    def add_var(name, dims):
        dim_size = get_dim_size(dims)
        ds[name] = (dims, np.empty(dim_size, dtype=dtype))

    for fld in fields:
        if fld not in ('u','v','w'):
            add_var(f'{fld}_west', ('time', 'z', 'y', 'xgw'))
            add_var(f'{fld}_east', ('time', 'z', 'y', 'xge'))
            add_var(f'{fld}_south', ('time', 'z', 'ygs', 'x'))
            add_var(f'{fld}_north', ('time', 'z', 'ygn', 'x'))

    if 'u' in fields:
        add_var('u_west', ('time', 'z', 'y', 'xhgw'))
        add_var('u_east', ('time', 'z', 'y', 'xhge'))
        add_var('u_south', ('time', 'z', 'ygs', 'xh'))
        add_var('u_north', ('time', 'z', 'ygn', 'xh'))

    if 'v' in fields:
        add_var('v_west', ('time', 'z', 'yh', 'xgw'))
        add_var('v_east', ('time', 'z', 'yh', 'xge'))
        add_var('v_south', ('time', 'z', 'yhgs', 'x'))
        add_var('v_north', ('time', 'z', 'yhgn', 'x'))

    if 'w' in fields:
        add_var('w_west', ('time', 'zh', 'y', 'xgw'))
        add_var('w_east', ('time', 'zh', 'y', 'xge'))
        add_var('w_south', ('time', 'zh', 'ygs', 'x'))
        add_var('w_north', ('time', 'zh', 'ygn', 'x'))

    return ds


def write_lbcs_as_binaries(
        lbc_ds,
        dtype,
        output_dir='.'):
    """
    Write all data variables from the Xarray Dataset `lbc_ds`
    in binary format as input to MicroHH.

    Parameters:
    -----------
    lbc_ds : Xarray.Dataset
        Input Dataset, created by `get_lbc_xr_dataset()`
    dtype : np.float32 || np.float64
        Datatype
    output_dir : str, optional
        Output directory of binary files.
    """

    for var in lbc_ds:
        lbc_ds[var].values.astype(dtype).tofile(f'{output_dir}/lbc_{var}.0000000')


class Domain:
    def __init__(
            self, 
            name,
            itot,
            jtot,
            dx,
            dy,
            work_path,
            end_time=None,
            start_offset=0,
            end_offset=0,
            i0_in_parent=None,
            j0_in_parent=None,
            parent=None,
            child=None,
            center_in_parent=False):
        """
        Domain and nesting specification.

        Parameters:
        -----------
        name : str
            Name of domain. Used to create `work_path/name` work directory.
        itot : int
            Number of grid points in x-direction.
        jtot : int
            Number of grid points in y-direction.
        dx : float
            Grid spacing in x-direction.
        dy : float
            Grid spacing in y-direction.
        work_path : str
            Work directory.
        end_time : int, optional
            End time of experiment (s). Mandatory for outer domain.
        start_offset : int, optional
            Time offset start of child domain relative to parent.
        end_offset : int, optional
            Time offset end of child domain relative to parent.
        i0_in_parent: int, optional
            Start index x-direction child in parent domain. Can be auto-calculated with `center_in_parent=True`.
        j0_in_parent: int, optional
            Start index y-direction child in parent domain. Can be auto-calculated with `center_in_parent=True`.
        parent : `Domain` instance, optional
            Parent of current domain.
        child : `Domain` instance, optional
            Child of current domain.
        center_in_parent : bool, optional
            Center current domain in parent domain.
        """

        self.name = name

        self.itot = itot
        self.jtot = jtot

        self.dx = dx
        self.dy = dy

        self.i0_in_parent = i0_in_parent
        self.j0_in_parent = j0_in_parent

        if parent is None:
            # Outer domain.
            if end_time is None:
                raise Exception('You must specify `end_time` for the outer domain!')

            self.start_offset = 0
            self.end_offset = 0

            self.start_time = 0
            self.end_time = end_time

        else:
            # Inner domains.
            if start_offset is None or end_offset is None:
                raise Exception('You must specify `start_offset` and `end_offset` for inner domains!')

            self.start_offset = start_offset
            self.end_offset = end_offset

            self.start_time = 0
            self.end_time = parent.end_time + end_offset - start_offset

        self.parent = parent
        self.child = child

        self.work_path = work_path
        self.work_dir = f'{work_path}/{name}'

        self.xsize = itot * dx
        self.ysize = jtot * dy

        if center_in_parent:
            self.center_in_parent()


    def center_in_parent(self):
        """
        Calculate `i0_in_parent` and `j0_in_parent`.
        """
        x0 = (self.parent.xsize - self.xsize) / 2.
        y0 = (self.parent.ysize - self.ysize) / 2.

        self.i0_in_parent = int(np.round(x0 / self.parent.dx))
        self.j0_in_parent = int(np.round(y0 / self.parent.dy))


class Projection:
    def __init__(
            self,
            xsize,
            ysize,
            itot,
            jtot,
            lon,
            lat,
            anchor='center',
            proj_str='+proj=utm +zone=31 +ellps=intl +towgs84=-87,-98,-121,0,0,0,0 +units=m +no_defs +type=crs'):
        """
        Projection class to transform LES (meters) to lon/lat (degrees) coordinate, and vice-versa.

        Several (lon,lat) pairs are defined:
        `(lon, lat)`: centers of each grid point (scalar location).
        `(lon_u, lat_u)`: middle-left edges of grid point (u location).
        `(lon_v, lat_v)`: lower-center edges of grid point (v location).
        `(lon_h, lat_h)`: lower-left edges of grid point (u,v location).

        Arguments:
            ----------
            xsize : float
                Domain size LES in x-direction (m).
            ysize : float
                Domain size LES in y-direction (m).
            itot : int
                Number of grid points in x-direction (-).
            jtot : int
                Number of grid points in y-direction (-).
            lon : float
                Longitude of LES domain. See `anchor` below (degrees).
            lat : float
                Latitude of LES domain. See `anchor` below (degrees).
            anchor : str, optional, default = 'center'
                Anchor point of (`lon,lat`), ∈ ('center', 'southwest')
            proj_str : str, optional, default = string for UTM31.
                Proj.4 / pyproj projection string.
        """

        self.proj_str = proj_str
        self.proj = pyproj.Proj(proj_str, preserve_units=True)

        self.xsize = xsize
        self.ysize = ysize

        self.itot = itot
        self.jtot = jtot

        self.dx = xsize / itot
        self.dy = ysize / jtot

        # Coordinates LES.
        self.x = np.arange(self.dx/2, self.xsize, self.dx)
        self.y = np.arange(self.dy/2, self.ysize, self.dy)

        self.xh = np.arange(0, self.xsize, self.dx)
        self.yh = np.arange(0, self.ysize, self.dy)

        if anchor == 'center':
            self.x_offset, self.y_offset = self.proj(lon, lat, inverse=False)
            self.x_offset -= self.xsize/2.
            self.y_offset -= self.ysize/2.

            self.central_lon = lon
            self.central_lat = lat

        elif anchor == 'southwest':
            self.x_offset, self.y_offset = self.proj(lon, lat, inverse=False)
            self.central_lon, self.central_lat = self.to_lonlat(self.xsize/2, self.ysize/2, meshgrid=False)

        else:
            raise Exception('Invalid anchor point domain!')

        # Coordinates at full and half levels.
        self.lon,   self.lat   = self.to_lonlat(self.x, self.y, meshgrid=True)
        self.lon_h, self.lat_h = self.to_lonlat(self.xh, self.yh, meshgrid=True)
        self.lon_u, self.lat_u = self.to_lonlat(self.xh, self.y, meshgrid=True)
        self.lon_v, self.lat_v = self.to_lonlat(self.x, self.yh, meshgrid=True)

        # Bounding box.
        x = np.array([0, self.xsize, self.xsize, 0, 0])
        y = np.array([0, 0, self.ysize, self.ysize, 0])

        self.bbox_lon, self.bbox_lat = self.to_lonlat(x, y, meshgrid=False)


    def to_lonlat(self, x, y, meshgrid=False):
        """
        Convert x/y (meters) LES coordinates to lon/lat (degrees).
        """
        if meshgrid:
            x, y = np.meshgrid(x, y)
        return self.proj(x+self.x_offset, y+self.y_offset, inverse=True)


    def to_xy(self, lon, lat):
        """
        Convert lon/lat (degrees) to LES (meters) coordinates.
        """
        x,y = self.proj(lon, lat, inverse=False)
        x -= self.x_offset
        y -= self.y_offset
        return x,y

