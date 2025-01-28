import numpy as np
import pyproj


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
