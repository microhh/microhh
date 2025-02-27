import numpy as np
from numba import njit

@njit
def _pow2(x):
    return x*x


@njit
def _add_source_kernel(
        field,
        t,
        strength,
        x0,
        y0,
        z0,
        sigma_x,
        sigma_y,
        sigma_z,
        x,
        y,
        z,
        dx,
        dy,
        dz,
        rho_ref):

    """
    Fast Numba kernel to add source as Gaussian blob to field,
    and normalize it to get a net `strength` as sum over blob.
    """

    nt, ktot, jtot, itot = field.shape

    raw_sum = 0.
    for k in range(ktot):

        # Limit vertical extent to +/- 4*sigma_z. This allows us
        # to clip the input fields, reducing their size.
        if (z[k] > z0-4*sigma_z) and (z[k] < z0+4*sigma_z):

            for j in range(jtot):
                for i in range(itot):

                    blob_norm = np.exp(
                            - _pow2(x[i]-x0)/_pow2(sigma_x)
                            - _pow2(y[j]-y0)/_pow2(sigma_y)
                            - _pow2(z[k]-z0)/_pow2(sigma_z))

                    raw_sum += blob_norm * dx * dy * dz[k] * rho_ref[k]

    scaling = strength / raw_sum

    for k in range(ktot):

        # Limit vertical extent to +/- 4*sigma_z. This allows us
        # to clip the input fields, reducing their size.
        if (z[k] > z0-4*sigma_z) and (z[k] < z0+4*sigma_z):

            for j in range(jtot):
                for i in range(itot):

                    field[t,k,j,i] += scaling * np.exp(
                            - _pow2(x[i]-x0)/_pow2(sigma_x)
                            - _pow2(y[j]-y0)/_pow2(sigma_y)
                            - _pow2(z[k]-z0)/_pow2(sigma_z))


class Emission_input:
    """
    """
    def __init__(
            self,
            fields,
            times,
            x,
            y,
            z,
            dz,
            rho_ref,
            TF=np.float64):
        """
        Help class to define point source emissions.
        Emissions can be added to a single grid point, or as Gaussian "blobs".

        After defining an `emiss = Emission_input(...)` object, emissions can be added as:
            `emiss.add_gaussian(field='s1', strength=1, time=0, x0=400, y0=800, z0=50, sigma_x=50, sigma_y=50, sigma_z=25)`
        or:
            `emiss.add_point(field='s1', strength=1, time=0, x0=400, y0=800, z0=50)`

        Parameters
        ----------
        fields : list
            List with scalar fields that have an emission.
        times : np.ndarray, shape (1,)
            Array with output times (s).
        x : np.ndarray, shape (1,)
            Array with x-coordinates grid (m).
        y : np.ndarray, shape (1,)
            Array with y-coordinates grid (m).
        z : np.ndarray (shape 1,)
            Array with z-coordinates grid (m).
        dz : np.ndarray shape (1,)
            Array with full level vertical grid spacing (m).
        rho_ref : np.ndarray shape (1,)
            Array with base state density (kg m-3).
        TF : np.float32 or np.float64
            Datatype used by MicroHH.
        """

        self.times = times
        self.fields = fields
        self.rho_ref = rho_ref

        self.nt = times.size

        self.x = x
        self.y = y
        self.z = z

        self.itot = x.size
        self.jtot = y.size
        self.ktot = z.size

        self.dx = x[1] - x[0]
        self.dy = y[1] - y[0]
        self.dz = dz

        # Create 3D emission fields.
        self.data = {}
        for field in self.fields:
            self.data[field] = np.zeros((self.nt, self.ktot, self.jtot, self.itot), dtype=TF)

        self.is_clipped = False


    def get_index(self, time):
        """
        Get index of `time` in global time array.
        """
        t = np.argmin(np.abs(self.times - time))

        if np.abs(self.times[t] - time) > 1e-6:
            raise Exception(f'Can not find time {time} in {self.times}.')

        return t


    def add_gaussian(self, field, strength, time, x0, y0, z0, sigma_x, sigma_y, sigma_z):
        """
        Add single point source emission, spread out over a Gaussian blob defined by `sigma_xyz`.

        Parameters
        ----------
        field : string
            Emission field name.
        strength : float
            Emission strength (kg s-1).
        time : int
            Emission time (s).
        x0 : float
            Center (x) of Gaussian blob (m).
        y0 : float
            Center (y) of Gaussian blob (m).
        z0 : float
            Center (z) of Gaussian blob (m) .
        sigma_x : float
            Std.dev of Gaussian blob in x-direction (m).
        sigma_y : float
            Std.dev of Gaussian blob in y-direction (m).
        sigma_z : float
            Std.dev of Gaussian blob in z-direction (m).
        """

        t = self.get_index(time)

        _add_source_kernel(
                self.data[field],
                t,
                strength,
                x0,
                y0,
                z0,
                sigma_x,
                sigma_y,
                sigma_z,
                self.x,
                self.y,
                self.z,
                self.dx,
                self.dy,
                self.dz,
                self.rho_ref)


    def add_point(self, field, strength, time, x0, y0, z0):
        """
        Add single point source emission, to single grid point.

        Parameters
        ----------
        field : string
            Emission field name.
        strength : float
            Emission strength (kg s-1).
        time : int
            Emission time (s).
        x0 : float
            Center (x) of emission (m).
        y0 : float
            Center (y) of emission (m).
        z0 : float
            Center (z) of emission (m) .
        """

        t = self.get_index(time)

        i = np.abs(self.x - x0).argmin()
        j = np.abs(self.y - y0).argmin()
        k = np.abs(self.z - z0).argmin()

        volume = self.dx * self.dy * self.dz[k]
        self.data[field][t,k,j,i] += strength / (volume * self.rho_ref[k])


    def clip(self):
        """
        Clip 3D fields to required vertical extent.
        This automatically determines the highest emission height.
        """
        self.kmax = 0

        # Find max height over all fields.
        for field, emission in self.data.items():
            # Take sum over time and x,y dimensions to get vertical profile.
            emiss_sum = emission[:,:,:,:].sum(axis=(0,2,3))

#            pl.figure()
#            pl.plot(emiss_sum, self.z)

            # Max height where total emission is positive.
            emiss_pos = np.where(emiss_sum > 0)[0]
            self.kmax = max(self.kmax, emiss_pos[-1]+1)

        # Clip fields.
        for field, emission in self.data.items():
            self.data[field] = emission[:,:self.kmax,:,:]

        print(f'Max emission height = {self.z[self.kmax]} m., ktot={self.kmax}')

        self.is_clipped = True


    def to_binary(self, path):
        """
        Save all fields in binary format for MicroHH.
        """

        if not self.is_clipped:
            print('WARNING: saving unclipped fields!')

        for name, fld in self.data.items():
            for t,time in enumerate(self.times):
                fld[t,:].tofile(f'{path}/{name}_emission.{time:07d}')



if __name__ == '__main__':
    """
    Just for testing...
    """
    import matplotlib.pyplot as pl
    pl.close('all')

    xsize = 3200
    ysize = 1600
    zsize = 2400

    itot = 64
    jtot = 32
    ktot = 96

    dx = xsize / itot
    dy = ysize / jtot
    dz0 = zsize / ktot

    x = np.arange(dx/2, xsize, dx)
    y = np.arange(dy/2, ysize, dy)
    z = np.arange(dz0/2, zsize, dz0)
    zh = np.arange(0, zsize+0.1, dz0)

    dz = zh[1:] - zh[:-1]
    rho_ref = np.ones(ktot)

    times = np.array([0])

    fields = ['s1']
    emiss = Emission_input(fields, times, x, y, z, dz, rho_ref)

    emiss.add_gaussian(field='s1', strength=1, time=0,    x0=400,  y0=800, z0=100, sigma_x=25, sigma_y=25, sigma_z=25)
    #emiss.add_gaussian(field='s1', strength=1, time=3600, x0=800,  y0=800, z0=100, sigma_x=25, sigma_y=25, sigma_z=25)
    #emiss.add_gaussian(field='s1', strength=2, time=7200, x0=1600, y0=800, z0=200, sigma_x=25, sigma_y=25, sigma_z=25)
    
    #emiss.add_point(field='s1', strength=1, time=0,    x0=2000, y0=800, z0=200)
    #emiss.add_point(field='s1', strength=2, time=3600, x0=2000, y0=800, z0=200)
    #emiss.add_point(field='s1', strength=3, time=7200, x0=2000, y0=800, z0=200)

    print(emiss.data['s1'].sum() * dx * dy * dz0)

    pl.figure()
    pl.pcolormesh(x, z, emiss.data['s1'][0, :, :, :].sum(axis=1))
    pl.colorbar()

    emiss.clip()
    emiss.to_binary('.')
