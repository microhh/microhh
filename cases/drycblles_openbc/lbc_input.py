import numpy as np
import netCDF4 as nc4


class LBC_input:
    """
    Help class to create/write lateral boundary conditions
    in the MicroHH-specific format
    """
    def __init__(self, fields, itot, jtot, ktot, time, dtype=np.float64):
        self.itot = itot
        self.jtot = jtot
        self.ktot = ktot
        self.time = time
        self.fields = fields
        self.dtype = dtype
        self.nt = time.size

        # Momentum
        if 'u' in fields:
            self.u_west = np.zeros((self.nt, ktot, jtot), dtype)
            self.u_east = np.zeros((self.nt, ktot, jtot), dtype)
            self.u_south = np.zeros((self.nt, ktot, itot), dtype)
            self.u_north = np.zeros((self.nt, ktot, itot), dtype)

        if 'v' in fields:
            self.v_west = np.zeros((self.nt, ktot, jtot), dtype)
            self.v_east = np.zeros((self.nt, ktot, jtot), dtype)
            self.v_south = np.zeros((self.nt, ktot, itot), dtype)
            self.v_north = np.zeros((self.nt, ktot, itot), dtype)

        #if 'w' in fields:
        #    self.w_west = np.zeros((self.nt, ktot+1, jtot), dtype)
        #    self.w_east = np.zeros((self.nt, ktot+1, jtot), dtype)
        #    self.w_south = np.zeros((self.nt, ktot+1, itot), dtype)
        #    self.w_north = np.zeros((self.nt, ktot+1, itot), dtype)

        for f in fields:
            if f not in ['u', 'v', 'w']:
                setattr(self, '{}_west'.format(f), np.zeros((self.nt, ktot, jtot), dtype))
                setattr(self, '{}_east'.format(f), np.zeros((self.nt, ktot, jtot), dtype))
                setattr(self, '{}_south'.format(f), np.zeros((self.nt, ktot, itot), dtype))
                setattr(self, '{}_north'.format(f), np.zeros((self.nt, ktot, itot), dtype))


    def to_netcdf(self, case_name):
        """
        Write LBCs in MicroHH-specific format.
        """

        nc = nc4.Dataset(
                '{}_lbc_input.nc'.format(case_name), mode='w', datamodel='NETCDF4', clobber=True)
        
        nc.createDimension('time', self.nt)

        nc.createDimension('x', self.itot)
        nc.createDimension('y', self.jtot)
        nc.createDimension('z', self.ktot)

        nc.createDimension('xh', self.itot)
        nc.createDimension('yh', self.jtot)
        nc.createDimension('zh', self.ktot)

        nc_t = nc.createVariable('time',  self.dtype, ('time'))
        nc_t[:] = self.time

        def add_variable(name, dims, data):
            var = nc.createVariable(name, self.dtype, dims)
            var[:] = data

        if 'u' in self.fields:
            add_variable('u_west', ('time', 'z', 'y'), self.u_west)
            add_variable('u_east', ('time', 'z', 'y'), self.u_east)
            add_variable('u_south', ('time', 'z', 'xh'), self.u_south)
            add_variable('u_north', ('time', 'z', 'xh'), self.u_north)

        if 'v' in self.fields:
            add_variable('v_west', ('time', 'z', 'yh'), self.v_west)
            add_variable('v_east', ('time', 'z', 'yh'), self.v_east)
            add_variable('v_south', ('time', 'z', 'x'), self.v_south)
            add_variable('v_north', ('time', 'z', 'x'), self.v_north)

        #if 'w' in self.fields:
        #    add_variable('w_west', ('time', 'zh', 'y'), self.w_west)
        #    add_variable('w_east', ('time', 'zh', 'y'), self.w_east)
        #    add_variable('w_south', ('time', 'zh', 'x'), self.w_south)
        #    add_variable('w_north', ('time', 'zh', 'x'), self.w_north)

        for f in self.fields:
            if f not in ['u', 'v', 'w']:
                add_variable('{}_west'.format(f), ('time', 'z', 'y'),  getattr(self, '{}_west'.format(f)))
                add_variable('{}_east'.format(f), ('time', 'z', 'y'),  getattr(self, '{}_east'.format(f)))
                add_variable('{}_south'.format(f), ('time', 'z', 'x'), getattr(self, '{}_south'.format(f)))
                add_variable('{}_north'.format(f), ('time', 'z', 'x'), getattr(self, '{}_north'.format(f)))
        
        nc.close()


if __name__ == '__main__':
    """
    Example / testing.
    """

    lbc = LBC_input(fields=['u','v','w'], itot=32, jtot=16, ktot=24, time=np.array([0, 3600]))

    lbc.u_west[:] = 2
    lbc.u_east[:] = 0
    lbc.u_south[:] = -1
    lbc.u_north[:] = 1

    # Other fields default to zero.

    lbc.to_netcdf('case_name')

