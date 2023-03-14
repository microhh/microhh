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

        if 'v' in fields:
            self.v_south = np.zeros((self.nt, ktot, itot), dtype)
            self.v_north = np.zeros((self.nt, ktot, itot), dtype)

        # Thermo/scalars/chemistry/...
        for f in fields:
            if f not in ['u', 'v']:
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

        nc_t = nc.createVariable('time',  self.dtype, ('time'))
        nc_t[:] = self.time

        if 'u' in self.fields:
            nc_uw = nc.createVariable('u_west',  self.dtype, ('time', 'z', 'y'))
            nc_ue = nc.createVariable('u_east',  self.dtype, ('time', 'z', 'y'))
            nc_uw[:] = self.u_west
            nc_ue[:] = self.u_east

        if 'v' in self.fields:
            nc_vs = nc.createVariable('v_south', self.dtype, ('time', 'z', 'x'))
            nc_vn = nc.createVariable('v_north', self.dtype, ('time', 'z', 'x'))
            nc_vs[:] = self.v_south
            nc_vn[:] = self.v_north

        for f in self.fields:
            if f not in ['u', 'v']:
                nc_sw = nc.createVariable('{}_west'.format(f),  self.dtype, ('time', 'z', 'y'))
                nc_se = nc.createVariable('{}_east'.format(f),  self.dtype, ('time', 'z', 'y'))
                nc_ss = nc.createVariable('{}_south'.format(f), self.dtype, ('time', 'z', 'x'))
                nc_sn = nc.createVariable('{}_north'.format(f), self.dtype, ('time', 'z', 'x'))

                nc_sw[:] = getattr(self, '{}_west'.format(f))
                nc_se[:] = getattr(self, '{}_east'.format(f))
                nc_ss[:] = getattr(self, '{}_south'.format(f))
                nc_sn[:] = getattr(self, '{}_north'.format(f))
        
        nc.close()

