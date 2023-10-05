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
import netCDF4 as nc4
import numpy as np
import os, sys

class LSM_input:
    def __init__(self, itot, jtot, ktot, TF=np.float64, debug=False, exclude_fields=[]):
        """
        Data structure for the required input for the MicroHH LSM.

        Arguments:
            itot (int) : Number of grid points in x-direction.
            jtot (int) : Number of grid points in y-direction.
            ktot (int) : Number of soil layers.
            TF (np.dtype) : Data type used by MicroHH (np.float64 for double precision,
                            np.float32 for single precision).
            debug (bool, default: False) : Switch to fill the emtpy fields with a large negative number,
                                           to ensure that every grid point is initialized before saving.
            exclude_fields (list) : Exclude default fields.
        """

        self.itot = itot
        self.jtot = jtot
        self.ktot = ktot
        self.debug = debug
        self.exclude_fields = exclude_fields

        # List of fields which are written to the binary input files for MicroHH
        self.fields_2d = [
                'c_veg', 'z0m', 'z0h', 'gD', 'lai',
                'rs_veg_min', 'rs_soil_min',
                'lambda_stable', 'lambda_unstable',
                'cs_veg', 'water_mask', 't_bot_water']
        self.fields_3d = [
                't_soil', 'theta_soil', 'index_soil', 'root_frac']

        # Horizonal grid (cell centers)
        self.x = np.zeros(itot, dtype=TF)
        self.y = np.zeros(jtot, dtype=TF)

        # Lat/lon coordinates of each grid point (not used by LES)
        self.lon = np.zeros((jtot, itot), dtype=TF)
        self.lat = np.zeros((jtot, itot), dtype=TF)

        # Create empty 2D/3D fields
        for fld in self.fields_2d:
            setattr(self, fld, np.zeros((jtot, itot), dtype=TF))

        for fld in self.fields_3d:
            setattr(self, fld, np.zeros((ktot, jtot, itot), dtype=TF))

        if debug:
            # Init all values at large negative number
            for field in self.fields_2d:
                data = getattr(self, field)
                data[:] = 1e12

            for field in self.fields_3d:
                data = getattr(self, field)
                data[:] = 1e12

    def check(self):
        """
        Check if all values have been set, in case of debug mode
        """

        if not self.debug:
            sys.exit('Can not check LSM input without debug mode...')
        else:
            for fld in self.fields_2d + self.fields_3d:
                if fld not in self.exclude_fields:
                    data = getattr(self, fld)
                    if np.any(data > 1e11):
                        print('WARNING: field \"{}\" is not initialised!'.format(fld))


    def save_binaries(self, path='.', allow_overwrite=False):
        """
        Write all required MicroHH input fields in binary format

        Arguments:
            path (str) : File path of the output.
            allow_overwrite (bool, default: False) : allow overwriting of existing files.
        """

        def save_bin(data, bin_file):
            if not allow_overwrite and os.path.exists(bin_file):
                sys.exit('Binary file \"{}\" already exists!'.format(bin_file))
            else:
                data.tofile(bin_file)

        for fld in self.fields_2d + self.fields_3d:
            if fld not in self.exclude_fields:
                data = getattr(self, fld)
                save_bin(data, '{}.0000000'.format(os.path.join(path,fld)))


    def save_netcdf(self, nc_file, allow_overwrite=False):
        """
        Save MicroHH input to NetCDF file, for visualisation et al.

        Arguments:
            nc_file (str) : name of output NetCDF file.
            allow_overwrite (bool, default: False) : allow overwriting of existing files.
        """

        if not allow_overwrite and os.path.exists(nc_file):
            sys.exit('NetCDF file \"{}\" already exists!'.format(nc_file))

        nc = nc4.Dataset(nc_file, 'w')

        dimx = nc.createDimension('x', self.itot)
        dimy = nc.createDimension('y', self.jtot)
        dimz = nc.createDimension('z', self.ktot)

        var_x = nc.createVariable('x', 'f8', 'x')
        var_y = nc.createVariable('y', 'f8', 'y')

        var_lon = nc.createVariable('lon', 'f8', ('y','x'))
        var_lat = nc.createVariable('lat', 'f8', ('y','x'))

        var_x[:] = self.x[:]
        var_y[:] = self.y[:]

        var_lon[:,:] = self.lon[:,:]
        var_lat[:,:] = self.lat[:,:]

        # Fields needed for offline LSM:
        for field in self.fields_2d:
            if field not in self.exclude_fields:
                data = getattr(self, field)
                var  = nc.createVariable(field, 'f8', ('y','x'))
                var[:] = data[:]

        for field in self.fields_3d:
            if field not in self.exclude_fields:
                data = getattr(self, field)
                var  = nc.createVariable(field, 'f8', ('z','y','x'))
                var[:] = data[:]

        nc.close()



if __name__ == '__main__':
    #
    # Example MicroHH input
    #
    itot = 64
    jtot = 48
    ktot_soil = 4

    lsm_input = LSM_input(itot, jtot, ktot_soil)

    # Set values
    # Required fields are:
    # 2D: c_veg, z0m, z0h, gD, lai, rs_veg_min, rs_soil_min, lambda_skin, index_veg, water_mask
    # 3D: t_soil, theta_soil, index_soil, root_frac
    dxy = 25
    lsm_input.x = np.arange(0.5*dxy, dxy*itot, dxy)
    lsm_input.y = np.arange(0.5*dxy, dxy*jtot, dxy)

    lsm_input.lai[:jtot//2,:] = 2.5
    lsm_input.lai[jtot//2:,:] = 4

    lsm_input.z0m[:,:] = 0.1

    lsm_input.theta_soil[0,:,:] = 0.2
    lsm_input.theta_soil[1,:,:] = 0.25
    lsm_input.theta_soil[2,:,:] = 0.3
    lsm_input.theta_soil[3,:,:] = 0.35
    # et cetera..

    # Save in binary and NetCDF format
    lsm_input.save_binaries(path='', allow_overwrite=True)
    lsm_input.save_netcdf('lsm_input.nc', allow_overwrite=True)
