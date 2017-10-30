import netCDF4 as nc4
import numpy   as np
import struct  as st
import glob
import re
import sys

# -------------------------
# General help functions
# -------------------------

def _int_or_float_or_str(value):
    """ Helper function: convert a string to int/float/str """
    try:
        if ('.' in value):
            return float(value)
        else:
            return int(float(value))
    except:
        return value.rstrip()


def _convert_value(value):
    """ Helper function: convert namelist value or list """
    if ',' in value:
        value = value.split(',')
        return [_int_or_float_or_str(val) for val in value]
    else:
        return _int_or_float_or_str(value)


def _find_namelist_file():
    """ Helper function: automatically find the .ini file in the current directory """
    namelist_file = glob.glob('*.ini')
    if len(namelist_file) == 0:
        raise RuntimeError('Can\'t find any .ini files in the current directory!')
    if len(namelist_file) > 1:
        raise RuntimeError('There are multiple .ini files: {}'.format(namelist_file))
    else:
        return namelist_file[0]


def _process_endian(endian):
    if endian is None:
        endian = sys.byteorder
    if endian not in ['little', 'big']:
        raise ValueError('endian has to be \"little\" or \"big\"!')
    end = '<' if endian == 'little' else '>'
    return end


class _Namelist_group:
    def __init__(self):
        self.vars = {}

    def get(self, name, default=None):
        if name in self.vars.keys():
            return self.vars[name]
        elif default is not None:
            return default
        else:
            raise ValueError('No item \"{}\" in namelist'.format(name))

    def __getitem__(self, name):
        if name in self.vars.keys():
            return self.vars[name]
        else:
            raise ValueError('No item \"{}\" in namelist'.format(name))

    def __repr__(self):
        return self.vars.__repr__()


# -------------------------
# Classes and functions to read and write MicroHH things
# -------------------------

class Read_namelist:
    """ Reads a MicroHH .ini file to memory
        All available variables are accessible as e.g.:
            nl = Read_namelist()    # with no name specified, it searches for a .ini file in the current dir
            itot = nl['grid']['itot']
            enttime = nl['time']['endtime']
            printing e.g. nl['grid'] provides an overview of the available variables in a group

        Arguments:
            namelist_file -- path to .ini file (optional)
    """
    def __init__(self, namelist_file=None):
        if (namelist_file is None):
            namelist_file = _find_namelist_file()

        self.groups = {}   # Dictionary holding all the data
        with open(namelist_file) as f:
            for line in f:
                lstrip = line.strip()
                if (len(lstrip) > 0 and lstrip[0] != "#"):
                    if lstrip[0] == '[' and lstrip[-1] == ']':
                        curr_group_name = lstrip[1:-1]
                        self.groups[curr_group_name] = _Namelist_group()
                    elif ("=" in line):
                        var_name = lstrip.split('=')[0]
                        value = _convert_value(lstrip.split('=')[1])
                        self.groups[curr_group_name].vars[var_name] = value

    def __getitem__(self, name):
        if name in self.groups.keys():
            return self.groups[name]
        else:
            raise RuntimeError('Can\'t find group \"{}\" in .ini file'.format(name))

    def __repr__(self):
        return 'Available groups:\n{}'.format(', '.join(self.groups.keys()))


def replace_namelist_value(variable, new_value, group=None, namelist_file=None):
    """
    Replace a variables value in an existing namelist
    If no group is provided, all occurances are replaced (e.g.
    sampletime, which is used in stats, cross, dump, ...)

    Arguments:
        variable -- name of variable to replace
        new_value -- value to set
        group -- which group to change (optional)
        namelist_file -- path to .ini file (optional)
    """
    if namelist_file is None:
        namelist_file = _find_namelist_file()

    # Read the entire namelist to memory
    with open(namelist_file, "r") as source:
        lines = source.readlines()

    # Loop over lines, and replace if match
    curr_group = None
    with open(namelist_file, "w") as source:
        for line in lines:
            lstrip = line.strip()
            if (len(lstrip) > 0 and lstrip[0] != "#"):
                if lstrip[0] == '[' and lstrip[-1] == ']':
                    curr_group = lstrip[1:-1]

            if curr_group == group or group is None:
                source.write(re.sub(r'({}).*'.format(variable), r'\1={}'.format(new_value), line))
            else:
                source.write(line)


class Read_statistics:
    """ Read all the NetCDF statistics
        Example:
        f = Read_statistics('drycblles.default.0000000.nc')
        print(f) prints a list with the available variables
        The data can be accessed as either f['th'] or f.th, which returns the numpy array with data
        The variable names can be accessed as f.names['th'], the units as f.units['th'], the dimensions as f.dimensions['th']

        Arguments:
            stat_file -- path to statistics file
        """
    def __init__(self, stat_file):
        f = nc4.Dataset(stat_file, 'r')

        # Dictionaries which hold the variable names, units, etc.
        self.data       = {}
        self.units      = {}
        self.names      = {}
        self.dimensions = {}

        # For each variable in the NetCDF file, read all the content and info
        for var in f.variables:
            self.data[var]       = f.variables[var].__array__()
            self.units[var]      = f.variables[var].units
            self.names[var]      = f.variables[var].long_name
            self.dimensions[var] = f.variables[var].dimensions

        f.close()

    def __getitem__(self, name):
        if name in self.data.keys():
            return self.data[name]
        else:
            raise RuntimeError('Can\'t find variable \"{}\" in statistics file'.format(name))

    def __getattr__(self, name):
        if name in self.data.keys():
            return self.data[name]
        else:
            raise RuntimeError('Can\'t find variable \"{}\" in statistics file'.format(name))

    def __repr__(self):
        return 'Available variables:\n{}'.format(', '.join(self.names.keys()))


class Read_grid:
    """
    Read the grid file from MicroHH.
    If no file name is provided, grid.0000000 from the current directory is read
    If no dimensions (itot=None, ..) are provided, they are read from the namelist from the
    current directory.

    Arguments:
        itot -- number of grid points in x direction
        jtot -- number of grid points in y direction
        ktot -- number of grid points in z direction
        zsize -- vertical size of domain
        file_name -- path to grid file (optional)
        endian -- Endianess file ('little'/'big'), with None it is automatically determined

    Provides:
        Full and half level x, y, and z (i.e. x, xh, y, yh, z, zh) as class members

    Example usage:
        grid = Read_grid()
        grid = Read_grid(itot, jtot, ktot, zsize)
        x = grid.x
    """
    def __init__(self, itot=None, jtot=None, ktot=None, zsize=None, filename=None, endian=None):
        self.en  = _process_endian(endian)
        filename = 'grid.0000000' if filename is None else filename

        if None in [itot, jtot, ktot, zsize]:
            nl    = Read_namelist()
            itot  = nl['grid']['itot']
            jtot  = nl['grid']['jtot']
            ktot  = nl['grid']['ktot']
            zsize = nl['grid']['zsize']

        self.fin = open(filename, 'rb')
        self.x  = self.read(itot)
        self.xh = self.read(itot)
        self.y  = self.read(jtot)
        self.yh = self.read(jtot)
        self.z  = self.read(ktot)
        self.zh = self.read(ktot)
        self.fin.close()
        del self.fin

        self.itot = self.x.size
        self.jtot = self.y.size
        self.ktot = self.z.size

        self.dx = self.x[1]-self.x[0] if itot > 1 else self.xh[1]
        self.dy = self.y[1]-self.y[0] if jtot > 1 else self.yh[1]

        self.xsize = self.itot * self.dx
        self.ysize = self.jtot * self.dy
        self.zsize = zsize

    def read(self, n):
        return np.array(st.unpack('{0}{1}d'.format(self.en, n), self.fin.read(n*8)))


class Create_grid:
    """
    Create the grid as used in MicroHH.
    If no dimensions (itot=None, ..) are provided, they are read from the namelist from the
    current directory.

    Arguments:
        vgrid_file -- file which contains the vertical grid (the .prof file)
        order -- numerical order of the model (only 2 for now)
        itot -- number of grid points in x direction
        jtot -- number of grid points in y direction
        ktot -- number of grid points in z direction
        xsize -- x size of domain
        ysize -- y size of domain
        zsize -- vertical size of domain

    Provides:
        Full and half level x, y, and z (i.e. x, xh, y, yh, z, zh) as class members

    Example usage:
        grid = Create_grid('drycblles.prof')
        x = grid.x
    """
    def __init__(self, vgrid_file, order=2, itot=None, jtot=None, ktot=None, xsize=None, ysize=None, zsize=None):

        if (order != 2):
            raise RuntimeError('Only 2nd order (order=2) supported for now..')

        # If no grid parameters are provided, try to read from namelist
        if None in [itot, jtot, ktot, xsize, ysize, zsize]:
            nl    = Read_namelist()
            itot  = nl['grid']['itot']
            jtot  = nl['grid']['jtot']
            ktot  = nl['grid']['ktot']
            xsize = nl['grid']['xsize']
            ysize = nl['grid']['ysize']
            zsize = nl['grid']['zsize']

        self.dx = xsize / itot
        self.dy = ysize / jtot

        # Full level grid
        self.x = np.arange(0.5*self.dx, xsize, self.dx)
        self.y = np.arange(0.5*self.dy, ysize, self.dy)

        # Read vertical grid from prof file
        f = np.genfromtxt(vgrid_file, names=True)
        self.z = f['z']

        # Half level grid
        self.xh = np.arange(0, xsize, self.dx)
        self.yh = np.arange(0, ysize, self.dy)
        self.zh = np.zeros(ktot)
        self.zh[1:] = 0.5*(self.z[1:] + self.z[:-1])

        # Save other parameters to make it equal to the Read_grid() class
        self.itot = itot
        self.jtot = jtot
        self.ktot = ktot
        self.xsize = xsize
        self.ysize = ysize
        self.zsize = zsize


def read_restart_file(path, itot, jtot, ktot, endian=None):
    """
    Read a MicroHH restart file as a Nd Numpy array

    Arguments:
        itot -- number of grid points in x direction
        jtot -- number of grid points in y direction
        ktot -- number of grid points in z direction
        endian -- Endianess file ('little'/'big'), with None it is automatically determined

    Return:
        Numpy array (dimensions: [z,y,x]) with restart file
    """

    en = _process_endian(endian)

    f  = open(path, 'rb')
    if (ktot > 1):
        field = np.zeros((ktot, jtot, itot))
        for k in range(ktot):
            raw = f.read(itot*jtot*8)
            tmp = np.array(st.unpack('{0}{1}d'.format(en, itot*jtot), raw))
            field[k,:,:] = tmp.reshape((jtot, itot))[:,:]
        f.close()
    else:
        raw = f.read(itot*jtot*8)
        tmp = np.array(st.unpack('{0}{1}d'.format(en, itot*jtot), raw))
        field = tmp.reshape((jtot, itot))

    return field


def write_restart_file(path, array, itot, jtot, ktot, per_slice=True, endian=None):
    """
    Write a MicroHH restart file from a Nd Numpy array

    Arguments:
        path -- file name/path of restart file
        array -- Numpy array (dimensions: [z,y,x]) containing restart file data
        itot -- number of grid points in x direction
        jtot -- number of grid points in y direction
        ktot -- number of grid points in z direction
        per_slice -- write 3D field at once (with False) or per xy-slice (with True) (default=True)
        endian -- Endianess file ('little'/'big'), with None it is automatically determined
    """

    en = _process_endian(endian)

    if(per_slice):
        # Write level by level (less memory hungry.....)
        fout  = open(path, "wb")
        for k in range(ktot):
            tmp  = array[k,:,:].reshape(itot*jtot)
            tmp2 = st.pack('{0}{1}d'.format(en, tmp.size), *tmp)
            fout.write(tmp2)
        fout.close()
    else:
        # Write entire field at once (memory hungry....)
        tmp  = array.reshape(array.size)
        tmp2 = st.pack('{0}{1}d'.format(en, tmp.size), *tmp)
        fout = open(path, "wb")
        fout.write(tmp2)
        fout.close()


def get_cross_indices(variable, mode):
    """
    Find the cross-section indices

    Arguments:
        variable -- name of variable (..)
        mode -- cross section mode (\"xy\", \"xz\", \"yz\")

    Returns:
        List with cross section indices
    """

    if mode not in ['xy','xz','yz']:
        raise ValueError('\"mode\" should be in {\"xy\", \"xz\", \"yz\"}')

    # Get a list of all the cross-section files
    files = glob.glob('{}.{}.*.*'.format(variable, mode))
    if len(files) == 0:
        raise Exception('Cannot find any cross-section')

    # Get a list with all the cross-section files for one time
    time  = files[0].split('.')[-1]
    files = glob.glob('{}.{}.*.{}'.format(variable, mode, time))

    # Get the indices
    indices = [int(f.split('.')[-2]) for f in files]
    indices.sort()
    return indices


def write_output(file_name, data):
    """
    Write input files for MicroHH; e.g. initial vertical profiles or time series

    Arguments:
        file_name -- file name/path to save the file
        data -- dictionary with column names as keys, and list/numpy arrays as data

    Example:
        z = np.arange(100)
        u = np.random.random(100)

        data = {'z':z, 'u':u}
        write_output('case.prof', data)

        Creates file case.prof like:

                  z                     u
        +0.00000000000000E+00 +2.44790854944653E-01
        +1.00000000000000E+00 +6.45343637327276E-01
        +2.00000000000000E+00 +5.28456216105296E-01
        et cetera...
    """

    f = open(file_name, 'w')

    # Get size of data
    size = data[list(data.keys())[0]].size

    # Write header
    for var in data.keys():
        f.write('{0:^21s} '.format(var))
    f.write('\n')

    # Write data
    for k in range(size):
        for var in data.keys():
            f.write('{0:+1.14E} '.format(data[var][k]))
        f.write('\n')

    f.close()
