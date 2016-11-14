import netCDF4 as nc4
import numpy   as np
import struct  as st
import matplotlib.pylab as pl

import glob
import fileinput
import re

class _Empty: pass

# Convert a string to int/float/str
def _int_or_float_or_str(value):
    try:
        if ('.' in value or 'e' in value):
            return float(value)
        else:
            return int(value)
    except:
        return value.rstrip()

# Convert namelist value or list
def _convert_value(value):
    if (',' in value):
        value = value.split(',') 
        return [_int_or_float_or_str(val) for val in value]
    else:
        return _int_or_float_or_str(value)

# Read namelist
class Read_namelist:
    def __init__(self, namelist_file=None):
        if (namelist_file is None):
            namelist_file = glob.glob('*.ini')[0]

        curr_group = None
        with open(namelist_file) as f:
            for line in f:
                if (len(line.strip()) > 0):
                    if (line.strip()[0] == '[' and line.strip()[-1] == ']'):
                        curr_group_name = line.strip()[1:-1]
                        curr_group = _Empty()
                        setattr(self, curr_group_name, curr_group)
                    elif (curr_group is not None):
                        setattr(curr_group, line.split('=')[0], _convert_value(line.split('=')[-1]))

# Read grid file MicroHH:
class Read_grid:
    def __init__(self, itot, jtot, ktot, zsize, filename=None, igc=0, jgc=0, kgc=0, endian='little'):
        if endian not in ['little', 'big']:
            raise ValueError('endian has to be \"little\" or \"big\"!')
        self.endian = '<' if endian == 'little' else '>'

        filename = 'grid.0000000' if filename is None else filename

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

        self.dx = self.x[1]-self.x[0]
        self.dy = self.y[1]-self.y[0] if jtot > 1 else 2*self.y[0]

        self.xsize = self.itot * self.dx
        self.ysize = self.jtot * self.dy
        self.zsize = zsize

        # Add ghost cells
        if (igc > 0):
            for i in range(igc):
                self.x  = np.insert(self.x,  0, self.x [0 ]-self.dx)
                self.x  = np.append(self.x,     self.x [-1]+self.dx)   
                self.xh = np.insert(self.xh, 0, self.xh[0 ]-self.dx) 
                self.xh = np.append(self.xh,    self.xh[-1]+self.dx) 
        if (jgc > 0):
            for i in range(jgc):
                self.y  = np.insert(self.y,  0, self.y [0 ]-self.dy)
                self.y  = np.append(self.y,     self.y [-1]+self.dy)   
                self.yh = np.insert(self.yh, 0, self.yh[0 ]-self.dy) 
                self.yh = np.append(self.yh,    self.yh[-1]+self.dy) 
        if (kgc > 0):
            self.z  = np.insert(self.z,  0, -self.z[0])
            self.z  = np.append(self.z,  2*zsize-self.z[-1])
            self.zh = np.insert(self.zh, 0, 0) # undefined
            self.zh = np.append(self.zh, zsize)

        self.istart = igc
        self.iend   = self.itot+igc
        self.jstart = jgc
        self.jend   = self.jtot+jgc
       
        ngz = 1 if kgc > 0 else 0
        self.kstart = ngz
        self.kend   = self.ktot+ngz

    def read(self, n):
        return np.array(st.unpack('{0}{1}d'.format(self.endian, n), self.fin.read(n*8)))

# Read all statistics
class Read_statistics:
    def __init__(self, stat_file):
        f = nc4.Dataset(stat_file, 'r')
        
        for var in f.variables:
            tmp      = _Empty()
            tmp.data = f.variables[var].__array__()
            tmp.unit = f.variables[var].getncattr('units')
            tmp.name = f.variables[var].getncattr('long_name')
            setattr(self, var, tmp)

        f.close()

# Read a single 3d restart file
def read_3d(file_in, itot, jtot, ktot, endian='little'):
    if endian not in ['little', 'big']:
        raise ValueError('endian has to be \"little\" or \"big\"!')
    en = '<' if endian == 'little' else '>'

    field = np.zeros((ktot, jtot, itot))

    f = open(file_in, 'rb')
    for k in range(ktot):
        raw = f.read(itot*jtot*8)
        tmp = np.array(st.unpack('{0}{1}d'.format(en, itot*jtot), raw))
        field[k,:,:] = tmp.reshape((jtot, itot))[:,:]
    f.close()

    return field

# Get the cross-section indices
def get_cross_indices(variable, mode):
    files = glob.glob('{}.{}.*.*'.format(variable, mode))
    time  = files[0].split('.')[-1]
    files = glob.glob('{}.{}.*.{}'.format(variable, mode, time))
    indices = [int(f.split('.')[-2]) for f in files]
    indices.sort()
    return indices

# Replace a variable value in the namelist
def replace_namelist_var(variable, new_value, namelist_file=None): 
    if (namelist_file is None):
        namelist_file = glob.glob('*.ini')[0]

    with open(namelist_file, "r") as source:
        lines = source.readlines()
    with open(namelist_file, "w") as source:
        for line in lines:
            source.write(re.sub(r'({}).*'.format(variable), r'\1={}'.format(new_value), line))
