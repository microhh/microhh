import netCDF4 as nc4
import numpy   as np
import struct  as st
import matplotlib.pylab as pl

import glob
import fileinput
import re

# Constants
kappa = 0.4;        # von Karman constant
grav  = 9.81;       # Gravitational acceleration [m s-2]
Rd    = 287.04;     # Gas constant for dry air [J K-1 kg-1] 
Rv    = 461.5;      # Gas constant for water vapor [J K-1 kg-1]
cp    = 1005;       # Specific heat of air at constant pressure [J kg-1 K-1]
Lv    = 2.5e6;      # Latent heat of condensation or vaporization [J kg-1]
T0    = 273.15;     # Freezing / melting temperature [K]
p0    = 1.e5;       # Reference pressure [pa]
ep    = Rd/Rv;

class _Empty: pass

# ------------------------------------
# ----- General helper functions -----
# ------------------------------------
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

# -------------------------
# ----- MicroHH tools -----
# -------------------------
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
    def __init__(self, itot, jtot, ktot, zsize, filename=None, n_ghost=0, endian='little'):
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
        if (n_ghost > 0):
            for i in range(n_ghost):
                self.x  = np.insert(self.x,  0, self.x [0 ]-self.dx)
                self.x  = np.append(self.x,     self.x [-1]+self.dx)   
                self.xh = np.insert(self.xh, 0, self.xh[0 ]-self.dx) 
                self.xh = np.append(self.xh,    self.xh[-1]+self.dx) 
           
                self.y  = np.insert(self.y,  0, self.y [0 ]-self.dy)
                self.y  = np.append(self.y,     self.y [-1]+self.dy)   
                self.yh = np.insert(self.yh, 0, self.yh[0 ]-self.dy) 
                self.yh = np.append(self.yh,    self.yh[-1]+self.dy) 
            
            self.z  = np.insert(self.z,  0, -self.z[0])
            self.z  = np.append(self.z,  2*zsize-self.z[-1])
            self.zh = np.insert(self.zh, 0, 0) # undefined
            self.zh = np.append(self.zh, zsize)

        self.istart = n_ghost
        self.iend   = self.itot+n_ghost 
        self.jstart = n_ghost
        self.jend   = self.jtot+n_ghost
       
        ngz = 1 if n_ghost > 0 else 0 
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

# Get the cross-section indices
def get_cross_indices(variable, mode):
    files = glob.glob('{}.{}.*.*'.format(variable, mode))
    time  = files[0].split('.')[-1]
    files = glob.glob('{}.{}.*.{}'.format(variable, mode, time))
    return [int(f.split('.')[-2]) for f in files]

# Replace a variable value in the namelist
def replace_namelist_var(variable, new_value, namelist_file=None): 
    if (namelist_file is None):
        namelist_file = glob.glob('*.ini')[0]

    with open(namelist_file, "r") as source:
        lines = source.readlines()
    with open(namelist_file, "w") as source:
        for line in lines:
            source.write(re.sub(r'({}).*'.format(variable), r'\1={}'.format(new_value), line))

# Moist thermodynamics, identical to MicroHH
# Calculation saturation vapor pressure (pa)
def esat(T):
    c0 = 0.6105851e+03 
    c1 = 0.4440316e+02 
    c2 = 0.1430341e+01 
    c3 = 0.2641412e-01 
    c4 = 0.2995057e-03 
    c5 = 0.2031998e-05 
    c6 = 0.6936113e-08 
    c7 = 0.2564861e-11 
    c8 = -.3704404e-13 

    if (np.size(T) == 1):
        x = max(-80, T-T0)
    else:
        x = T-T0
        x[T<-80] = -80

    return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

# Calculation saturation specific humidity (kg/kg)
def qsat(p, T):
    return ep*esat(T)/(p-(1-ep)*esat(T))

# Exner function
def exner(p):
    return pow((p/p0),(Rd/cp));

# --------------------------
# ----- Plotting tools -----
# --------------------------
def remove_top_right_ax(ax=None):
    if ax is None:
        ax = pl.gca()

    ax.spines['right'].set_visible(False)
    ax.get_yaxis().tick_left()
    ax.spines['top'].set_visible(False)
    ax.get_xaxis().tick_bottom()

# Define the ColorBrewer (set 1) colors
c1    = '#333333'
c2    = '#e41a1c'
c3    = '#377eb8'
c4    = '#4daf4a'
c5    = '#984ea3'
c6    = '#ff7f00'
cc    = [c1,c2,c3,c4,c5,c6]

