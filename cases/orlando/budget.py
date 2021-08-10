
from netCDF4 import Dataset
from pylab import *
from numpy import *

x = Dataset("orlando.chemistry.0000000.nc")
cb = x.groups['default'].variables['chem_budget'][:]
