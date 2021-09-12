from netCDF4 import Dataset
from pylab import *
import numpy as np
import sys
whitelist='oh_1_ho2_1 oh_1_o3_1 ho2_1_ho2_1 ho2_1_o3_1 no_1_o3_1 no2_1_oh_1 no2_1_o3_1 no_1_ho2_1 no2_1_no3_1 no_1_no3_1 ho2_1_no3_1 ch3o2_1_no_1 ch3o2_1_ho2_1 co_1_oh_1 ch4_1_oh_1 hcho_1_oh_1 isop_1_oh_1 isopao2_1_ch3o2_1 isopao2_1_ho2_1 isopao2_1_no_1 isopbo2_1_ch3o2_1 isopbo2_1_ho2_1 isopbo2_1_no_1 isopooh_1_oh_1 hald_1_oh_1 xo2_1_ho2_1 xo2_1_no_1 mvkmacr_1_oh_1 isop_1_o3_1 isop_1_no3_1'
wl = whitelist.split()
profs = []
with Dataset("orlando.default.0000000.nc") as x:
    zax = x.variables['z'][:]
    zaxh = x.variables['zh'][:]
    z = x.groups['covariances']
    for wls  in wl:
        profs.append(z.variables[wls][:])
