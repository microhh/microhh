import numpy as np
import matplotlib.pylab as pl

import microhh_tools as mht
import ib_tools as ibt

# 1. Read the namelist, and grid info
nl = mht.Read_namelist()
gr = mht.Read_grid(nl['grid']['itot'], nl['grid']['jtot'], nl['grid']['ktot'], nl['grid']['zsize'])

# 2. Read the polygon input
poly = read_polygon_structured('container_corners.txt')


