import numpy as np
import pylab as pl
with open('input_chem_jvals') as f:
    line = f.readline()
    line = f.readline()
    line = f.readline()
    line = f.readline()
    header = line.split(' ')
    lines = f.readlines()
vals = []
for line in lines:
    vals.append([float(i) for i in line.split()])
vals = np.array(vals)
for i in range(8):
    f,ax = pl.subplots()
    ax.plot(vals[:,0],vals[:,i+2])




