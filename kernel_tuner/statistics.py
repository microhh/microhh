import matplotlib.pyplot as pl
import numpy as np
import json
import glob

pl.close('all')
pl.ion()

def get_timings(kernel_name, gridsizes):

    dt = np.zeros_like(gridsizes, dtype=float)
    
    for i,gridsize in enumerate(gridsizes):
        with open( '{0}_{1:03d}.json'.format(kernel_name, gridsize) ) as f:
            data = json.load(f)
            timings = data[0]
    
            fastest = 1e9
            for timing in timings:
                fastest = min(fastest, timing['time'])
            dt[i] = fastest

    return dt


if __name__ == '__main__':
    gridsize = np.arange(32, 513, 32)
    normal = get_timings('diff_c_g',      gridsize)
    smem   = get_timings('diff_c_g_smem', gridsize)

    fac = gridsize**3

    pl.figure(figsize=(8,4))
    pl.subplot(121)
    pl.plot(gridsize, normal/fac, 'k-x', label='non smem')
    pl.plot(gridsize, smem  /fac, 'r-x', label='smem')
    pl.ylim(0, 2e-7)
    pl.ylabel('time/gridpoint (s)')
    pl.xlabel('gridpoints (-)')
    pl.legend()
    pl.grid()

    pl.subplot(122)
    pl.plot(gridsize, normal/smem, 'k-x')
    pl.ylabel('non_smem/smem (-)')
    pl.xlabel('gridpoints (-)')
    pl.grid()

    pl.tight_layout()
