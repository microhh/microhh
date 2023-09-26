import matplotlib.pyplot as pl
import json
import numpy as np

pl.close('all')

cache_file = '../captures/diff_tke2__evisc_heat_float_128x128x128.cache'

with open (cache_file, 'r') as f:
    data = json.load(f)
timings = data['cache']

time = np.zeros(len(timings), dtype=float)
for i,it in enumerate(timings.values()):
    time[i] = it['time']

i_best = time.argmin()
iters = np.arange(time.size)

pl.figure()
pl.plot(iters, time)
pl.scatter(iters[i_best], time[i_best], label='Best')
pl.plot(iters, np.ones_like(iters)*time[i_best], 'k:')
pl.ylabel('Kernel execution time (s)')
pl.xlabel('Optimization iteration (-)')
pl.legend()
pl.ylim(0, pl.ylim()[1])
