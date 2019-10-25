import sys

sys.path.append('../python/')
import microhh_tools as mht

import taylorgreen.taylorgreen_test as taylorgreen
import drycblles.drycblles_test as drycblles

for prec in ['sp', 'dp']:
    for mode in ['cpu']:
        microhh_exec = 'microhh_{}_{}'.format(prec, mode)

        taylorgreen.run_test(microhh_exec, prec, 'taylorgreen')
        taylorgreen.plot_test(microhh_exec, prec, 'taylorgreen')

        drycblles.run_test(microhh_exec, prec, 'drycblles')

