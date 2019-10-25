import sys

sys.path.append('../python/')
import microhh_tools as mht

import taylorgreen.taylorgreen_test as taylorgreen
import drycblles.drycblles_test as drycblles
import bomex.bomex_test as bomex


precs = ['sp', 'dp']
modes = ['cpu']

for prec in precs:
    for mode in modes:
        microhh_exec = 'microhh_{}_{}'.format(prec, mode)
        experiment = '{}_{}'.format(prec, mode)

        taylorgreen.run_test(microhh_exec, prec, 'taylorgreen', experiment)
        taylorgreen.plot_test(microhh_exec, prec, 'taylorgreen', experiment)

        drycblles.run_test(microhh_exec, prec, 'drycblles', experiment)

        bomex.run_test(microhh_exec, prec, 'bomex', experiment)
