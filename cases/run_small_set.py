import sys

sys.path.append('../python/')
import microhh_tools as mht

import taylorgreen.taylorgreen_test as tg
import taylorgreen.taylorgreen_test_plot as tgp

for prec in ['sp', 'dp']:
    for mode in ['cpu']:
        microhh_exec = 'microhh_{}_{}'.format(prec, mode)
        tg.run_test(microhh_exec, prec, 'taylorgreen')
        tgp.plot_test(microhh_exec, prec, 'taylorgreen')

