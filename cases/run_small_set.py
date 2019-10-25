import sys

sys.path.append('../python/')
import microhh_tools as mht

import taylorgreen.taylorgreen_test as taylorgreen
import drycblles.drycblles_test as drycblles
import bomex.bomex_test as bomex


#modes = ['cpu', 'cpumpi', 'gpu']

precs = ['dp', 'sp']
modes = ['cpu', 'cpumpi']

for prec in precs:
    for mode in modes:
        microhh_exec = 'microhh_{}_{}'.format(prec, mode)
        experiment = '{}_{}'.format(prec, mode)

        # taylorgreen.run_test(microhh_exec, prec, mode, 'taylorgreen', experiment)
        # taylorgreen.plot_test(microhh_exec, prec, mode, 'taylorgreen', experiment)

        # drycblles.run_test(microhh_exec, prec, mode, 'drycblles', experiment)
        drycblles.run_restart_test(microhh_exec, prec, mode, 'drycblles', experiment)

        # bomex.run_test(microhh_exec, prec, mode, 'bomex', experiment)
        bomex.run_restart_test(microhh_exec, prec, mode, 'bomex', experiment)

