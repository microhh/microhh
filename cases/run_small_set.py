import sys

sys.path.append('../python/')
import microhh_tools as mht

import taylorgreen.taylorgreen_test as taylorgreen
import moser180.moser180_test as moser180
import drycblles.drycblles_test as drycblles
import bomex.bomex_test as bomex


modes = ['cpu', 'cpumpi', 'gpu']
precs = ['dp', 'sp']

for prec in precs:
    for mode in modes:
        microhh_exec = 'microhh_{}_{}'.format(prec, mode)
        experiment = '{}_{}'.format(prec, mode)

        # taylorgreen.run(microhh_exec, prec, mode, 'taylorgreen', experiment)

        # moser180.run(microhh_exec, mode, 'moser180', experiment)
        moser180.run_restart(microhh_exec, mode, 'moser180', experiment)

        # drycblles.run(microhh_exec, mode, 'drycblles', experiment)
        drycblles.run_restart(microhh_exec, mode, 'drycblles', experiment)

        # bomex.run(microhh_exec, mode, 'bomex', experiment)
        bomex.run_restart(microhh_exec, mode, 'bomex', experiment)

