import sys

sys.path.append('../python/')
import microhh_tools as mht

import taylorgreen.taylorgreen_test as taylorgreen
import moser180.moser180_test as moser180
import drycblles.drycblles_test as drycblles
import bomex.bomex_test as bomex
import gabls1.gabls1_test as gabls1

#modes = ['cpu', 'cpumpi', 'gpu']
modes = ['cpu', 'cpumpi']
precs = ['dp', 'sp']

for prec in precs:
    for mode in modes:
        microhh_exec = 'microhh_{}_{}'.format(prec, mode)
        experiment = '{}_{}'.format(prec, mode)

        #taylorgreen.run(microhh_exec, prec, mode, 'taylorgreen', experiment)

        # Moser180 neutral channel flow
        mht.run_case('moser180',
                moser180.opt_small, moser180.opt_mpi,
                microhh_exec, mode, 'moser180', experiment)

        mht.run_restart('moser180',
                moser180.opt_small_restart, moser180.opt_mpi, moser180.dict_opts,
                microhh_exec, mode, 'moser180', experiment)


        # Dry convective boundary layer
        mht.run_case('drycblles',
                drycblles.opt_small, drycblles.opt_mpi,
                microhh_exec, mode, 'drycblles', experiment)

        mht.run_restart('drycblles',
                drycblles.opt_small, drycblles.opt_mpi, None,
                microhh_exec, mode, 'drycblles', experiment)


        # BOMEX LES intercomparison
        mht.run_case('bomex',
                bomex.opt_small, bomex.opt_mpi,
                microhh_exec, mode, 'bomex', experiment)

        mht.run_restart('bomex',
                bomex.opt_small_nostats, bomex.opt_mpi, bomex.dict_opts,
                microhh_exec, mode, 'bomex', experiment)


        # GABLS1 LES intercomparison
        mht.run_case('gabls1',
                gabls1.opt_small, gabls1.opt_mpi,
                microhh_exec, mode, 'gabls1', experiment)

        mht.run_restart('gabls1',
                gabls1.opt_small, gabls1.opt_mpi, gabls1.dict_opts,
                microhh_exec, mode, 'gabls1', experiment)
