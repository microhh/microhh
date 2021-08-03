import sys

sys.path.append('../python/')
import microhh_tools as mht

#import taylorgreen.taylorgreenconv as tg
#import conservation.run_conservation as conv

modes = ['cpu', 'cpumpi', 'gpu']
#modes = ['cpu', 'cpumpi']
precs = ['dp', 'sp']

les_cases   = ['arm', 'bomex', 'drycblles', 'eady', 'gabls1', 'rico', 'sullivan2011']  # dycoms+lasso+rcemip+lasso missing
dns_cases   = ['drycbl', 'ekman', 'drycblslope', 'moser180', 'moser600']    # prandtlslope missing

les_options = {
        'grid': {'itot': 16, 'jtot': 16, 'xsize': 1600, 'ysize': 1600},
        'time': {'endtime': 300}}

dns_options = {
        'grid': {'itot': 8, 'jtot': 8},
        'time': {'endtime': 1}}

mpi_options = {
        'master': {'npx': 2, 'npy': 2}}

for prec in precs:
    for mode in modes:
        microhh_exec = 'microhh_{}_{}'.format(prec, mode)
        experiment   = '{}_{}'.format(prec, mode)

        for case in les_cases:
            mht.run_case(case,
                    les_options, mpi_options,
                    microhh_exec, mode, case, experiment)

        for case in dns_cases:
            mht.run_case(case,
                    dns_options, mpi_options,
                    microhh_exec, mode, case, experiment)

#        # 3) Do conservation and taylorgreen test
#        tg.main(exec, prec)
#        conv.main()
