import sys, os

sys.path.append('../python/')
import microhh_tools as mht

#import taylorgreen.taylorgreenconv as tg
#import conservation.run_conservation as conv

modes = ['cpu', 'cpumpi', 'gpu']
precs = ['dp', 'sp']

les_cases   = ['bomex', 'drycblles', 'eady', 'gabls1', 'rico', 'sullivan2011', 'rcemip']  # dycoms+lasso missing
dns_cases   = ['drycbl', 'ekman', 'drycblslope', 'moser180', 'moser600']    # prandtlslope missing

# Cases which require an additional preprocessing script.
additional_pre = {'rcemip': {'link_coefficients.py': None}}

les_options = {
        'grid': {'itot': 16, 'jtot': 16, 'xsize': 1600, 'ysize': 1600},
        'time': {'endtime': 300}}

dns_options = {
        'grid': {'itot': 8, 'jtot': 8},
        'time': {'endtime': 1}}

mpi_options = {
        'master': {'npx': 2, 'npy': 2}}

failed = 0

for prec in precs:
    for mode in modes:
        # Likely MicroHH binary locations:
        ex1 = 'microhh_{}_{}'.format(prec, mode)
        ex2 = '../build_{}_{}/microhh'.format(prec, mode)

        if os.path.exists(ex1):
            microhh_exec = ex1
        elif os.path.exists(ex2):
            microhh_exec = ex2
        else:
            raise Exception('Can not find an executable for \"{}\" + \"{}\"'.format(prec, mode))

        experiment   = '{}_{}'.format(prec, mode)

        for case in les_cases:
            pre = {}
            if case in additional_pre:
                pre = additional_pre[case]

            failed += mht.run_case(case,
                    les_options, mpi_options,
                    microhh_exec, mode, case, experiment,
                    additional_pre_py=pre)

        for case in dns_cases:
            failed += mht.run_case(case,
                    dns_options, mpi_options,
                    microhh_exec, mode, case, experiment)

#        # 3) Do conservation and taylorgreen test
#        tg.main(exec, prec)
#        conv.main()

if failed == 0:
    print('Hurray! Zero cases failed...')
else:
    print('Uh oh! {} cases failed!'.format(failed))
