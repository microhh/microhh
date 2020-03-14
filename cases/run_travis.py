import sys
import os

sys.path.append('../python/')
import microhh_tools as mht

modes = ['cpu']
precs = ['dp'] # ['dp', 'sp']

#
# Link executables to working directory
#
for prec in precs:
    for mode in modes:
        src = '../build_{}_{}/microhh'.format(prec, mode)
        dst = 'microhh_{}_{}'.format(prec, mode)
        if not os.path.exists(dst):
            os.symlink(src, dst)

#
# Short tests of init + run phase, to test if cases start
#
les_cases   = ['arm', 'bomex', 'drycblles', 'eady', 'gabls1', 'rico', 'sullivan2011']  # dycoms+lasso+rcemip missing
dns_cases   = ['drycbl', 'ekman', 'drycblslope', 'moser180', 'moser600']    # prandtlslope missing

les_options = {
        'grid': {'itot': 8, 'jtot': 8, 'xsize': 800, 'ysize': 800},
        'time': {'endtime': 200, 'savetime': 100}}

dns_options = {
        'grid': {'itot': 8, 'jtot': 8},
        'time': {'endtime': 2, 'savetime': 1}}

mpi_options = {
        'master': {'npx': 2, 'npy': 2}}

print('-----------------')
print('Running all cases')
print('-----------------')

err = 0

for prec in precs:
    for mode in modes:
        microhh_exec = 'microhh_{}_{}'.format(prec, mode)
        experiment   = '{}_{}'.format(prec, mode)

        for case in les_cases:
            err += mht.run_case(case,
                    les_options, mpi_options,
                    microhh_exec, mode, case, experiment)

        for case in dns_cases:
            err += mht.run_case(case,
                    dns_options, mpi_options,
                    microhh_exec, mode, case, experiment)

#
# Restart tests for a DNS + LES case
#
no_stats = {
        'stats': {'swstats': 0}, 'cross': {'swcross': 0}, 'column': {'swcolumn': 0}}

dns_perturbations = {
        '4th':    {},
        '4th-4m': {'grid': {'swspatialorder': 4}, 'advec': {'swadvec': '4m'}, 'diff': {'swdiff': 4}},
        '2nd':    {'grid': {'swspatialorder': 2}, 'advec': {'swadvec': 2   }, 'diff': {'swdiff': 2}}}

dns_options.update(no_stats)

les_perturbations = {
        'default': {},
        'vapor': {'thermo': {'swthermo': 'vapor'}},
        'fixed_basestate': {'thermo': {'swupdatebasestate': 0}},
        'advec_2i3': {'advec': {'swadvec': '2i3'}},
        'advec_2i4': {'advec': {'swadvec': '2i4'}}}

les_options.update(no_stats)

print('---------------------')
print('Running restart tests')
print('---------------------')

for prec in precs:
    for mode in modes:
        microhh_exec = 'microhh_{}_{}'.format(prec, mode)
        experiment   = '{}_{}'.format(prec, mode)

        err += mht.run_restart('drycbl',
                dns_options, mpi_options, dns_perturbations,
                microhh_exec, mode, 'drycbl', experiment)

        err += mht.run_restart('bomex',
                les_options, mpi_options, les_perturbations,
                microhh_exec, mode, 'bomex', experiment)

if err > 0:
    sys.exit('One of more travis case tests failed!')
