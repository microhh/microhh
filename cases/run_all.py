import sys
sys.path.append('../python/')
import microhh_tools as mht
import taylorgreen.taylorgreenconv as tg
import conservation.run_conservation as conv



options = {'itot' : 64, 'jtot' : 64, 'ktot' : 64}
mode_options = {'cpu' : {}, 'gpu' : {}, 'mpi' : {'npx' : 2, 'npy' : 4 }}
for prec in ['sp', 'dp']:
    for mode in ['cpu', 'mpi', 'gpu']:
        cases=[]
        rundir = mode + '_' + prec
        exec   = 'microhh_' + mode + '_' + prec

        # 1) Run all real cases at coarse resolution
        for case in ['arm', 'bomex', 'drycbl','drycblles', 'drycblles_subs', 'dryslope', 'dycoms', 'ekman', 'lasso', 'prandtlslope','rico', 'moser180']:
            cases.append(mht.Case(case,options={**options, **mode_options},rundir=rundir,keep=True))
        
        # 2) Test a RICO restart
        cases.append(mht.generator_restart(mht.Case('rico',options={**options, **mode_options, **{'endtime' : 60}})))

        # 3) Do conservation and taylorgreen test
        tg.main(exec, prec)
        conv.main()
        # 4) Run tests
        mht.test_cases(cases, exec)
