import sys
import shutil
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
sys.path.append('../../python/')
import microhh_tools as mht

no_opts = {}

opt_mpi = {
        'master': {'npx': 2, 'npy': 2}}

dict_rk = {
    'rk3': {'time': {'rkorder': 3}},
    'rk4': {'time': {'rkorder': 4}}}

dict_dt = {
    'dt1000': {'time': {'dt': 10.  }},
    'dt0500': {'time': {'dt':  5.  }},
    'dt0250': {'time': {'dt':  2.5 }},
    'dt0125': {'time': {'dt':  1.25}}}


if __name__ == '__main__':

    case_name = 'conservation'

    kwargs = dict([arg.split('=') for arg in sys.argv[2:]])

    if len(sys.argv) > 1:
        function_name = sys.argv[1]

        if function_name == 'run_case':
            mht.run_case(case_name, no_opts, opt_mpi, **kwargs)
        else:
            raise Exception('\"{}\" is an invalid option'.format(function_name))

    else:
        mht.run_case(case_name, no_opts, opt_mpi)

