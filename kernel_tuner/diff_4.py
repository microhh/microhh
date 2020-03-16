import kernel_tuner
import numpy as np
from collections import OrderedDict
from numba import jit
import json

from helpers import Grid, Fields

@jit(nopython=True, nogil=True, fastmath=True)
def diff_c(
        at, a, dzi4, dzhi4, dxi, dyi, visc, jj, kk,
        istart, jstart, kstart,
        iend, jend, kend):

    cg0  =   1./24.
    cg1  = -27./24.
    cg2  =  27./24.
    cg3  =  -1./24.

    bg0  = -23./24.
    bg1  =  21./24.
    bg2  =   3./24.
    bg3  =  -1./24.

    tg0  =   1./24.
    tg1  =  -3./24.
    tg2  = -21./24.
    tg3  =  23./24.

    cdg0 = -1460./576.
    cdg1 =   783./576.
    cdg2 =   -54./576.
    cdg3 =     1./576.

    ii1 = 1
    ii2 = 2
    ii3 = 3
    jj1 = 1*jj
    jj2 = 2*jj
    jj3 = 3*jj
    kk1 = 1*kk
    kk2 = 2*kk
    kk3 = 3*kk

    dxidxi = dxi*dxi
    dyidyi = dyi*dyi

    for k in range(kstart, kend):
        for j in range(jstart, jend):
            for i in range(istart, iend):
                ijk = i + j*jj + k*kk

                # bottom boundary
                if k == kstart:
                    at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
                    at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
                    at[ijk] += visc * ( cg0*(bg0*a[ijk-kk2] + bg1*a[ijk-kk1] + bg2*a[ijk    ] + bg3*a[ijk+kk1]) * dzhi4[k-1]
                                      + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzhi4[k  ]
                                      + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[k+1]
                                      + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzhi4[k+2] ) * dzi4[k];

                # top boundary
                elif k == kend-1:
                    at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
                    at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
                    at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzhi4[k-1]
                                      + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzhi4[k  ]
                                      + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[k+1]
                                      + cg3*(tg0*a[ijk-kk1] + tg1*a[ijk    ] + tg2*a[ijk+kk1] + tg3*a[ijk+kk2]) * dzhi4[k+2] ) * dzi4[k];

                # interior
                else:
                    at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
                    at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
                    at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzhi4[k-1]
                                      + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzhi4[k  ]
                                      + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[k+1]
                                      + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzhi4[k+2] ) * dzi4[k];


if __name__ == '__main__':

    double = np.float64
    kernel_name = 'diff_c_g_smem'     # `diff_c_g` or `diff_c_g_smem`

    # Load CUDA source code
    with open('diff_4.cu', 'r') as f:
        kernel_string = f.read()

    # Test different grid dimensions (all as N^3)
    for gridsize in np.arange(32, 513, 32):

        print('===========================')
        print('{0}x{0}x{0} grid points'.format(gridsize))
        print('===========================')

        for manual_check in [True, False]:
    
            grid       = Grid(3200, 3200, 3200, gridsize, gridsize, gridsize, 3, 3, double)
            fields     = Fields(['s'], grid.ncells, grid.ijcells, grid.kcells, double)
            fields_ref = Fields(['s'], grid.ncells, grid.ijcells, grid.kcells, double)

            visc = double(1.e-5)
            ngc  = np.int32(3)

            # Argument list of CUDA kernel
            args = [
                fields.s.tend, fields.s.fld,
                grid.dzi4, grid.dzhi4,
                grid.dxi, grid.dyi, visc,
                grid.icells, grid.ijcells, 
                grid.istart, grid.jstart, grid.kstart, 
                grid.iend, grid.jend, grid.kend,
                grid.icells, grid.jcells, ngc]

            # Validate results with Python version of code
            diff_c(
                fields_ref.s.tend, fields_ref.s.fld,
                grid.dzi4, grid.dzhi4,
                grid.dxi, grid.dyi, visc,
                grid.icells, grid.ijcells, 
                grid.istart, grid.jstart, grid.kstart, 
                grid.iend, grid.jend, grid.kend)
            
            if manual_check:

                # CUDA
                params = { 'block_size_x': 32, 'block_size_y': 4 }
                results = kernel_tuner.run_kernel(
                    kernel_name, kernel_string, (grid.itot, grid.jtot, grid.ktot), args, params)

                # Compare results
                if np.allclose(results[0], fields_ref.s.tend, atol=1e-14):
                    print('SUCCESS, CUDA and Python results are identical!')
                    print('Max diff = {}'.format(np.absolute(results[0]-fields_ref.s.tend).max()))
                else:
                    print('ERROR: CUDA and Python results not identical')
                    print('Max diff = {}'.format(np.absolute(results[0]-fields_ref.s.tend).max()))

            else:

                # Tune parameters
                tune_params = OrderedDict()
                tune_params["block_size_x"] = [32*i for i in range(1,9)]
                tune_params["block_size_y"] = [2**i for i in range(0,5)]

                # True answer from the Python code
                answer = len(args)*[None]
                answer[0] = fields_ref.s.tend

                tune = kernel_tuner.tune_kernel(
                    kernel_name, kernel_string, (grid.itot, grid.jtot, grid.ktot), args, tune_params, answer=answer, atol=1e-14)

                # Write results to log file
                with open('{0}_{1:03d}.json'.format(kernel_name, gridsize), 'w') as fp:
                    json.dump(tune, fp)
