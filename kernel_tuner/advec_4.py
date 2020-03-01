import kernel_tuner
import numpy as np
from collections import OrderedDict

from helpers import Grid, Fields

kernel_string = """

    // 4th order interpolation
    const double ci0  = -1./16.;
    const double ci1  =  9./16.;
    const double ci2  =  9./16.;
    const double ci3  = -1./16.;

    // 4th order gradient
    const double cg0  =   1.;
    const double cg1  = -27.;
    const double cg2  =  27.;
    const double cg3  =  -1.;
    const double cgi  =   1./24.;

    __global__ 
    void advec_u_g(double* __restrict__ ut, double* __restrict__ u, 
                   double* __restrict__ v,  double* __restrict__ w,
                   double* __restrict__ dzi4, double dxi, double dyi, 
                   int jj, int kk,
                   int istart, int jstart, int kstart,
                   int iend,   int jend,   int kend)
        {
            const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
            const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
            const int k = blockIdx.z + kstart;

            const int ii1 = 1;
            const int ii2 = 2;
            const int ii3 = 3;
            const int jj1 = 1*jj;
            const int jj2 = 2*jj;
            const int jj3 = 3*jj;
            const int kk1 = 1*kk;
            const int kk2 = 2*kk;
            const int kk3 = 3*kk;

            if (i < iend && j < jend && k > kstart && k < kend-1)
            {
                const int ijk = i + j*jj + k*kk;
                ut[ijk] -= ( cg0*((ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]) * (ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]))
                           + cg1*((ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]) * (ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]))
                           + cg2*((ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]) * (ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]))
                           + cg3*((ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3]) * (ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3])) ) * cgi*dxi;

                ut[ijk] -= ( cg0*((ci0*v[ijk-ii2-jj1] + ci1*v[ijk-ii1-jj1] + ci2*v[ijk-jj1] + ci3*v[ijk+ii1-jj1]) * (ci0*u[ijk-jj3] + ci1*u[ijk-jj2] + ci2*u[ijk-jj1] + ci3*u[ijk    ]))
                           + cg1*((ci0*v[ijk-ii2    ] + ci1*v[ijk-ii1    ] + ci2*v[ijk    ] + ci3*v[ijk+ii1    ]) * (ci0*u[ijk-jj2] + ci1*u[ijk-jj1] + ci2*u[ijk    ] + ci3*u[ijk+jj1]))
                           + cg2*((ci0*v[ijk-ii2+jj1] + ci1*v[ijk-ii1+jj1] + ci2*v[ijk+jj1] + ci3*v[ijk+ii1+jj1]) * (ci0*u[ijk-jj1] + ci1*u[ijk    ] + ci2*u[ijk+jj1] + ci3*u[ijk+jj2]))
                           + cg3*((ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii1+jj2] + ci2*v[ijk+jj2] + ci3*v[ijk+ii1+jj2]) * (ci0*u[ijk    ] + ci1*u[ijk+jj1] + ci2*u[ijk+jj2] + ci3*u[ijk+jj3])) ) * cgi*dyi;

                ut[ijk] -= ( cg0*((ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk-kk1] + ci3*w[ijk+ii1-kk1]) * (ci0*u[ijk-kk3] + ci1*u[ijk-kk2] + ci2*u[ijk-kk1] + ci3*u[ijk    ]))
                           + cg1*((ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk    ] + ci3*w[ijk+ii1    ]) * (ci0*u[ijk-kk2] + ci1*u[ijk-kk1] + ci2*u[ijk    ] + ci3*u[ijk+kk1]))
                           + cg2*((ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk+kk1] + ci3*w[ijk+ii1+kk1]) * (ci0*u[ijk-kk1] + ci1*u[ijk    ] + ci2*u[ijk+kk1] + ci3*u[ijk+kk2]))
                           + cg3*((ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk+kk2] + ci3*w[ijk+ii1+kk2]) * (ci0*u[ijk    ] + ci1*u[ijk+kk1] + ci2*u[ijk+kk2] + ci3*u[ijk+kk3])) ) * dzi4[k];
            }
        }
"""

if __name__ == '__main__':

    TF = np.float64
    
    for gridsize in [32,64,128,256,512]:
        
        print('===========================')
        print('{0}x{0}x{0} grid points'.format(gridsize))
        print('===========================')
    
        grid   = Grid(3200, 3200, 3200, gridsize, gridsize, gridsize, 3, 3, TF)
        fields = Fields(['u','v','w'], grid.ncells, grid.ijcells, grid.kcells, TF)

        args = [
            fields.u.tend, fields.u.fld, fields.v.fld, fields.w.fld, 
            grid.dzi4, grid.dxi, grid.dyi, 
            grid.icells, grid.ijcells, 
            grid.istart, grid.jstart, grid.kstart, 
            grid.iend, grid.jend, grid.kend]
        
        tune_params = OrderedDict()
        tune_params["block_size_x"] = [32*i for i in range(1,9)]
        tune_params["block_size_y"] = [2**i for i in range(3)]

        tune = kernel_tuner.tune_kernel("advec_u_g", kernel_string, grid.ncells, args, tune_params)
