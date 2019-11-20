import kernel_tuner
import numpy as np
from collections import OrderedDict

from helpers import Grid, Fields

kernel_string = """
    __global__
    void diff_c_g(double* __restrict__ const at, const double* __restrict__ const a,
                  const double* __restrict__ const dzi4, const double* __restrict__ const dzhi4,
                  const double dxi, const double dyi, const double visc,
                  const int jj,     const int kk,
                  const int istart, const int jstart, const int kstart,
                  const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        constexpr double cg0  =   1./24.;
        constexpr double cg1  = -27./24.;
        constexpr double cg2  =  27./24.;
        constexpr double cg3  =  -1./24.;

        constexpr double bg0  = -23./24.;
        constexpr double bg1  =  21./24.;
        constexpr double bg2  =   3./24.;
        constexpr double bg3  =  -1./24.;

        constexpr double tg0  =   1./24.;
        constexpr double tg1  =  -3./24.;
        constexpr double tg2  = -21./24.;
        constexpr double tg3  =  23./24.;

        constexpr double cdg0 = -1460./576.;
        constexpr double cdg1 =   783./576.;
        constexpr double cdg2 =   -54./576.;
        constexpr double cdg3 =     1./576.;
 
        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;

            const int ii1 = 1;
            const int ii2 = 2;
            const int ii3 = 3;
            const int jj1 = 1*jj;
            const int jj2 = 2*jj;
            const int jj3 = 3*jj;
            const int kk1 = 1*kk;
            const int kk2 = 2*kk;
            const int kk3 = 3*kk;

            const double dxidxi = dxi*dxi;
            const double dyidyi = dyi*dyi;

            // bottom boundary
            if (k == kstart)
            {
                at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
                at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
                at[ijk] += visc * ( cg0*(bg0*a[ijk-kk2] + bg1*a[ijk-kk1] + bg2*a[ijk    ] + bg3*a[ijk+kk1]) * dzhi4[k-1]
                            + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzhi4[k  ]
                            + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[k+1]
                            + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzhi4[k+2] ) * dzi4[k];
            }
            // top boundary
            else if (k == kend-1)
            {
                at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
                at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
                at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzhi4[k-1]
                             + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzhi4[k  ]
                             + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[k+1]
                             + cg3*(tg0*a[ijk-kk1] + tg1*a[ijk    ] + tg2*a[ijk+kk1] + tg3*a[ijk+kk2]) * dzhi4[k+2] ) * dzi4[k];
            }
            // interior
            else
            {
                at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
                at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
                at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzhi4[k-1]
                             + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzhi4[k  ]
                             + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[k+1]
                             + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzhi4[k+2] ) * dzi4[k];
            }
        }
    }


"""

if __name__ == '__main__':

    double = np.float64
    
    for gridsize in [32,64,128,256,512]:
        
        print('===========================')
        print('{0}x{0}x{0} grid points'.format(gridsize))
        print('===========================')
    
        grid   = Grid(3200, 3200, 3200, gridsize, gridsize, gridsize, 3, 3, double)
        fields = Fields(['s'], grid.ncells, grid.ijcells, grid.kcells, double)

        visc = double(1.e-5)

        args = [
            fields.s.tend, fields.s.fld,
            grid.dzi4, grid.dzhi4,
            grid.dxi, grid.dyi, visc,
            grid.icells, grid.ijcells, 
            grid.istart, grid.jstart, grid.kstart, 
            grid.iend, grid.jend, grid.kend]
        
        tune_params = OrderedDict()
        tune_params["block_size_x"] = [32*i for i in range(1,9)]
        tune_params["block_size_y"] = [2**i for i in range(3)]

        tune = kernel_tuner.tune_kernel("diff_c_g", kernel_string, grid.ncells, args, tune_params)
