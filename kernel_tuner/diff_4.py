import kernel_tuner
import numpy as np
from collections import OrderedDict
from numba import jit

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


kernel_string = """
    __global__
    void diff_c_g(double* __restrict__ const at, const double* __restrict__ const a,
                  const double* __restrict__ const dzi4, const double* __restrict__ const dzhi4,
                  const double dxi, const double dyi, const double visc,
                  const int jj,     const int kk,
                  const int istart, const int jstart, const int kstart,
                  const int iend,   const int jend,   const int kend,
                  const int ngc)
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

kernel_string_smem = """
    __global__
    void diff_c_g(double* __restrict__ const at, const double* __restrict__ const a,
                  const double* __restrict__ const dzi4, const double* __restrict__ const dzhi4,
                  const double dxi, const double dyi, const double visc,
                  const int jj,     const int kk,
                  const int istart, const int jstart, const int kstart,
                  const int iend,   const int jend,   const int kend,
                  const int ngc)
    {
        const int tx = threadIdx.x;
        const int ty = threadIdx.y;
        const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k  = blockIdx.z + kstart;
        const int blockxpad = blockDim.x+2*ngc;

        // DANGER DANGER ghost cells harcoded..
        __shared__ double as[(block_size_y+6)*(block_size_x+6)];

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
            const int ijk  = i + j*jj + k*kk;              // index in global memory
            const int ijks = (tx+ngc)+(ty+ngc)*blockxpad;  // Same location in 2d shared mem

            const int ii1 = 1;
            const int ii2 = 2;
            const int ii3 = 3;
            const int jj3 = 3*jj;
            const int kk1 = 1*kk;
            const int kk2 = 2*kk;
            const int kk3 = 3*kk;

            const int jjs1 = 1*blockxpad;
            const int jjs2 = 2*blockxpad;
            const int jjs3 = 3*blockxpad;

            const double dxidxi = dxi*dxi;
            const double dyidyi = dyi*dyi;

            // Read horizontal slice to shared memory
            as[ijks] = a[ijk];

            if(ty < ngc)
                as[ijks-jjs3] = a[ijk-jj3];
            if(ty >= blockDim.y-ngc)
                as[ijks+jjs3] = a[ijk+jj3];

            if(tx < ngc)
                as[ijks-ii3] = a[ijk-ii3];
            if(tx >= blockDim.x-ngc)
                as[ijks+ii3] = a[ijk+ii3];

            __syncthreads();

            // bottom boundary
            if (k == kstart)
            {
                at[ijk] += visc * (cdg3*as[ijks-ii3 ] + cdg2*as[ijks-ii2 ] + cdg1*as[ijks-ii1 ] + cdg0*as[ijks] + cdg1*as[ijks+ii1 ] + cdg2*as[ijks+ii2 ] + cdg3*as[ijks+ii3 ])*dxidxi;
                at[ijk] += visc * (cdg3*as[ijks-jjs3] + cdg2*as[ijks-jjs2] + cdg1*as[ijks-jjs1] + cdg0*as[ijks] + cdg1*as[ijks+jjs1] + cdg2*as[ijks+jjs2] + cdg3*as[ijks+jjs3])*dyidyi;
                at[ijk] += visc * ( cg0*(bg0*a[ijk-kk2] + bg1*a[ijk-kk1] + bg2*as[ijks  ] + bg3*a[ijk+kk1]) * dzhi4[k-1]
                                  + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*as[ijks  ] + cg3*a[ijk+kk1]) * dzhi4[k  ]
                                  + cg2*(cg0*a[ijk-kk1] + cg1*as[ijks  ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[k+1]
                                  + cg3*(cg0*as[ijks  ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzhi4[k+2] ) * dzi4[k];
            }
            // top boundary
            else if (k == kend-1)
            {
                at[ijk] += visc * (cdg3*as[ijks-ii3 ] + cdg2*as[ijks-ii2 ] + cdg1*as[ijks-ii1 ] + cdg0*as[ijks] + cdg1*as[ijks+ii1 ] + cdg2*as[ijks+ii2 ] + cdg3*as[ijks+ii3 ])*dxidxi;
                at[ijk] += visc * (cdg3*as[ijks-jjs3] + cdg2*as[ijks-jjs2] + cdg1*as[ijks-jjs1] + cdg0*as[ijks] + cdg1*as[ijks+jjs1] + cdg2*as[ijks+jjs2] + cdg3*as[ijks+jjs3])*dyidyi;
                at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*as[ijks  ]) * dzhi4[k-1]
                                  + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*as[ijks  ] + cg3*a[ijk+kk1]) * dzhi4[k  ]
                                  + cg2*(cg0*a[ijk-kk1] + cg1*as[ijks  ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[k+1]
                                  + cg3*(tg0*a[ijk-kk1] + tg1*as[ijks  ] + tg2*a[ijk+kk1] + tg3*a[ijk+kk2]) * dzhi4[k+2] ) * dzi4[k];
            }
            // interior
            else
            {
                at[ijk] += visc * (cdg3*as[ijks-ii3 ] + cdg2*as[ijks-ii2 ] + cdg1*as[ijks-ii1 ] + cdg0*as[ijks] + cdg1*as[ijks+ii1 ] + cdg2*as[ijks+ii2 ] + cdg3*as[ijks+ii3 ])*dxidxi;
                at[ijk] += visc * (cdg3*as[ijks-jjs3] + cdg2*as[ijks-jjs2] + cdg1*as[ijks-jjs1] + cdg0*as[ijks] + cdg1*as[ijks+jjs1] + cdg2*as[ijks+jjs2] + cdg3*as[ijks+jjs3])*dyidyi;
                at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*as[ijks  ]) * dzhi4[k-1]
                                  + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*as[ijks  ] + cg3*a[ijk+kk1]) * dzhi4[k  ]
                                  + cg2*(cg0*a[ijk-kk1] + cg1*as[ijks  ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[k+1]
                                  + cg3*(cg0*as[ijks  ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzhi4[k+2] ) * dzi4[k];
            }
        }
    }

"""

if __name__ == '__main__':

    double = np.float64
    kernel = kernel_string
    
    for gridsize in [32,64,128,256,512]:
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
                grid.iend, grid.jend, grid.kend, ngc]

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
                params = { "block_size_x": 4, "block_size_y": 4 }
                results = kernel_tuner.run_kernel(
                    "diff_c_g", kernel, (grid.itot, grid.jtot, grid.ktot), args, params)

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
                tune_params["block_size_y"] = [2**i for i in range(2,5)]

                # True answer from the Python code
                answer = len(args)*[None]
                answer[0] = fields_ref.s.tend

                tune = kernel_tuner.tune_kernel(
                    "diff_c_g", kernel, (grid.itot, grid.jtot, grid.ktot), args, tune_params, answer=answer, atol=1e-14)
