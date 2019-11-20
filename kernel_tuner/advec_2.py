import kernel_tuner
import numpy as np
from collections import OrderedDict

from helpers import Grid, Fields

kernel_string = """
    __global__
    void advec_s_g(double* __restrict__ st, double* __restrict__ s,
                   double* __restrict__ u,  double* __restrict__ v, double* __restrict__ w,
                   double* __restrict__ rhoref, double* __restrict__ rhorefh,
                   double* __restrict__ dzi, double dxi, double dyi,
                   int jj, int kk,
                   int istart, int jstart, int kstart,
                   int iend,   int jend,   int kend)
    {
        const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k  = blockIdx.z + kstart;
        const int ii = 1;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            st[ijk] +=
                - (  u[ijk+ii] * 0.5*(s[ijk   ] + s[ijk+ii])
                   - u[ijk   ] * 0.5*(s[ijk-ii] + s[ijk   ]) ) * dxi

                - (  v[ijk+jj] * 0.5*(s[ijk   ] + s[ijk+jj])
                   - v[ijk   ] * 0.5*(s[ijk-jj] + s[ijk   ]) ) * dyi

                - (  rhorefh[k+1] * w[ijk+kk] * 0.5*(s[ijk   ] + s[ijk+kk])
                   - rhorefh[k  ] * w[ijk   ] * 0.5*(s[ijk-kk] + s[ijk   ]) ) / rhoref[k] * dzi[k];
        }
    }
"""

if __name__ == '__main__':

    TF = np.float64
    
    for gridsize in [32,64,128,256,512]:
        
        print('===========================')
        print('{0}x{0}x{0} grid points'.format(gridsize))
        print('===========================')
    
        grid   = Grid(3200, 3200, 3200, gridsize, gridsize, gridsize, 2, 1, TF)
        fields = Fields(['u','v','w','s'], grid.ncells, grid.ijcells, grid.kcells, TF)

        args = [
            fields.s.tend, fields.s.fld, fields.u.fld, fields.v.fld, fields.w.fld, 
            fields.rhoref, fields.rhorefh, grid.dzi, 
            grid.dxi, grid.dyi, 
            grid.icells, grid.ijcells, 
            grid.istart, grid.jstart, grid.kstart, 
            grid.iend, grid.jend, grid.kend]
        
        tune_params = OrderedDict()
        tune_params["block_size_x"] = [32*i for i in range(1,9)]
        tune_params["block_size_y"] = [2**i for i in range(3)]

        tune = kernel_tuner.tune_kernel("advec_s_g", kernel_string, grid.ncells, args, tune_params)
