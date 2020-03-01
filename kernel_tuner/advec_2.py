import kernel_tuner
import numpy as np
from collections import OrderedDict
from numba import jit

from helpers import Grid, Fields


@jit(nopython=True, nogil=True, fastmath=True)
def advec_s(
        st, s, u, v, w,
        rhoref, rhorefh, dzi, dxi, dyi, jj, kk,
        istart, jstart, kstart,
        iend,   jend,   kend):

    ii = 1
    
    for k in range(kstart, kend):
        for j in range(jstart, jend):
            for i in range(istart, iend):
                ijk = i + j*jj + k*kk
    
                st[ijk] += \
                    - (  u[ijk+ii] * 0.5*(s[ijk   ] + s[ijk+ii])
                       - u[ijk   ] * 0.5*(s[ijk-ii] + s[ijk   ]) ) * dxi \
                    - (  v[ijk+jj] * 0.5*(s[ijk   ] + s[ijk+jj])
                       - v[ijk   ] * 0.5*(s[ijk-jj] + s[ijk   ]) ) * dyi \
                    - (  rhorefh[k+1] * w[ijk+kk] * 0.5*(s[ijk   ] + s[ijk+kk])
                       - rhorefh[k  ] * w[ijk   ] * 0.5*(s[ijk-kk] + s[ijk   ]) ) / rhoref[k] * dzi[k];


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

        for manual_check in [True, False]:
    
            grid       = Grid(3200, 3200, 3200, gridsize, gridsize, gridsize, 2, 1, TF)
            fields     = Fields(['u','v','w','s'], grid.ncells, grid.ijcells, grid.kcells, TF)
            fields_ref = Fields(['u','v','w','s'], grid.ncells, grid.ijcells, grid.kcells, TF)

            # Argument list of CUDA kernel
            args = [
                fields.s.tend, fields.s.fld, fields.u.fld, fields.v.fld, fields.w.fld, 
                fields.rhoref, fields.rhorefh, grid.dzi, 
                grid.dxi, grid.dyi, 
                grid.icells, grid.ijcells, 
                grid.istart, grid.jstart, grid.kstart, 
                grid.iend, grid.jend, grid.kend]

            # Validate results with Python version of code
            advec_s(
                fields_ref.s.tend, fields_ref.s.fld, fields_ref.u.fld, fields_ref.v.fld, fields_ref.w.fld, 
                fields_ref.rhoref, fields_ref.rhorefh, grid.dzi, 
                grid.dxi, grid.dyi, 
                grid.icells, grid.ijcells, 
                grid.istart, grid.jstart, grid.kstart, 
                grid.iend, grid.jend, grid.kend)

            if manual_check:

                # CUDA
                params = { "block_size_x": 4, "block_size_y": 4 }
                results = kernel_tuner.run_kernel(
                    "advec_s_g", kernel_string, (grid.itot, grid.jtot, grid.ktot), args, params)

                # Compare results
                if np.allclose(results[0], fields_ref.s.tend, atol=1e-14):
                    print('SUCCESS, CUDA and Python results are identical!')
                else:
                    print('ERROR: CUDA and Python results not identical')
                    print('Max diff = {}'.format(np.absolute(results[0]-fields_ref.s.tend).max()))

            else:

                # Tune parameters
                tune_params = OrderedDict()
                tune_params["block_size_x"] = [32*i for i in range(1,9)]
                tune_params["block_size_y"] = [2**i for i in range(3)]

                # True answer from the Python code
                answer = len(args)*[None]
                answer[0] = fields_ref.s.tend

                kernel_tuner.tune_kernel(
                    "advec_s_g", kernel_string, (grid.itot, grid.jtot, grid.ktot), args, tune_params, answer=answer, atol=1e-14)
