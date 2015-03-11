/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "advec_2.h"
#include "grid.h"
#include "fields.h"
#include "tools.h"
#include "constants.h"
#include "tools.h"
#include "finite_difference.h"

using namespace Finite_difference::O2;

namespace
{
    __global__ 
    void advec_uvw_g(double* __restrict__ ut, double* __restrict__ vt, double * __restrict__ wt, 
                     double* __restrict__ u,  double* __restrict__ v,  double * __restrict__ w,
                     double* __restrict__ rhoref, double* __restrict__ rhorefh,
                     double* __restrict__ dzi,    double* __restrict__ dzhi, double dxi, double dyi, 
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
            ut[ijk] += 
                - (  interp2(u[ijk   ], u[ijk+ii]) * interp2(u[ijk   ], u[ijk+ii])
                   - interp2(u[ijk-ii], u[ijk   ]) * interp2(u[ijk-ii], u[ijk   ]) ) * dxi

                - (  interp2(v[ijk-ii+jj], v[ijk+jj]) * interp2(u[ijk   ], u[ijk+jj])
                   - interp2(v[ijk-ii   ], v[ijk   ]) * interp2(u[ijk-jj], u[ijk   ]) ) * dyi 

                - (  rhorefh[k+1] * interp2(w[ijk-ii+kk], w[ijk+kk]) * interp2(u[ijk   ], u[ijk+kk])
                   - rhorefh[k  ] * interp2(w[ijk-ii   ], w[ijk   ]) * interp2(u[ijk-kk], u[ijk   ]) ) / rhoref[k] * dzi[k];

            vt[ijk] += 
                - (  interp2(u[ijk+ii-jj], u[ijk+ii]) * interp2(v[ijk   ], v[ijk+ii])
                   - interp2(u[ijk   -jj], u[ijk   ]) * interp2(v[ijk-ii], v[ijk   ]) ) * dxi

                - (  interp2(v[ijk   ], v[ijk+jj]) * interp2(v[ijk   ], v[ijk+jj])
                   - interp2(v[ijk-jj], v[ijk   ]) * interp2(v[ijk-jj], v[ijk   ]) ) * dyi

                - (  rhorefh[k+1] * interp2(w[ijk-jj+kk], w[ijk+kk]) * interp2(v[ijk   ], v[ijk+kk])
                   - rhorefh[k  ] * interp2(w[ijk-jj   ], w[ijk   ]) * interp2(v[ijk-kk], v[ijk   ]) ) / rhoref[k] * dzi[k];

            if (k > kstart)
            {
                wt[ijk] += 
                    - (  interp2(u[ijk+ii-kk], u[ijk+ii]) * interp2(w[ijk   ], w[ijk+ii])
                       - interp2(u[ijk   -kk], u[ijk   ]) * interp2(w[ijk-ii], w[ijk   ]) ) * dxi

                    - (  interp2(v[ijk+jj-kk], v[ijk+jj]) * interp2(w[ijk   ], w[ijk+jj])
                       - interp2(v[ijk   -kk], v[ijk   ]) * interp2(w[ijk-jj], w[ijk   ]) ) * dyi

                    - (  rhoref[k  ] * interp2(w[ijk   ], w[ijk+kk]) * interp2(w[ijk   ], w[ijk+kk])
                       - rhoref[k-1] * interp2(w[ijk-kk], w[ijk   ]) * interp2(w[ijk-kk], w[ijk   ]) ) / rhorefh[k] * dzhi[k];
            }
        }
    }

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
                - (  u[ijk+ii] * interp2(s[ijk   ], s[ijk+ii])
                   - u[ijk   ] * interp2(s[ijk-ii], s[ijk   ]) ) * dxi

                - (  v[ijk+jj] * interp2(s[ijk   ], s[ijk+jj])
                   - v[ijk   ] * interp2(s[ijk-jj], s[ijk   ]) ) * dyi 

                - (  rhorefh[k+1] * w[ijk+kk] * interp2(s[ijk   ], s[ijk+kk])
                   - rhorefh[k  ] * w[ijk   ] * interp2(s[ijk-kk], s[ijk   ]) ) / rhoref[k] * dzi[k];
        }
    }

    __global__ 
    void calc_cfl_g(double* __restrict__ u, double* __restrict__ v, double* __restrict__ w, 
                    double* __restrict__ cfl, double* __restrict__ dzi, double dxi, double dyi,
                    int jj, int kk,
                    int istart, int jstart, int kstart,
                    int iend, int jend, int kend)
    {
        const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart; 
        const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
        const int k  = blockIdx.z + kstart; 
        const int ii = 1;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            cfl[ijk] = std::abs(interp2(u[ijk], u[ijk+ii]))*dxi + 
                       std::abs(interp2(v[ijk], v[ijk+jj]))*dyi + 
                       std::abs(interp2(w[ijk], w[ijk+kk]))*dzi[k];
        }
    }
}

#ifdef USECUDA
unsigned long Advec_2::get_time_limit(unsigned long idt, double dt)
{
    // Calculate cfl and prevent zero divisons.
    double cfl = get_cfl(dt);
    cfl = std::max(cflmin, cfl);
    const unsigned long idtlim = idt * cflmax / cfl;

    return idtlim;
}

double Advec_2::get_cfl(const double dt)
{
    const int blocki = grid->ithread_block;
    const int blockj = grid->jthread_block;
    const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kcells);
    dim3 blockGPU(blocki, blockj, 1);

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    const int offs = grid->memoffset;

    calc_cfl_g<<<gridGPU, blockGPU>>>(
        &fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs], 
        &fields->atmp["tmp1"]->data_g[offs], grid->dzi_g, dxi, dyi,
        grid->icellsp, grid->ijcellsp,
        grid->istart,  grid->jstart, grid->kstart,
        grid->iend,    grid->jend,   grid->kend);
    cuda_check_error(); 

    double cfl = grid->get_max_g(&fields->atmp["tmp1"]->data_g[offs], fields->atmp["tmp2"]->data_g); 
    grid->get_max(&cfl); 
    cfl = cfl*dt;

    return cfl;
}

void Advec_2::exec()
{
    const int blocki = grid->ithread_block;
    const int blockj = grid->jthread_block;
    const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    const int offs = grid->memoffset;

    advec_uvw_g<<<gridGPU, blockGPU>>>(
        &fields->ut->data_g[offs], &fields->vt->data_g[offs], &fields->wt->data_g[offs], 
        &fields->u->data_g[offs],  &fields->v->data_g[offs],  &fields->w->data_g[offs], 
        fields->rhoref_g, fields->rhorefh_g, grid->dzi_g, grid->dzhi_g, dxi, dyi, 
        grid->icellsp, grid->ijcellsp,
        grid->istart,  grid->jstart, grid->kstart,
        grid->iend,    grid->jend,   grid->kend);
    cuda_check_error(); 

    for (FieldMap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
        advec_s_g<<<gridGPU, blockGPU>>>(
            &it->second->data_g[offs], &fields->sp[it->first]->data_g[offs], 
            &fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs],
            fields->rhoref_g, fields->rhorefh_g, grid->dzi_g, dxi, dyi,
            grid->icellsp, grid->ijcellsp,
            grid->istart,  grid->jstart, grid->kstart,
            grid->iend,    grid->jend,   grid->kend);
    cuda_check_error(); 
}
#endif
