/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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
#include "stats.h"
#include "finite_difference.h"
#include "field3d_operators.h"

using namespace Finite_difference::O2;

namespace
{
    template<typename TF>__global__
    void advec_uvw_g(TF* __restrict__ ut, TF* __restrict__ vt, TF * __restrict__ wt,
                     TF* __restrict__ u,  TF* __restrict__ v,  TF * __restrict__ w,
                     const TF* __restrict__ rhoref, const TF* __restrict__ rhorefh,
                     const TF* __restrict__ dzi,    const TF* __restrict__ dzhi, TF dxi, TF dyi,
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

    template<typename TF>__global__
    void advec_s_g(TF* __restrict__ st, TF* __restrict__ s,
                   TF* __restrict__ u,  TF* __restrict__ v, TF* __restrict__ w,
                   const TF* __restrict__ rhoref, const TF* __restrict__ rhorefh,
                   const TF* __restrict__ dzi, TF dxi, TF dyi,
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

    template<typename TF>__global__
    void calc_cfl_g(TF* __restrict__ u, TF* __restrict__ v, TF* __restrict__ w,
                    TF* __restrict__ cfl, const TF* __restrict__ dzi, TF dxi, TF dyi,
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
template<typename TF>
unsigned long Advec_2<TF>::get_time_limit(unsigned long idt, double dt)
{
    // Calculate cfl and prevent zero divisons.
    double cfl = get_cfl(dt);
    cfl = std::max(cflmin, cfl);
    const unsigned long idtlim = idt * cflmax / cfl;

    return idtlim;
}

template<typename TF>
double Advec_2<TF>::get_cfl(const double dt)
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    const TF dxi = 1./gd.dx;
    const TF dyi = 1./gd.dy;

    auto tmp1 = fields.get_tmp_g();

    calc_cfl_g<TF><<<gridGPU, blockGPU>>>(
        fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
        tmp1->fld_g, gd.dzi_g, dxi, dyi,
        gd.icells, gd.ijcells,
        gd.istart,  gd.jstart, gd.kstart,
        gd.iend,    gd.jend,   gd.kend);
    cuda_check_error();

    TF cfl = field3d_operators.calc_max_g(tmp1->fld_g);
    fields.release_tmp_g(tmp1);

    cfl = cfl*dt;

    return static_cast<double>(cfl);
}

template<typename TF>
void Advec_2<TF>::exec(Stats<TF>& stats)
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const TF dxi = 1./gd.dx;
    const TF dyi = 1./gd.dy;

    advec_uvw_g<TF><<<gridGPU, blockGPU>>>(
        fields.mt.at("u")->fld_g, fields.mt.at("v")->fld_g, fields.mt.at("w")->fld_g,
        fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
        fields.rhoref_g, fields.rhorefh_g, gd.dzi_g, gd.dzhi_g, dxi, dyi,
        gd.icells, gd.ijcells,
        gd.istart,  gd.jstart, gd.kstart,
        gd.iend,    gd.jend,   gd.kend);
    cuda_check_error();

    for (auto& it : fields.st)
        advec_s_g<TF><<<gridGPU, blockGPU>>>(
            it.second->fld_g, fields.sp.at(it.first)->fld_g,
            fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
            fields.rhoref_g, fields.rhorefh_g, gd.dzi_g, dxi, dyi,
            gd.icells, gd.ijcells,
            gd.istart,  gd.jstart, gd.kstart,
            gd.iend,    gd.jend,   gd.kend);
    cuda_check_error();

    cudaDeviceSynchronize();
    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);
    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);
}
#endif


#ifdef FLOAT_SINGLE
template class Advec_2<float>;
#else
template class Advec_2<double>;
#endif
