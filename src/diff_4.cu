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

#include "grid.h"
#include "fields.h"
#include "diff_4.h"
#include "defines.h"
#include "stats.h"
#include "constants.h"
#include "tools.h"
#include "finite_difference.h"

using namespace Finite_difference::O4;

namespace
{
    template<typename TF> __global__
    void diff_c_g(TF* __restrict__ const at, const TF* __restrict__ const a,
                  const TF* __restrict__ const dzi4, const TF* __restrict__ const dzhi4,
                  const TF dx, const TF dy, const TF visc,
                  const int jj,     const int kk,
                  const int istart, const int jstart, const int kstart,
                  const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

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

            const TF dxidxi = TF(1.)/(dx*dx);
            const TF dyidyi = TF(1.)/(dy*dy);

            // bottom boundary
            if (k == kstart)
            {
                at[ijk] += visc * (cdg3<TF>*a[ijk-ii3] + cdg2<TF>*a[ijk-ii2] + cdg1<TF>*a[ijk-ii1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+ii1] + cdg2<TF>*a[ijk+ii2] + cdg3<TF>*a[ijk+ii3])*dxidxi;
                at[ijk] += visc * (cdg3<TF>*a[ijk-jj3] + cdg2<TF>*a[ijk-jj2] + cdg1<TF>*a[ijk-jj1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+jj1] + cdg2<TF>*a[ijk+jj2] + cdg3<TF>*a[ijk+jj3])*dyidyi;
                at[ijk] += visc * ( cg0<TF>*(bg0<TF>*a[ijk-kk2] + bg1<TF>*a[ijk-kk1] + bg2<TF>*a[ijk    ] + bg3<TF>*a[ijk+kk1]) * dzhi4[k-1]
                            + cg1<TF>*(cg0<TF>*a[ijk-kk2] + cg1<TF>*a[ijk-kk1] + cg2<TF>*a[ijk    ] + cg3<TF>*a[ijk+kk1]) * dzhi4[k  ]
                            + cg2<TF>*(cg0<TF>*a[ijk-kk1] + cg1<TF>*a[ijk    ] + cg2<TF>*a[ijk+kk1] + cg3<TF>*a[ijk+kk2]) * dzhi4[k+1]
                            + cg3<TF>*(cg0<TF>*a[ijk    ] + cg1<TF>*a[ijk+kk1] + cg2<TF>*a[ijk+kk2] + cg3<TF>*a[ijk+kk3]) * dzhi4[k+2] ) * dzi4[k];
            }
            // top boundary
            else if (k == kend-1)
            {
                at[ijk] += visc * (cdg3<TF>*a[ijk-ii3] + cdg2<TF>*a[ijk-ii2] + cdg1<TF>*a[ijk-ii1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+ii1] + cdg2<TF>*a[ijk+ii2] + cdg3<TF>*a[ijk+ii3])*dxidxi;
                at[ijk] += visc * (cdg3<TF>*a[ijk-jj3] + cdg2<TF>*a[ijk-jj2] + cdg1<TF>*a[ijk-jj1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+jj1] + cdg2<TF>*a[ijk+jj2] + cdg3<TF>*a[ijk+jj3])*dyidyi;
                at[ijk] += visc * ( cg0<TF>*(cg0<TF>*a[ijk-kk3] + cg1<TF>*a[ijk-kk2] + cg2<TF>*a[ijk-kk1] + cg3<TF>*a[ijk    ]) * dzhi4[k-1]
                             + cg1<TF>*(cg0<TF>*a[ijk-kk2] + cg1<TF>*a[ijk-kk1] + cg2<TF>*a[ijk    ] + cg3<TF>*a[ijk+kk1]) * dzhi4[k  ]
                             + cg2<TF>*(cg0<TF>*a[ijk-kk1] + cg1<TF>*a[ijk    ] + cg2<TF>*a[ijk+kk1] + cg3<TF>*a[ijk+kk2]) * dzhi4[k+1]
                             + cg3<TF>*(tg0<TF>*a[ijk-kk1] + tg1<TF>*a[ijk    ] + tg2<TF>*a[ijk+kk1] + tg3<TF>*a[ijk+kk2]) * dzhi4[k+2] ) * dzi4[k];
            }
            // interior
            else
            {
                at[ijk] += visc * (cdg3<TF>*a[ijk-ii3] + cdg2<TF>*a[ijk-ii2] + cdg1<TF>*a[ijk-ii1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+ii1] + cdg2<TF>*a[ijk+ii2] + cdg3<TF>*a[ijk+ii3])*dxidxi;
                at[ijk] += visc * (cdg3<TF>*a[ijk-jj3] + cdg2<TF>*a[ijk-jj2] + cdg1<TF>*a[ijk-jj1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+jj1] + cdg2<TF>*a[ijk+jj2] + cdg3<TF>*a[ijk+jj3])*dyidyi;
                at[ijk] += visc * ( cg0<TF>*(cg0<TF>*a[ijk-kk3] + cg1<TF>*a[ijk-kk2] + cg2<TF>*a[ijk-kk1] + cg3<TF>*a[ijk    ]) * dzhi4[k-1]
                             + cg1<TF>*(cg0<TF>*a[ijk-kk2] + cg1<TF>*a[ijk-kk1] + cg2<TF>*a[ijk    ] + cg3<TF>*a[ijk+kk1]) * dzhi4[k  ]
                             + cg2<TF>*(cg0<TF>*a[ijk-kk1] + cg1<TF>*a[ijk    ] + cg2<TF>*a[ijk+kk1] + cg3<TF>*a[ijk+kk2]) * dzhi4[k+1]
                             + cg3<TF>*(cg0<TF>*a[ijk    ] + cg1<TF>*a[ijk+kk1] + cg2<TF>*a[ijk+kk2] + cg3<TF>*a[ijk+kk3]) * dzhi4[k+2] ) * dzi4[k];
            }
        }
    }

    template<typename TF> __global__
    void diff_w_g(TF* __restrict__ const at, const TF* __restrict__ const a,
                  const TF* __restrict__ const dzi4, const TF* __restrict__ const dzhi4,
                  const TF dx, const TF dy, const TF visc,
                  const int jj,     const int kk,
                  const int istart, const int jstart, const int kstart,
                  const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k > kstart && k < kend)
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

            const TF dxidxi = 1./(dx*dx);
            const TF dyidyi = 1./(dy*dy);

            if (k == kstart+1)
            {
                // bottom boundary
                at[ijk] += visc * (cdg3<TF>*a[ijk-ii3] + cdg2<TF>*a[ijk-ii2] + cdg1<TF>*a[ijk-ii1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+ii1] + cdg2<TF>*a[ijk+ii2] + cdg3<TF>*a[ijk+ii3])*dxidxi;
                at[ijk] += visc * (cdg3<TF>*a[ijk-jj3] + cdg2<TF>*a[ijk-jj2] + cdg1<TF>*a[ijk-jj1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+jj1] + cdg2<TF>*a[ijk+jj2] + cdg3<TF>*a[ijk+jj3])*dyidyi;
                at[ijk] += visc * ( cg0<TF>*(bg0<TF>*a[ijk-kk2] + bg1<TF>*a[ijk-kk1] + bg2<TF>*a[ijk    ] + bg3<TF>*a[ijk+kk1]) * dzi4[k-1]
                             + cg1<TF>*(cg0<TF>*a[ijk-kk2] + cg1<TF>*a[ijk-kk1] + cg2<TF>*a[ijk    ] + cg3<TF>*a[ijk+kk1]) * dzi4[k  ]
                             + cg2<TF>*(cg0<TF>*a[ijk-kk1] + cg1<TF>*a[ijk    ] + cg2<TF>*a[ijk+kk1] + cg3<TF>*a[ijk+kk2]) * dzi4[k+1]
                             + cg3<TF>*(cg0<TF>*a[ijk    ] + cg1<TF>*a[ijk+kk1] + cg2<TF>*a[ijk+kk2] + cg3<TF>*a[ijk+kk3]) * dzi4[k+2] ) * dzhi4[k];
            }
            else if (k == kend-1)
            {
                // top boundary
                at[ijk] += visc * (cdg3<TF>*a[ijk-ii3] + cdg2<TF>*a[ijk-ii2] + cdg1<TF>*a[ijk-ii1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+ii1] + cdg2<TF>*a[ijk+ii2] + cdg3<TF>*a[ijk+ii3])*dxidxi;
                at[ijk] += visc * (cdg3<TF>*a[ijk-jj3] + cdg2<TF>*a[ijk-jj2] + cdg1<TF>*a[ijk-jj1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+jj1] + cdg2<TF>*a[ijk+jj2] + cdg3<TF>*a[ijk+jj3])*dyidyi;
                at[ijk] += visc * ( cg0<TF>*(cg0<TF>*a[ijk-kk3] + cg1<TF>*a[ijk-kk2] + cg2<TF>*a[ijk-kk1] + cg3<TF>*a[ijk    ]) * dzi4[k-2]
                             + cg1<TF>*(cg0<TF>*a[ijk-kk2] + cg1<TF>*a[ijk-kk1] + cg2<TF>*a[ijk    ] + cg3<TF>*a[ijk+kk1]) * dzi4[k-1]
                             + cg2<TF>*(cg0<TF>*a[ijk-kk1] + cg1<TF>*a[ijk    ] + cg2<TF>*a[ijk+kk1] + cg3<TF>*a[ijk+kk2]) * dzi4[k  ]
                             + cg3<TF>*(tg0<TF>*a[ijk-kk1] + tg1<TF>*a[ijk    ] + tg2<TF>*a[ijk+kk1] + tg3<TF>*a[ijk+kk2]) * dzi4[k+1] ) * dzhi4[k];
            }
            else
            {
                // interior
                at[ijk] += visc * (cdg3<TF>*a[ijk-ii3] + cdg2<TF>*a[ijk-ii2] + cdg1<TF>*a[ijk-ii1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+ii1] + cdg2<TF>*a[ijk+ii2] + cdg3<TF>*a[ijk+ii3])*dxidxi;
                at[ijk] += visc * (cdg3<TF>*a[ijk-jj3] + cdg2<TF>*a[ijk-jj2] + cdg1<TF>*a[ijk-jj1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+jj1] + cdg2<TF>*a[ijk+jj2] + cdg3<TF>*a[ijk+jj3])*dyidyi;
                at[ijk] += visc * ( cg0<TF>*(cg0<TF>*a[ijk-kk3] + cg1<TF>*a[ijk-kk2] + cg2<TF>*a[ijk-kk1] + cg3<TF>*a[ijk    ]) * dzi4[k-2]
                             + cg1<TF>*(cg0<TF>*a[ijk-kk2] + cg1<TF>*a[ijk-kk1] + cg2<TF>*a[ijk    ] + cg3<TF>*a[ijk+kk1]) * dzi4[k-1]
                             + cg2<TF>*(cg0<TF>*a[ijk-kk1] + cg1<TF>*a[ijk    ] + cg2<TF>*a[ijk+kk1] + cg3<TF>*a[ijk+kk2]) * dzi4[k  ]
                             + cg3<TF>*(cg0<TF>*a[ijk    ] + cg1<TF>*a[ijk+kk1] + cg2<TF>*a[ijk+kk2] + cg3<TF>*a[ijk+kk3]) * dzi4[k+1] ) * dzhi4[k];
            }
        }
    }
}

#ifdef USECUDA
template<typename TF>
void Diff_4<TF>::exec(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    const int blocki = 128;
    const int blockj = 2;
    const int gridi = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    diff_c_g<TF><<<gridGPU, blockGPU>>>(
            fields.mt.at("u")->fld_g, fields.mp.at("u")->fld_g,
            gd.dzi4_g, gd.dzhi4_g,
            gd.dx, gd.dy, fields.visc,
            gd.icells, gd.ijcells,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    diff_c_g<TF><<<gridGPU, blockGPU>>>(
            fields.mt.at("v")->fld_g, fields.mp.at("v")->fld_g,
            gd.dzi4_g, gd.dzhi4_g,
            gd.dx, gd.dy, fields.visc,
            gd.icells, gd.ijcells,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    diff_w_g<TF><<<gridGPU, blockGPU>>>(
            fields.mt.at("w")->fld_g, fields.mp.at("w")->fld_g,
            gd.dzi4_g, gd.dzhi4_g,
            gd.dx, gd.dy, fields.visc,
            gd.icells, gd.ijcells,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    for (auto& it : fields.st)
        diff_c_g<TF><<<gridGPU, blockGPU>>>(
                it.second->fld_g, fields.sp.at(it.first)->fld_g,
                gd.dzi4_g, gd.dzhi4_g,
                gd.dx, gd.dy, fields.sp.at(it.first)->visc,
                gd.icells, gd.ijcells,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend);
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
template class Diff_4<float>;
#else
template class Diff_4<double>;
#endif
