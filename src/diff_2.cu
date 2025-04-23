/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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
#include "diff_2.h"
#include "stats.h"
#include "defines.h"
#include "constants.h"
#include "tools.h"

namespace
{
    template<typename TF> __global__
    void diff_c_g(TF* __restrict__ const at, const TF* __restrict__ const a,
                  const TF* __restrict__ const dzi, const TF* __restrict__ const dzhi,
                  const TF dxidxi, const TF dyidyi, const TF visc,
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
            const int ii = 1;

            at[ijk] += visc * (
                + (  (a[ijk+ii] - a[ijk   ])
                   - (a[ijk   ] - a[ijk-ii]) ) * dxidxi
                + (  (a[ijk+jj] - a[ijk   ])
                   - (a[ijk   ] - a[ijk-jj]) ) * dyidyi
                + (  (a[ijk+kk] - a[ijk   ]) * dzhi[k+1]
                   - (a[ijk   ] - a[ijk-kk]) * dzhi[k]   ) * dzi[k]);
        }
    }

    template<typename TF> __global__
    void diff_w_g(TF* __restrict__ const at, const TF* __restrict__ const a,
                  const TF* __restrict__ const dzi, const TF* __restrict__ const dzhi,
                  const TF dxidxi, const TF dyidyi, const TF visc,
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
            const int ii = 1;

            at[ijk] += visc * (
                + (  (a[ijk+ii] - a[ijk   ])
                    - (a[ijk   ] - a[ijk-ii]) ) * dxidxi
                + (  (a[ijk+jj] - a[ijk   ])
                    - (a[ijk   ] - a[ijk-jj]) ) * dyidyi
                + (  (a[ijk+kk] - a[ijk   ]) * dzi[k]
                    - (a[ijk   ] - a[ijk-kk]) * dzi[k-1] ) * dzhi[k]);
        }
    }


    template<typename TF>__global__
    void calc_diff_flux_g(TF* __restrict__ st, const TF* const __restrict__ s, const TF visc, const TF* __restrict__ const dzhi,
                   const int icells, const int jcells, const int kstart, const int kend, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;
        const int k = blockIdx.z + kstart;

        if (i < icells && j < jcells && k < kend)
        {
            const int ijk = i + j*icells + k*ijcells;
            st[ijk] = - visc*(s[ijk] - s[ijk-ijcells])*dzhi[k];
        }
    }

}

#ifdef USECUDA
template<typename TF>
void Diff_2<TF>::exec(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const TF dxidxi = 1./(gd.dx*gd.dx);
    const TF dyidyi = 1./(gd.dy*gd.dy);


    diff_c_g<TF><<<gridGPU, blockGPU>>>(
        fields.mt.at("u")->fld_g, fields.mp.at("u")->fld_g,
        gd.dzi_g, gd.dzhi_g,
        dxidxi, dyidyi, fields.visc,
        gd.icells, gd.ijcells,
        gd.istart,  gd.jstart, gd.kstart,
        gd.iend,    gd.jend,   gd.kend);
    cuda_check_error();

    diff_c_g<TF><<<gridGPU, blockGPU>>>(
        fields.mt.at("v")->fld_g, fields.mp.at("v")->fld_g,
        gd.dzi_g, gd.dzhi_g,
        dxidxi, dyidyi, fields.visc,
        gd.icells, gd.ijcells,
        gd.istart,  gd.jstart, gd.kstart,
        gd.iend,    gd.jend,   gd.kend);
    cuda_check_error();

    diff_w_g<TF><<<gridGPU, blockGPU>>>(
        fields.mt.at("w")->fld_g, fields.mp.at("w")->fld_g,
        gd.dzi_g, gd.dzhi_g,
        dxidxi, dyidyi, fields.visc,
        gd.icells, gd.ijcells,
        gd.istart,  gd.jstart, gd.kstart,
        gd.iend,    gd.jend,   gd.kend);
    cuda_check_error();

    for (auto &it : fields.st)
        diff_c_g<TF><<<gridGPU, blockGPU>>>(
            it.second->fld_g, fields.sp.at(it.first)->fld_g,
            gd.dzi_g, gd.dzhi_g,
            dxidxi, dyidyi, fields.sp.at(it.first)->visc,
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

template<typename TF>
void Diff_2<TF>::get_diff_flux(
        Field3d<TF>& diff_flux, const Field3d<TF>& fld)
{
    using namespace Tools_g;

    auto& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);
    const int ngrid  = gd.ncells/blocki + (gd.ncells%blocki > 0);
    const int nblock = gd.ithread_block;
    dim3 gridGPU(gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);
    calc_diff_flux_g<TF><<<gridGPU, blockGPU>>>(
        diff_flux.fld_g.data(), fld.fld_g.data(), fld.visc, gd.dzhi_g.data(), gd.icells, gd.jcells, gd.kstart, gd.kend, gd.ijcells);
}

#endif


#ifdef FLOAT_SINGLE
template class Diff_2<float>;
#else
template class Diff_2<double>;
#endif
