/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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

#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "thermo_buoy.h"
#include "master.h"
#include "finite_difference.h"
#include "tools.h"
#include "stats.h"

namespace
{
    template<typename TF> __global__
    void calc_buoyancy_g(TF* __restrict__ b,
                         TF* __restrict__ bin)
    {
        b[threadIdx.x] = bin[threadIdx.x];
    }

    template<typename TF> __global__
    void calc_buoyancy_bot_g(TF* __restrict__ b,     TF* __restrict__ bbot,
                             TF* __restrict__ bin,    TF* __restrict__ bbotin,
                             int kstart, int icells, int jcells,
                             int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            bbot[ij] = bbotin[ij];
            b[ijk]   = bin[ijk];
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_flux_bot_g(TF* __restrict__ bfluxbot, TF* __restrict__ bfluxbotin,
                                  int kstart, int icells, int jcells,
                                  int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            bfluxbot[ij] = bfluxbotin[ij];
        }
    }

    template<typename TF> __global__
    void calc_N2_g(TF* __restrict__ N2,    TF* __restrict__ b,
                   const TF bg_n2, TF* __restrict__ dzi,
                   int istart, int jstart, int kstart,
                   int iend,   int jend,   int kend,
                   int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            N2[ijk] = static_cast<TF>(0.5)*(b[ijk+kk] - b[ijk-kk])*dzi[k] + bg_n2;
        }
    }


    template<typename TF> __global__
    void calc_buoyancy_tend_2nd_g(TF* __restrict__ wt, TF* __restrict__ b,
                                  int istart, int jstart, int kstart,
                                  int iend,   int jend,   int kend,
                                  int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        using Finite_difference::O2::interp2;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            wt[ijk] += interp2(b[ijk-kk], b[ijk]);
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_tend_u_2nd_g(TF* const __restrict__ ut, const TF* const __restrict__ b,
                                    const TF sinalpha,
                                    const int istart, const int jstart, const int kstart,
                                    const int iend,   const int jend,   const int kend,
                                    const int jj,     const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int ii1 = 1;

        using Finite_difference::O2::interp2;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            ut[ijk] += sinalpha * interp2(b[ijk-ii1], b[ijk]);
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_tend_w_2nd_g(TF* __restrict__ wt, const TF* const __restrict__ b,
                                    const TF cosalpha,
                                    const int istart, const int jstart, const int kstart,
                                    const int iend, const int jend, const int kend,
                                    const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int kk1 = 1*kk;

        using Finite_difference::O2::interp2;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            wt[ijk] += cosalpha * interp2(b[ijk-kk1], b[ijk]);
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_tend_b_2nd_g(TF* const __restrict__ bt,
                                    const TF* const __restrict__ u, const TF* const __restrict__ w,
                                    const TF utrans, const TF n2, const TF sinalpha, const TF cosalpha,
                                    const int istart, const int jstart, const int kstart,
                                    const int iend, const int jend, const int kend,
                                    const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int ii1 = 1;
        const int kk1 = 1*kk;

        using Finite_difference::O2::interp2;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            bt[ijk] -= n2 * ( sinalpha * ( interp2(u[ijk], u[ijk+ii1]) + utrans )
                            + cosalpha * ( interp2(w[ijk], w[ijk+kk1]) ) );
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_tend_4th_g(TF* __restrict__ wt, TF* __restrict__ b,
                                  int istart, int jstart, int kstart,
                                  int iend,   int jend,   int kend,
                                  int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        using namespace Finite_difference::O4;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            wt[ijk] += ci0<TF>*b[ijk-kk2] + ci1<TF>*b[ijk-kk1] + ci2<TF>*b[ijk] + ci3<TF>*b[ijk+kk1];
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_tend_u_4th_g(TF* const __restrict__ ut, const TF* const __restrict__ b,
                                    const TF sinalpha,
                                    const int istart, const int jstart, const int kstart,
                                    const int iend,   const int jend,   const int kend,
                                    const int jj,     const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int ii1 = 1;
        const int ii2 = 2;

        using namespace Finite_difference::O4;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            ut[ijk] += sinalpha * (ci0<TF>*b[ijk-ii2] + ci1<TF>*b[ijk-ii1] + ci2<TF>*b[ijk] + ci3<TF>*b[ijk+ii1]);
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_tend_w_4th_g(TF* __restrict__ wt, const TF* const __restrict__ b,
                                    const TF cosalpha,
                                    const int istart, const int jstart, const int kstart,
                                    const int iend, const int jend, const int kend,
                                    const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        using namespace Finite_difference::O4;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            wt[ijk] += cosalpha * (ci0<TF>*b[ijk-kk2] + ci1<TF>*b[ijk-kk1] + ci2<TF>*b[ijk] + ci3<TF>*b[ijk+kk1]);
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_tend_b_4th_g(TF* const __restrict__ bt,
                                    const TF* const __restrict__ u, const TF* const __restrict__ w,
                                    const TF utrans, const TF n2, const TF sinalpha, const TF cosalpha,
                                    const int istart, const int jstart, const int kstart,
                                    const int iend, const int jend, const int kend,
                                    const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int ii1 = 1;
        const int ii2 = 2;

        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        using namespace Finite_difference::O4;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            bt[ijk] -= n2 * ( sinalpha * ( (ci0<TF>*u[ijk-ii1] + ci1<TF>*u[ijk] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii2]) + utrans )
                            + cosalpha * (  ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2]) );
        }
    }

    template<typename TF> __global__
    void calc_baroclinic_2nd_g(TF* __restrict__ bt, const TF* __restrict__ v,
                               const TF dbdy_ls,
                               int istart, int jstart, int kstart,
                               int iend,   int jend,   int kend,
                               int jj, int kk)
    {
        using Finite_difference::O2::interp2;

        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            bt[ijk] -= dbdy_ls * interp2(v[ijk], v[ijk+jj]);
        }
    }

    template<typename TF> __global__
    void calc_baroclinic_4th_g(TF* __restrict__ bt, const TF* __restrict__ v,
                               const TF dbdy_ls,
                               int istart, int jstart, int kstart,
                               int iend,   int jend,   int kend,
                               int jj, int kk)
    {
        using Finite_difference::O4::interp4c;

        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int jj1 = 1*jj;
        const int jj2 = 2*jj;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            bt[ijk] -= dbdy_ls * interp4c(v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2]);
        }
    }
} // End namespace.

#ifdef USECUDA
template<typename TF>
void Thermo_buoy<TF>::exec(const double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax-1);
    dim3 blockGPU(blocki, blockj, 1);

    if (grid.get_spatial_order() == Grid_order::Second)
    {
        if (bs.has_slope || bs.has_N2)
        {
            const TF sinalpha = std::sin(bs.alpha);
            const TF cosalpha = std::cos(bs.alpha);

            calc_buoyancy_tend_u_2nd_g<<<gridGPU, blockGPU>>>(
                fields.mt.at("u")->fld_g, fields.sp.at("b")->fld_g,
                sinalpha,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();

            calc_buoyancy_tend_w_2nd_g<<<gridGPU, blockGPU>>>(
                fields.mt.at("w")->fld_g, fields.sp.at("b")->fld_g,
                cosalpha,
                gd.istart, gd.jstart, gd.kstart+1,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();

            calc_buoyancy_tend_b_2nd_g<<<gridGPU, blockGPU>>>(
                fields.st.at("b")->fld_g,
                fields.mp.at("u")->fld_g, fields.mp.at("w")->fld_g,
                gd.utrans, bs.n2, sinalpha, cosalpha,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();
        }
        else
        {
            calc_buoyancy_tend_2nd_g<<<gridGPU, blockGPU>>>(
            fields.mt.at("w")->fld_g, fields.sp.at("b")->fld_g,
            gd.istart, gd.jstart, gd.kstart+1,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
            cuda_check_error();
        }

        if (swbaroclinic)
        {
            calc_baroclinic_2nd_g<<<gridGPU, blockGPU>>>(
                fields.st.at("th")->fld_g, fields.mp.at("v")->fld_g,
                dbdy_ls,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();
        }
    }
    else if (grid.get_spatial_order() == Grid_order::Fourth)
    {
        if (bs.has_slope || bs.has_N2)
        {
            const TF sinalpha = std::sin(bs.alpha);
            const TF cosalpha = std::cos(bs.alpha);
            calc_buoyancy_tend_u_4th_g<<<gridGPU, blockGPU>>>(
                fields.mt.at("u")->fld_g, fields.sp.at("b")->fld_g,
                sinalpha,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();

            calc_buoyancy_tend_w_4th_g<<<gridGPU, blockGPU>>>(
                fields.mt.at("w")->fld_g, fields.sp.at("b")->fld_g,
                cosalpha,
                gd.istart, gd.jstart, gd.kstart+1,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();

            calc_buoyancy_tend_b_4th_g<<<gridGPU, blockGPU>>>(
                fields.st.at("b")->fld_g,
                fields.mp.at("u")->fld_g, fields.mp.at("w")->fld_g,
                gd.utrans, bs.n2, sinalpha, cosalpha,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();
        }
        else
        {
            calc_buoyancy_tend_4th_g<<<gridGPU, blockGPU>>>(
            fields.mt.at("w")->fld_g, fields.sp.at("b")->fld_g,
            gd.istart, gd.jstart, gd.kstart+1,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
            cuda_check_error();
        }

        if (swbaroclinic)
        {
            calc_baroclinic_4th_g<<<gridGPU, blockGPU>>>(
                fields.st.at("b")->fld_g, fields.mp.at("v")->fld_g,
                dbdy_ls,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();
        }
    }
    cudaDeviceSynchronize();
    stats.calc_tend(*fields.mt.at("w"), tend_name);

}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_buoy<TF>::get_thermo_field_g(
        Field3d<TF>& fld, const std::string& name, const bool cyclic)
{
    auto& gd = grid.get_grid_data();

    int blocksize = min(256, 16 * ((gd.ncells / 16) + (gd.ncells % 16 > 0)));

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    if (name == "b")
    {
        calc_buoyancy_g<<<gd.ncells, blocksize>>>(
            fld.fld_g, fields.sp.at("b")->fld_g);
        cuda_check_error();
    }
    else if (name == "N2")
    {
        calc_N2_g<<<gridGPU, blockGPU>>>(
            fld.fld_g, fields.sp.at("b")->fld_g, bs.n2, gd.dzi_g,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
        cuda_check_error();
    }
}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_buoy<TF>::get_buoyancy_fluxbot_g(Field3d<TF>& b)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU (gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);

    calc_buoyancy_flux_bot_g<<<gridGPU, blockGPU>>>(
        b.flux_bot_g, fields.sp.at("b")->flux_bot_g,
        gd.kstart, gd.icells, gd.jcells,
        gd.icells, gd.ijcells);
    cuda_check_error();
}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_buoy<TF>::get_buoyancy_surf_g(Field3d<TF>& b)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU (gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);

    calc_buoyancy_bot_g<<<gridGPU, blockGPU>>>(
        b.fld_g, b.flux_bot_g,
        fields.sp.at("b")->fld_g, fields.sp.at("b")->fld_bot_g,
        gd.kstart, gd.icells, gd.jcells,
        gd.icells, gd.ijcells);
    cuda_check_error();

    calc_buoyancy_flux_bot_g<<<gridGPU, blockGPU>>>(
        b.flux_bot_g, fields.sp.at("b")->flux_bot_g,
        gd.kstart, gd.icells, gd.jcells,
        gd.icells, gd.ijcells);
    cuda_check_error();
}
#endif


#ifdef FLOAT_SINGLE
template class Thermo_buoy<float>;
#else
template class Thermo_buoy<double>;
#endif
