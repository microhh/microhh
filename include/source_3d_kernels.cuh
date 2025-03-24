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

#ifndef SOURCE_GAUSSIAN_KERNELS_G_H
#define SOURCE_GAUSSIAN_KERNELS_G_H

namespace Source_3d_kernels_g
{
    template<typename TF> __global__
    void add_source_tend_g(
            TF* const __restrict__ st_out,
            const TF* const __restrict__ st_in,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k  = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk_in = i + j*jstride + (k-kstart)*kstride;
            const int ijk_out = i + j*jstride + k*kstride;

            st_out[ijk_out] += st_in[ijk_in];
        }
    }


    template<typename TF> __global__
    void add_source_tend_heat_g(
            TF* const __restrict__ st_out,
            const TF* const __restrict__ Te,    // Absolute emission temperature (K).
            const TF* const __restrict__ Qe,    // Volume flux (m3 s-1).
            const TF* const __restrict__ T,     // Absolute temperature LES.
            const TF* const __restrict__ dz,
            const TF dx,
            const TF dy,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k  = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk_in = i + j*jstride + (k-kstart)*kstride;
            const int ijk = i + j*jstride + k*kstride;

            const TF Vi = TF(1) / (dx * dy * dz[k]);
            st_out[ijk] += Qe[ijk_in] * (Te[ijk_in] - T[ijk]) * Vi;
        }
    }


    template<typename TF> __global__
    void interpolate_emission_g(
            TF* const __restrict__ emission_out,
            const TF* const __restrict__ emission_prev,
            const TF* const __restrict__ emission_next,
            const TF fac0,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k  = blockIdx.z;

        const TF fac1 = TF(1) - fac0;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jstride + k*kstride;
            emission_out[ijk] = fac0 * emission_prev[ijk] + fac1 * emission_next[ijk];
        }
    }
}
#endif