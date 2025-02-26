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

#ifndef SOURCE_3D_KERNELS_H
#define SOURCE_3D_KERNELS_H


namespace Source_3d_kernels
{
    template<typename TF>
    void add_source_tend(
            TF* const __restrict__ st_out,
            const TF* const __restrict__ st_in,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        for(int k = kstart; k<kend; ++k)
            for(int j = jstart; j<jend; ++j)
                #pragma ivdep
                for(int i = istart; i<iend; ++i)
                {
                    const int ijk_in = i + j*jstride + (k-kstart)*kstride;
                    const int ijk_out = i + j*jstride + k*kstride;

                    st_out[ijk_out] += st_in[ijk_in];
                }
    }


    template<typename TF>
    void interpolate_emission(
            TF* const __restrict__ emission_out,
            const TF* const __restrict__ emission_prev,
            const TF* const __restrict__ emission_next,
            const TF fac0,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        const TF fac1 = TF(1) - fac0;

        for(int k = kstart; k<kend; ++k)
            for(int j = jstart; j<jend; ++j)
                #pragma ivdep
                for(int i = istart; i<iend; ++i)
                {
                    const int ijk = i + j*jstride + k*kstride;
                    emission_out[ijk] = fac0 * emission_prev[ijk] + fac1 * emission_next[ijk];
                }
    }
}
#endif
