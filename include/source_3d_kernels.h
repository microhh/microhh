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

#include "constants.h"

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
    void add_source_tend_heat(
            TF* const __restrict__ th_tend,
            const TF* const __restrict__ Te,      // Absolute emission temperature (K).
            const TF* const __restrict__ Me,      // Emission mass flux (kg s-1).
            const TF* const __restrict__ T,       // Absolute temperature LES (K).
            const TF* const __restrict__ rhoref,  // Base state density (kg m-3).
            const TF* const __restrict__ dz,
            const TF* const __restrict__ exner,
            const TF dx,
            const TF dy,
            const TF subdti,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        for(int k = kstart; k<kend; ++k)
        {
            const TF mi = TF(1) / (rhoref[k] * dx * dy * dz[k]);
            const TF exneri = TF(1) / exner[k];

            for(int j = jstart; j<jend; ++j)
                #pragma ivdep
                for(int i = istart; i<iend; ++i)
                {
                    const int ijk_in = i + j*jstride + (k-kstart)*kstride;
                    const int ijk = i + j*jstride + k*kstride;

                    const TF f = Me[ijk_in] * mi;
                    const TF fac = std::min(TF(1), subdti / std::max(f, TF(Constants::dtiny)));

                    th_tend[ijk] += fac * f * (Te[ijk_in] - T[ijk]) * exneri;
                }
        }
    }


    template<typename TF>
    void add_source_tend_moisture(
            TF* const __restrict__ qt_tend,
            const TF* const __restrict__ qe,      // Specific humidity of emission (kg kg-1).
            const TF* const __restrict__ Me,      // Emission mass flux (kg s-1).
            const TF* const __restrict__ qt,      // Specific humidity LES.
            const TF* const __restrict__ rhoref,  // Base state density (kg m-3).
            const TF* const __restrict__ dz,
            const TF dx,
            const TF dy,
            const TF subdti,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        for(int k = kstart; k<kend; ++k)
        {
            const TF mi = TF(1) / (rhoref[k] * dx * dy * dz[k]);

            for(int j = jstart; j<jend; ++j)
                #pragma ivdep
                for(int i = istart; i<iend; ++i)
                {
                    const int ijk_in = i + j*jstride + (k-kstart)*kstride;
                    const int ijk = i + j*jstride + k*kstride;

                    const TF f = Me[ijk_in] * mi;
                    const TF fac = std::min(TF(1), subdti / std::max(f, TF(Constants::dtiny)));

                    qt_tend[ijk] += fac * f * (qe[ijk_in] - qt[ijk]);
                }
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
