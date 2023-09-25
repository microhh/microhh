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

#ifndef MICROHHC_DIFF_TKE2_KL_KERNELS_CUH
#define MICROHHC_DIFF_TKE2_KL_KERNELS_CUH

#include "cuda_tiling.h"
#include "fast_math.h"
#include "monin_obukhov.h"

namespace diff_tke2
{
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;

    template<typename TF, bool sw_surface_model, bool sw_mason>
    struct evisc_g {
        DEFINE_GRID_KERNEL("diff_tke2::evisc", sw_surface_model ? 1 : 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout gd,
                int i, int j, int k,
                Level level,
                TF* const __restrict__ evisc,
                const TF* const __restrict__ sgstke,
                const TF* const __restrict__ N2,
                const TF* const __restrict__ bgradbot,
                const TF* const __restrict__ z,
                const TF* const __restrict__ z0m,
                const TF* const __restrict__ mlen0,
                const TF cn, const TF cm)
        {
            const int ij  = i + j*gd.jj;
            const int ijk = i + j*gd.jj + k*gd.kk;

            // Variables for the wall damping and length scales
            const TF n_mason = TF(2.);
            TF mlen;
            TF fac;
    
            if (!sw_surface_model)
                asm("trap;");
            else
            {
                mlen = mlen0[k];
    
                if (level.distance_to_start() == 0)
                {
                    if ( bgradbot[ij] > 0 ) // Only if stably stratified, adapt length scale
                        mlen = cn * sqrt(sgstke[ijk]) / sqrt(bgradbot[ij]);
                }
                else
                {
                    if ( N2[ijk] > 0 ) // Only if stably stratified, adapt length scale
                        mlen = cn * sqrt(sgstke[ijk]) / sqrt(N2[ijk]);
                }
    
                fac  = min(mlen0[k], mlen);
    
                if (sw_mason) // Apply Mason's wall correction here
                    fac = pow(TF(1.)/(TF(1.)/pow(fac, n_mason) + TF(1.)/
                            (pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);
    
                // Calculate eddy diffusivity for momentum.
                evisc[ijk] = cm * fac * sqrt(sgstke[ijk]);
            }
        }
    };


    template<typename TF, bool sw_surface_model, bool sw_mason>
    struct evisc_heat_g {
        DEFINE_GRID_KERNEL("diff_tke2::evisc_heat", sw_surface_model ? 1 : 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout gd,
                int i, int j, int k,
                Level level,
                TF* const __restrict__ evisch,
                const TF* __restrict__ evisc,
                const TF* __restrict__ sgstke,
                const TF* __restrict__ N2,
                const TF* __restrict__ bgradbot,
                const TF* __restrict__ z,
                const TF* __restrict__ z0m,
                const TF* __restrict__ mlen0,
                const TF cn, const TF ch1, const TF ch2)
        {
            const int ij  = i + j*gd.jj;
            const int ijk = i + j*gd.jj + k*gd.kk;

            // Variables for the wall damping and length scales
            const TF n_mason = TF(2.);
            TF mlen;
            TF fac;

            if (!sw_surface_model)
                asm("trap;");
            else
            {
                mlen = mlen0[k];

                if (level.distance_to_start() == 0)
                {
                    if ( bgradbot[ij] > 0 ) // Only if stably stratified, adapt length scale
                        mlen = cn * sqrt(sgstke[ijk]) / sqrt(bgradbot[ij]);
                }
                else
                {
                    if ( N2[ijk] > 0 ) // Only if stably stratified, adapt length scale
                        mlen = cn * sqrt(sgstke[ijk]) / sqrt(N2[ijk]);
                }

                fac  = min(mlen0[k], mlen);

                if (sw_mason) // Apply Mason's wall correction here
                    fac = pow(TF(1.)/(TF(1.)/pow(fac, n_mason) + TF(1.)/
                                (pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);

                // Calculate eddy diffusivity for momentum.
                evisch[ijk] = (ch1 + ch2 * fac / mlen0[k]) * evisc[ijk];
            }
        }
    };
}
#endif
