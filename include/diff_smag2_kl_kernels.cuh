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

#ifndef MICROHHC_DIFF_SMAG2_KL_KERNELS_CUH
#define MICROHHC_DIFF_SMAG2_KL_KERNELS_CUH

#include "cuda_tiling.h"
#include "fast_math.h"
#include "monin_obukhov.h"

namespace diff_smag2 {
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;

    template<typename TF, bool surface_model_enabled>
    struct evisc_g {
        DEFINE_GRID_KERNEL("diff_smag2::evisc", surface_model_enabled ? 1 : 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
            Grid_layout gd, int i, int j, int k, Level level,
            TF* __restrict__ evisc,
            const TF* __restrict__ N2,
            const TF* __restrict__ bgradbot,
            const TF* __restrict__ mlen0,
            const TF* __restrict__ z0m,
            const TF* __restrict__ z,
            const TF tPri)
        {
//            const TF n_mason = TF(2);

            const int jj = gd.jj;
            const int kk = gd.kk;
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + k*kk;

            if (surface_model_enabled)
            {
                TF RitPrratio;

                if (level.distance_to_start() == 0)
                {
                    // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
                    RitPrratio = bgradbot[ij] / evisc[ijk] * tPri;
                }
                else
                {
                    // Add the buoyancy production to the TKE
                    RitPrratio = N2[ijk] / evisc[ijk] * tPri;
                }

                RitPrratio = fmin(RitPrratio, TF(1.-Constants::dsmall));

                // Mason mixing length
                const TF m = mlen0[k];
                const TF r = fm::pow2(Constants::kappa<TF>*(z[k] + z0m[ij]));
                const TF mlen_squared = (m * r) / (m + r); // == 1/(1/r + 1/m)
                evisc[ijk] = mlen_squared * sqrt(evisc[ijk] * (TF(1.) - RitPrratio));

            }
            else
            {
                // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
                TF RitPrratio = N2[ijk] / evisc[ijk] * tPri;
                RitPrratio = fmin(RitPrratio, TF(1.-Constants::dsmall));
                evisc[ijk] = fm::pow2(mlen0[k]) * sqrt(evisc[ijk] * (TF(1.)-RitPrratio));
            }
        }
    };

    template<typename TF> __global__
    void evisc_neutral_g(
            TF* __restrict__ evisc,
            const TF* __restrict__ z0m,
            const TF* __restrict__ z,
            const TF* __restrict__ mlen0,
            const int istart, const int jstart, const int kstart,
            const int iend, const int jend, const int kend,
            const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const TF n_mason = TF(2);

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            const int ij = i + j*jj;

            const TF mlen = pow(TF(1.)/(TF(1.)/mlen0[k] + TF(1.)/(pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);
            evisc[ijk] = fm::pow2(mlen) * sqrt(evisc[ijk]);
        }
    }

    template<typename TF> __global__
    void evisc_neutral_vandriest_g(
            TF* __restrict__ evisc,
            const TF* __restrict__ u, const TF* __restrict__ v,
            const TF* __restrict__ mlen_smag,
            const TF* __restrict__ z, const TF* __restrict__ dzhi,
            const TF zsize, const TF visc,
            const int istart, const int jstart, const int kstart,
            const int iend, const int jend, const int kend,
            const int jj, const int kk)

    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const TF A_vandriest = TF(26.);

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            const int ijk_bot = i + j*jj + kstart*kk;
            const int ijk_top = i + j*jj + kend*kk;

            const TF u_tau_bot = pow(
                    fm::pow2( visc*(u[ijk_bot] - u[ijk_bot-kk] )*dzhi[kstart] )
                    + fm::pow2( visc*(v[ijk_bot] - v[ijk_bot-kk] )*dzhi[kstart] ), TF(0.25) );
            const TF u_tau_top = pow(
                    fm::pow2( visc*(u[ijk_top] - u[ijk_top-kk] )*dzhi[kend] )
                    + fm::pow2( visc*(v[ijk_top] - v[ijk_top-kk] )*dzhi[kend] ), TF(0.25) );

            const TF fac_bot = TF(1.) - exp( -(       z[k] *u_tau_bot) / (A_vandriest*visc) );
            const TF fac_top = TF(1.) - exp( -((zsize-z[k])*u_tau_top) / (A_vandriest*visc) );

            const TF fac = min(fac_bot, fac_top);

            evisc[ijk] = fm::pow2(fac * mlen_smag[k]) * sqrt(evisc[ijk]);
        }
    }
}
#endif //MICROHHC_DIFF_SMAG2_KL_KERNELS_CUH
