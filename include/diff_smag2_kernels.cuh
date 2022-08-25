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

#ifndef MICROHHC_DIFF_SMAG2_KERNELS_CUH
#define MICROHHC_DIFF_SMAG2_KERNELS_CUH

#include "cuda_tiling.h"
#include "fast_math.h"
#include "monin_obukhov.h"

namespace diff_smag2 {
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;

    enum class Surface_model {Enabled, Disabled};


    template<typename TF, Surface_model surface_model>
    struct calc_strain2_g {
        DEFINE_GRID_KERNEL("diff_smag2__calc_strain2", 1)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout gd, int i, int j, int k, Level level,
                TF* __restrict__ strain2,
                const TF* __restrict__ u,
                const TF* __restrict__ v,
                const TF* __restrict__ w,
                const TF* __restrict__ dudz,
                const TF* __restrict__ dvdz,
                const TF* __restrict__ dzi,
                const TF* __restrict__ dzhi,
                const TF dxi, const TF dyi)
        {
            const int ii = 1;
            const int jj = gd.jj;
            const int kk = gd.kk;
            const int ij  = i*ii + j*jj;
            const int ijk = i*ii + j*jj + k*kk;

            if (level.distance_to_start() == 0 && surface_model == Surface_model::Enabled)
            {
                strain2[ijk] = TF(2.)*(
                        // du/dx + du/dx
                        + fm::pow2((u[ijk+ii]-u[ijk])*dxi)

                        // dv/dy + dv/dy
                        + fm::pow2((v[ijk+jj]-v[ijk])*dyi)

                        // dw/dz + dw/dz
                        + fm::pow2((w[ijk+kk]-w[ijk])*dzi[k])

                        // du/dy + dv/dx
                        + TF(0.125)*fm::pow2((u[ijk      ]-u[ijk   -jj])*dyi  + (v[ijk      ]-v[ijk-ii   ])*dxi)
                        + TF(0.125)*fm::pow2((u[ijk+ii   ]-u[ijk+ii-jj])*dyi  + (v[ijk+ii   ]-v[ijk      ])*dxi)
                        + TF(0.125)*fm::pow2((u[ijk   +jj]-u[ijk      ])*dyi  + (v[ijk   +jj]-v[ijk-ii+jj])*dxi)
                        + TF(0.125)*fm::pow2((u[ijk+ii+jj]-u[ijk+ii   ])*dyi  + (v[ijk+ii+jj]-v[ijk   +jj])*dxi)

                        // du/dz + dw/dx
                        + TF(0.5)*fm::pow2(dudz[ij])
                        + TF(0.125)*fm::pow2((w[ijk      ]-w[ijk-ii   ])*dxi)
                        + TF(0.125)*fm::pow2((w[ijk+ii   ]-w[ijk      ])*dxi)
                        + TF(0.125)*fm::pow2((w[ijk   +kk]-w[ijk-ii+kk])*dxi)
                        + TF(0.125)*fm::pow2((w[ijk+ii+kk]-w[ijk   +kk])*dxi)

                        // dv/dz + dw/dy
                        + TF(0.5)*fm::pow2(dvdz[ij])
                        + TF(0.125)*fm::pow2((w[ijk      ]-w[ijk-jj   ])*dyi)
                        + TF(0.125)*fm::pow2((w[ijk+jj   ]-w[ijk      ])*dyi)
                        + TF(0.125)*fm::pow2((w[ijk   +kk]-w[ijk-jj+kk])*dyi)
                        + TF(0.125)*fm::pow2((w[ijk+jj+kk]-w[ijk   +kk])*dyi) );

                // add a small number to avoid zero divisions
                strain2[ijk] += Constants::dsmall;
            }
            else
            {
                strain2[ijk] = TF(2.)*(
                        // du/dx + du/dx
                        + fm::pow2((u[ijk+ii]-u[ijk])*dxi)
                        // dv/dy + dv/dy
                        + fm::pow2((v[ijk+jj]-v[ijk])*dyi)
                        // dw/dz + dw/dz
                        + fm::pow2((w[ijk+kk]-w[ijk])*dzi[k])
                        // du/dy + dv/dx
                        + TF(0.125)*fm::pow2((u[ijk      ]-u[ijk   -jj])*dyi  + (v[ijk      ]-v[ijk-ii   ])*dxi)
                        + TF(0.125)*fm::pow2((u[ijk+ii   ]-u[ijk+ii-jj])*dyi  + (v[ijk+ii   ]-v[ijk      ])*dxi)
                        + TF(0.125)*fm::pow2((u[ijk   +jj]-u[ijk      ])*dyi  + (v[ijk   +jj]-v[ijk-ii+jj])*dxi)
                        + TF(0.125)*fm::pow2((u[ijk+ii+jj]-u[ijk+ii   ])*dyi  + (v[ijk+ii+jj]-v[ijk   +jj])*dxi)
                        // du/dz + dw/dx
                        + TF(0.125)*fm::pow2((u[ijk      ]-u[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-ii   ])*dxi)
                        + TF(0.125)*fm::pow2((u[ijk+ii   ]-u[ijk+ii-kk])*dzhi[k  ] + (w[ijk+ii   ]-w[ijk      ])*dxi)
                        + TF(0.125)*fm::pow2((u[ijk   +kk]-u[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-ii+kk])*dxi)
                        + TF(0.125)*fm::pow2((u[ijk+ii+kk]-u[ijk+ii   ])*dzhi[k+1] + (w[ijk+ii+kk]-w[ijk   +kk])*dxi)
                        // dv/dz + dw/dy
                        + TF(0.125)*fm::pow2((v[ijk      ]-v[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-jj   ])*dyi)
                        + TF(0.125)*fm::pow2((v[ijk+jj   ]-v[ijk+jj-kk])*dzhi[k  ] + (w[ijk+jj   ]-w[ijk      ])*dyi)
                        + TF(0.125)*fm::pow2((v[ijk   +kk]-v[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-jj+kk])*dyi)
                        + TF(0.125)*fm::pow2((v[ijk+jj+kk]-v[ijk+jj   ])*dzhi[k+1] + (w[ijk+jj+kk]-w[ijk   +kk])*dyi) );

                // add a small number to avoid zero divisions
                strain2[ijk] += Constants::dsmall;
            }
        }
    };

    template<typename TF, Surface_model surface_model>
    struct evisc_g {
        DEFINE_GRID_KERNEL("diff_smag2__evisc", surface_model == Surface_model::Enabled ? 1 : 0)

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
            const TF n_mason = TF(2);

            const int kstart = gd.kstart;
            const int jj = gd.jj;
            const int kk = gd.kk;
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + k*kk;

            if (level.distance_to_start() == 0 && surface_model == Surface_model::Enabled)
            {
                // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
                TF RitPrratio = bgradbot[ij] / evisc[ijk] * tPri;
                RitPrratio = fmin(RitPrratio, TF(1.-Constants::dsmall));

                const TF mlen = pow(TF(1.)/(TF(1.)/mlen0[k] + TF(1.)/(pow(Constants::kappa<TF>*(z[kstart]+z0m[ij]), n_mason))), TF(1.)/n_mason);
                evisc[ijk] = fm::pow2(mlen) * sqrt(evisc[ijk] * (TF(1.)-RitPrratio));
            }
            else if (surface_model == Surface_model::Enabled)
            {
                // Add the buoyancy production to the TKE
                TF RitPrratio = N2[ijk] / evisc[ijk] * tPri;
                RitPrratio = fmin(RitPrratio, TF(1.-Constants::dsmall));

                // Mason mixing length
                const TF mlen = pow(TF(1.)/(TF(1.)/mlen0[k] + TF(1.)/(pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);
                evisc[ijk] = fm::pow2(mlen) * sqrt(evisc[ijk] * (TF(1.)-RitPrratio));
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

    template<typename TF> __global__
    void calc_ghostcells_evisc(
            TF* __restrict__ evisc,
            const int icells, const int jcells,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int kb = kstart;
            const int kt = kend-1;

            const int ijkb = i + j*jj + kb*kk;
            const int ijkt = i + j*jj + kt*kk;

            evisc[ijkb-kk] = evisc[ijkb];
            evisc[ijkt+kk] = evisc[ijkt];
        }
    }

    template<typename TF, Surface_model surface_model>
    struct diff_uvw_g {
        DEFINE_GRID_KERNEL("diff_smag2__diff_uvw", surface_model == Surface_model::Enabled ? 1 : 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout gd, int i, int j, int k, Level level,
                TF* __restrict__ ut, TF* __restrict__ vt, TF* __restrict__ wt,
                const TF* __restrict__ evisc,
                const TF* __restrict__ u, const TF* __restrict__ v, const TF* __restrict__ w,
                const TF* __restrict__ fluxbotu, const TF* __restrict__ fluxtopu,
                const TF* __restrict__ fluxbotv, const TF* __restrict__ fluxtopv,
                const TF* __restrict__ dzi, const TF* __restrict__ dzhi, const TF dxi, const TF dyi,
                const TF* __restrict__ rhoref, const TF* __restrict__ rhorefh,
                const TF* __restrict__ rhorefi, const TF* __restrict__ rhorefhi,
                const TF visc)
        {
            const int ii = gd.ii;
            const int jj = gd.jj;
            const int kk = gd.kk;
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + k*kk;

            // U
            const TF evisceu = evisc[ijk   ] + visc;
            const TF eviscwu = evisc[ijk-ii] + visc;
            const TF eviscnu = TF(0.25)*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]) + visc;
            const TF eviscsu = TF(0.25)*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]) + visc;
            const TF evisctu = TF(0.25)*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]) + visc;
            const TF eviscbu = TF(0.25)*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]) + visc;

            // V
            const TF eviscev = TF(0.25)*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]) + visc;
            const TF eviscwv = TF(0.25)*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]) + visc;
            const TF eviscnv = evisc[ijk   ] + visc;
            const TF eviscsv = evisc[ijk-jj] + visc;
            const TF evisctv = TF(0.25)*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]) + visc;
            const TF eviscbv = TF(0.25)*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]) + visc;

            // W
            const TF eviscew = TF(0.25)*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+ii-kk] + evisc[ijk+ii]) + visc;
            const TF eviscww = TF(0.25)*(evisc[ijk-ii-kk] + evisc[ijk-ii] + evisc[ijk   -kk] + evisc[ijk   ]) + visc;
            const TF eviscnw = TF(0.25)*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+jj-kk] + evisc[ijk+jj]) + visc;
            const TF eviscsw = TF(0.25)*(evisc[ijk-jj-kk] + evisc[ijk-jj] + evisc[ijk   -kk] + evisc[ijk   ]) + visc;
            const TF evisctw = evisc[ijk   ] + visc;
            const TF eviscbw = evisc[ijk-kk] + visc;

            if (level.distance_to_start() == 0 && surface_model == Surface_model::Enabled)
            {
                ut[ijk] +=
                        // du/dx + du/dx
                        + (  evisceu*(u[ijk+ii]-u[ijk   ])*dxi
                             - eviscwu*(u[ijk   ]-u[ijk-ii])*dxi ) * TF(2.)*dxi
                        // du/dy + dv/dx
                        + (  eviscnu*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                             - eviscsu*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                        // du/dz + dw/dx
                        + (  rhorefh[k+1] * evisctu*((u[ijk+kk]-u[ijk   ])* dzhi[k+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
                             + rhorefh[k  ] * fluxbotu[ij] ) * rhorefi[k] * dzi[k];

                vt[ijk] +=
                        // dv/dx + du/dy
                        + (  eviscev*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                             - eviscwv*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                        // dv/dy + dv/dy
                        + (  eviscnv*(v[ijk+jj]-v[ijk   ])*dyi
                             - eviscsv*(v[ijk   ]-v[ijk-jj])*dyi ) * TF(2.)*dyi
                        // dv/dz + dw/dy
                        + (  rhorefh[k+1] * evisctv*((v[ijk+kk]-v[ijk   ])*dzhi[k+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
                             + rhorefh[k  ] * fluxbotv[ij] ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_end() == 0 && surface_model == Surface_model::Enabled)
            {
                ut[ijk] +=
                        // du/dx + du/dx
                        + (  evisceu*(u[ijk+ii]-u[ijk   ])*dxi
                             - eviscwu*(u[ijk   ]-u[ijk-ii])*dxi ) * TF(2.)*dxi
                        // du/dy + dv/dx
                        + (  eviscnu*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                             - eviscsu*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                        // du/dz + dw/dx
                        + (- rhorefh[k+1] * fluxtopu[ij]
                           - rhorefh[k] * eviscbu*((u[ijk   ]-u[ijk-kk])* dzhi[k] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) * rhorefi[k] * dzi[k];

                vt[ijk] +=
                        // dv/dx + du/dy
                        + (  eviscev*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                             - eviscwv*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                        // dv/dy + dv/dy
                        + (  eviscnv*(v[ijk+jj]-v[ijk   ])*dyi
                             - eviscsv*(v[ijk   ]-v[ijk-jj])*dyi ) * TF(2.)*dyi
                        // dv/dz + dw/dy
                        + (- rhorefh[k  ] * fluxtopv[ij]
                           - rhorefh[k-1] * eviscbv*((v[ijk   ]-v[ijk-kk])*dzhi[k-1] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) * rhorefi[k-1] * dzi[k-1];

                wt[ijk] +=
                        // dw/dx + du/dz
                        + (  eviscew*((w[ijk+ii]-w[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-kk])*dzhi[k])
                             - eviscww*((w[ijk   ]-w[ijk-ii])*dxi + (u[ijk   ]-u[ijk+  -kk])*dzhi[k]) ) * dxi
                        // dw/dy + dv/dz
                        + (  eviscnw*((w[ijk+jj]-w[ijk   ])*dyi + (v[ijk+jj]-v[ijk+jj-kk])*dzhi[k])
                             - eviscsw*((w[ijk   ]-w[ijk-jj])*dyi + (v[ijk   ]-v[ijk+  -kk])*dzhi[k]) ) * dyi
                        // dw/dz + dw/dz
                        + (  rhoref[k  ] * evisctw*(w[ijk+kk]-w[ijk   ])*dzi[k  ]
                             - rhoref[k-1] * eviscbw*(w[ijk   ]-w[ijk-kk])*dzi[k-1] ) * rhorefhi[k] * TF(2.) * dzhi[k];
            }
            else
            {
                ut[ijk] +=
                        // du/dx + du/dx
                        + (  evisceu*(u[ijk+ii]-u[ijk   ])*dxi
                             - eviscwu*(u[ijk   ]-u[ijk-ii])*dxi ) * TF(2.)*dxi
                        // du/dy + dv/dx
                        + (  eviscnu*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                             - eviscsu*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                        // du/dz + dw/dx
                        + (  rhorefh[k+1] * evisctu*((u[ijk+kk]-u[ijk   ])* dzhi[k+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
                             - rhorefh[k  ] * eviscbu*((u[ijk   ]-u[ijk-kk])* dzhi[k  ] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) * rhorefi[k] * dzi[k];

                vt[ijk] +=
                        // dv/dx + du/dy
                        + (  eviscev*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                             - eviscwv*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                        // dv/dy + dv/dy
                        + (  eviscnv*(v[ijk+jj]-v[ijk   ])*dyi
                             - eviscsv*(v[ijk   ]-v[ijk-jj])*dyi ) * TF(2.)*dyi
                        // dv/dz + dw/dy
                        + (  rhorefh[k+1] * evisctv*((v[ijk+kk]-v[ijk   ])*dzhi[k+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
                             - rhorefh[k  ] * eviscbv*((v[ijk   ]-v[ijk-kk])*dzhi[k  ] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) * rhorefi[k] * dzi[k];

                wt[ijk] +=
                        // dw/dx + du/dz
                        + (  eviscew*((w[ijk+ii]-w[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-kk])*dzhi[k])
                             - eviscww*((w[ijk   ]-w[ijk-ii])*dxi + (u[ijk   ]-u[ijk+  -kk])*dzhi[k]) ) * dxi
                        // dw/dy + dv/dz
                        + (  eviscnw*((w[ijk+jj]-w[ijk   ])*dyi + (v[ijk+jj]-v[ijk+jj-kk])*dzhi[k])
                             - eviscsw*((w[ijk   ]-w[ijk-jj])*dyi + (v[ijk   ]-v[ijk+  -kk])*dzhi[k]) ) * dyi
                        // dw/dz + dw/dz
                        + (  rhoref[k  ] * evisctw*(w[ijk+kk]-w[ijk   ])*dzi[k  ]
                             - rhoref[k-1] * eviscbw*(w[ijk   ]-w[ijk-kk])*dzi[k-1] ) * rhorefhi[k] * TF(2.) * dzhi[k];
            }
        }
    };

    template<typename TF, Surface_model surface_model>
    struct diff_c_g {
        DEFINE_GRID_KERNEL("diff_smag2__diff_c", surface_model == Surface_model::Enabled ? 1 : 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout gd, int i, int j, int k, Level level,
                TF* __restrict__ at,
                const TF* __restrict__ a, const TF* __restrict__ evisc,
                const TF* __restrict__ fluxbot, const TF* __restrict__ fluxtop,
                const TF* __restrict__ dzi, const TF* __restrict__ dzhi, const TF dxidxi, const TF dyidyi,
                const TF* __restrict__ rhorefi, const TF* __restrict__ rhorefh,
                const TF tPri, const TF visc)
        {
            const int ii = gd.ii;
            const int jj = gd.jj;
            const int kk = gd.kk;
            const int ij  = i*ii + j*jj;
            const int ijk = i*ii + j*jj + k*kk;

            if (level.distance_to_start() == 0 && surface_model == Surface_model::Enabled)
            {
                const TF evisce = TF(0.5)*(evisc[ijk   ]+evisc[ijk+ii])*tPri + visc;
                const TF eviscw = TF(0.5)*(evisc[ijk-ii]+evisc[ijk   ])*tPri + visc;
                const TF eviscn = TF(0.5)*(evisc[ijk   ]+evisc[ijk+jj])*tPri + visc;
                const TF eviscs = TF(0.5)*(evisc[ijk-jj]+evisc[ijk   ])*tPri + visc;
                const TF evisct = TF(0.5)*(evisc[ijk   ]+evisc[ijk+kk])*tPri + visc;

                at[ijk] +=
                        + (  evisce*(a[ijk+ii]-a[ijk   ])
                             - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                        + (  eviscn*(a[ijk+jj]-a[ijk   ])
                             - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                        + (  rhorefh[k+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[k+1]
                             + rhorefh[k  ] * fluxbot[ij] ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_end() == 0 && surface_model == Surface_model::Enabled)
            {
                const TF evisce = TF(0.5)*(evisc[ijk   ]+evisc[ijk+ii])*tPri + visc;
                const TF eviscw = TF(0.5)*(evisc[ijk-ii]+evisc[ijk   ])*tPri + visc;
                const TF eviscn = TF(0.5)*(evisc[ijk   ]+evisc[ijk+jj])*tPri + visc;
                const TF eviscs = TF(0.5)*(evisc[ijk-jj]+evisc[ijk   ])*tPri + visc;
                const TF eviscb = TF(0.5)*(evisc[ijk-kk]+evisc[ijk   ])*tPri + visc;

                at[ijk] +=
                        + (  evisce*(a[ijk+ii]-a[ijk   ])
                             - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                        + (  eviscn*(a[ijk+jj]-a[ijk   ])
                             - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                        + (- rhorefh[k  ] * fluxtop[ij]
                           - rhorefh[k-1] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[k-1] ) * rhorefi[k-1] * dzi[k-1];
            }
            else
            {
                const TF evisce = TF(0.5)*(evisc[ijk   ]+evisc[ijk+ii])*tPri + visc;
                const TF eviscw = TF(0.5)*(evisc[ijk-ii]+evisc[ijk   ])*tPri + visc;
                const TF eviscn = TF(0.5)*(evisc[ijk   ]+evisc[ijk+jj])*tPri + visc;
                const TF eviscs = TF(0.5)*(evisc[ijk-jj]+evisc[ijk   ])*tPri + visc;
                const TF evisct = TF(0.5)*(evisc[ijk   ]+evisc[ijk+kk])*tPri + visc;
                const TF eviscb = TF(0.5)*(evisc[ijk-kk]+evisc[ijk   ])*tPri + visc;

                at[ijk] +=
                        + (  evisce*(a[ijk+ii]-a[ijk   ])
                             - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                        + (  eviscn*(a[ijk+jj]-a[ijk   ])
                             - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                        + (  rhorefh[k+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[k+1]
                             - rhorefh[k  ] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[k]  ) * rhorefi[k] * dzi[k];
            }
        }
    };

    template<typename TF> __global__
    void calc_dnmul_g(TF* __restrict__ dnmul, TF* __restrict__ evisc,
                      const TF* __restrict__ dzi, TF tPrfac_i, const TF dxidxi, const TF dyidyi,
                      const int istart, const int jstart, const int kstart,
                      const int iend,   const int jend,   const int kend,
                      const int jj,     const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            dnmul[ijk] = fabs(evisc[ijk]*tPrfac_i*(dxidxi + dyidyi + dzi[k]*dzi[k]));
        }
    }
}

#endif //MICROHHC_DIFF_SMAG2_KERNELS_CUH
