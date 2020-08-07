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
#include "thermo_moist.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "master.h"
#include "tools.h"
#include "column.h"
#include "stats.h"
#include "thermo_moist_functions.h"
#include <iostream>

namespace
{
    using namespace Constants;
    using namespace Finite_difference::O2;
    using namespace Thermo_moist_functions;

    template<typename TF>
    inline __device__ Struct_sat_adjust<TF> sat_adjust_g(
            const TF thl, const TF qt, const TF p, const TF exn)
    {
        using Fast_math::pow2;

        int niter = 0;
        const int nitermax = 10;

        TF tnr_old = TF(1.e9);

        const TF tl = thl * exn;
        TF qs = qsat_liq(p, tl);

        Struct_sat_adjust<TF> ans =
        {
            TF(0.), // ql
            TF(0.), // qi
            tl, // t
            qs, // qs
        };

        // Calculate if q-qs(Tl) <= 0. If so, return 0. Else continue with saturation adjustment.
        if (qt-ans.qs <= TF(0.))
            return ans;

        /* Saturation adjustment solver.
         * Root finding function is f(T) = T - tnr - Lv/cp*qt + alpha_w * Lv/cp*qs(T) + alpha_i*Ls/cp*qs(T)
         * dq_sat/dT derivatives can be rewritten using Claussius-Clapeyron (desat/dT = L{v,s}*esat / (Rv*T^2)).
         */

        TF tnr = tl;

        // Warm adjustment.
        if (tl >= T0<TF>)
        {
            while (fabs(tnr-tnr_old)/tnr_old > TF(1.e-5) && niter < nitermax)
            {
                tnr_old = tnr;
                qs = qsat_liq(p, tnr);
                const TF f =
                    tnr - tl - Lv<TF>/cp<TF>*(qt - qs);

                const TF f_prime = TF(1.) + Lv<TF>/cp<TF>*dqsatdT_liq(p, tnr);

                tnr -= f / f_prime;

                niter += 1;
            }

            qs = qsat_liq(p, tnr);
            ans.ql = fmax(TF(0.), qt - qs);
            ans.t  = tnr;
            ans.qs = qs;
        }
        // Cold adjustment.
        else
        {
            while (fabs(tnr-tnr_old)/tnr_old > TF(1.e-5) && niter < nitermax)
            {
                tnr_old = tnr;
                qs = qsat(p, tnr);
                const TF alpha_w = water_fraction(tnr);
                const TF alpha_i = TF(1.) - alpha_w;
                const TF dalphadT = (alpha_w > TF(0.) && alpha_w < TF(1.)) ? TF(0.025) : TF(0.);
                const TF dqsatdT_w = dqsatdT_liq(p, tnr);
                const TF dqsatdT_i = dqsatdT_ice(p, tnr);

                const TF f =
                    tnr - tl - alpha_w*Lv<TF>/cp<TF>*qt - alpha_i*Ls<TF>/cp<TF>*qt
                             + alpha_w*Lv<TF>/cp<TF>*qs + alpha_i*Ls<TF>/cp<TF>*qs;

                const TF f_prime = TF(1.)
                    - dalphadT*Lv<TF>/cp<TF>*qt + dalphadT*Ls<TF>/cp<TF>*qt
                    + dalphadT*Lv<TF>/cp<TF>*qs - dalphadT*Ls<TF>/cp<TF>*qs
                    + alpha_w*Lv<TF>/cp<TF>*dqsatdT_w
                    + alpha_i*Ls<TF>/cp<TF>*dqsatdT_i;

                tnr -= f / f_prime;

                niter += 1;
            }

            const TF alpha_w = water_fraction(tnr);
            const TF alpha_i = TF(1.) - alpha_w;

            qs = qsat(p, tnr);
            const TF ql_qi = fmax(TF(0.), qt - qs);

            ans.ql = alpha_w*ql_qi;
            ans.qi = alpha_i*ql_qi;
            ans.t  = tnr;
            ans.qs = qs;
        }

        // Raise exception if nitermax is reached.
        if (niter == nitermax)
        {
            printf("ERROR: saturation adjustment did not converge!\n");
            asm("trap;");
        }

        return ans;
    }

    template<typename TF> __global__
    void calc_buoyancy_tend_2nd_g(TF* __restrict__ wt, TF* __restrict__ th, TF* __restrict__ qt,
                                  TF* __restrict__ thvrefh, TF* __restrict__ exnh, TF* __restrict__ ph,
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

            // Half level temperature and moisture content
            const TF thh = static_cast<TF>(0.5) * (th[ijk-kk] + th[ijk]); // Half level liq. water pot. temp.
            const TF qth = static_cast<TF>(0.5) * (qt[ijk-kk] + qt[ijk]); // Half level specific hum.

            Struct_sat_adjust<TF> ssa = sat_adjust_g(thh, qth, ph[k], exnh[k]);

            // Calculate tendency.
            if (ssa.ql > 0)
                wt[ijk] += buoyancy(exnh[k], thh, qth, ssa.ql, ssa.qi, thvrefh[k]);
            else
                wt[ijk] += buoyancy_no_ql(thh, qth, thvrefh[k]);
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_g(TF* __restrict__ b,  TF* __restrict__ th,
                         TF* __restrict__ qt, TF* __restrict__ thvref,
                         TF* __restrict__ p,  TF* __restrict__ exn,
                         int istart, int jstart, int kstart,
                         int iend,   int jend,   int kcells,
                         int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z;

        if (i < iend && j < jend && k < kstart)
        {
            const int ijk   = i + j*jj + k*kk;
            b[ijk] = buoyancy_no_ql(th[ijk], qt[ijk], thvref[k]);
        }
        else if (i < iend && j < jend && k < kcells)
        {
            const int ijk = i + j*jj + k*kk;

            Struct_sat_adjust<TF> ssa = sat_adjust_g(th[ijk], qt[ijk], p[k], exn[k]);

            if (ssa.ql > 0)
                b[ijk] = buoyancy(exn[k], th[ijk], qt[ijk], ssa.ql, ssa.qi, thvref[k]);
            else
                b[ijk] = buoyancy_no_ql(th[ijk], qt[ijk], thvref[k]);
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_h_g(TF* __restrict__ bh,  TF* __restrict__ th,
                         TF* __restrict__ qt, TF* __restrict__ thvrefh,
                         TF* __restrict__ ph,  TF* __restrict__ exnh,
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
            const int kk  = i + j*jj;

            // Half level temperature and moisture content
            const TF thh = static_cast<TF>(0.5) * (th[ijk-kk] + th[ijk]); // Half level liq. water pot. temp.
            const TF qth = static_cast<TF>(0.5) * (qt[ijk-kk] + qt[ijk]); // Half level specific hum.

            Struct_sat_adjust<TF> ssa = sat_adjust_g(thh, qth, ph[k], exnh[k]);

            // Calculate tendency
            if (ssa.ql > 0)
                bh[ijk] += buoyancy(exnh[k], thh, qth, ssa.ql, ssa.qi, thvrefh[k]);
            else
                bh[ijk] += buoyancy_no_ql(thh, qth, thvrefh[k]);
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_bot_g(TF* __restrict__ b,      TF* __restrict__ bbot,
                             TF* __restrict__ th,     TF* __restrict__ thbot,
                             TF* __restrict__ qt,     TF* __restrict__ qtbot,
                             TF* __restrict__ thvref, TF* __restrict__ thvrefh,
                             int kstart, int icells, int jcells,
                             int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            bbot[ij ] = buoyancy_no_ql(thbot[ij], qtbot[ij], thvrefh[kstart]);
            b   [ijk] = buoyancy_no_ql(th[ijk],   qt[ijk],   thvref[kstart]);
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_flux_bot_g(TF* __restrict__ bfluxbot,
                                  TF* __restrict__ th, TF* __restrict__ thfluxbot,
                                  TF* __restrict__ qt, TF* __restrict__ qtfluxbot,
                                  TF* __restrict__ thvrefh,
                                  int kstart, int icells, int jcells,
                                  int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        // Calculate the surface buoyancy flux using the first model level temperature and humidity
        // to ensure bitwise identical restarts.
        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            bfluxbot[ij] = buoyancy_flux_no_ql(th[ijk], thfluxbot[ij], qt[ijk], qtfluxbot[ij], thvrefh[kstart]);
        }
    }

    template<typename TF> __global__
    void calc_N2_g(TF* __restrict__ N2, TF* __restrict__ th,
                   TF* __restrict__ thvref, TF* __restrict__ dzi,
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
            N2[ijk] = grav<TF>/thvref[k]*static_cast<TF>(0.5)*(th[ijk+kk] - th[ijk-kk])*dzi[k];
        }
    }

    template<typename TF> __global__
    void calc_liquid_water_g(TF* __restrict__ ql, TF* __restrict__ th, TF* __restrict__ qt,
                             TF* __restrict__ exn, TF* __restrict__ p,
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
            ql[ijk] = sat_adjust_g(th[ijk], qt[ijk], p[k], exn[k]).ql;
        }
    }

    template<typename TF> __global__
    void calc_liquid_water_h_g(TF* __restrict__ qlh, TF* __restrict__ th, TF* __restrict__ qt,
                             TF* __restrict__ exnh, TF* __restrict__ ph,
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
            const int kk  = i + j*jj;

            const TF thh = static_cast<TF>(0.5) * (th[ijk-kk] + th[ijk]); // Half level liq. water pot. temp.
            const TF qth = static_cast<TF>(0.5) * (qt[ijk-kk] + qt[ijk]); // Half level specific hum.
            qlh[ijk] = sat_adjust_g(thh, qth, ph[k], exnh[k]).ql; // Half level liquid water content
        }
    }

    /*
    // BvS: no longer used, base state is calculated at the host
    // CvH: This unused code does not take into account ice
    template <typename TF> __global__
    void calc_base_state_g(TF* __restrict__ pref,     TF* __restrict__ prefh,
                           TF* __restrict__ rho,      TF* __restrict__ rhoh,
                           TF* __restrict__ thv,      TF* __restrict__ thvh,
                           TF* __restrict__ ex,       TF* __restrict__ exh,
                           TF* __restrict__ thlmean,  TF* __restrict__ qtmean,
                           TF* __restrict__ z,        TF* __restrict__ dz,
                           TF* __restrict__ dzh,
                           TF pbot, int kstart, int kend)
    {
        TF ql, si, qti, qli;
        TF rdcp = Rd<TF>/cp<TF>;

        const TF ssurf  = interp2(thlmean[kstart-1], thlmean[kstart]);
        const TF qtsurf = interp2(qtmean[kstart-1],  qtmean[kstart]);

        // Calculate surface (half=kstart) values
        exh[kstart]   = exner(pbot);
        ql            = sat_adjust_g(ssurf,qtsurf,pbot,exh[kstart]).ql;
        thvh[kstart]  = (ssurf + Lv<TF>*ql/(cp<TF>*exh[kstart])) * (1. - (1. - Rv<TF>/Rd<TF>)*qtsurf - Rv<TF>/Rd<TF>*ql);
        prefh[kstart] = pbot;
        rhoh[kstart]  = pbot / (Rd<TF> * exh[kstart] * thvh[kstart]);

        // First full grid level pressure
        pref[kstart] = pow((pow(pbot,rdcp) - grav<TF> * pow(p0<TF>,rdcp) * z[kstart] / (cp<TF> * thvh[kstart])),(1./rdcp));

        for (int k=kstart+1; k<kend+1; k++)
        {
            // 1. Calculate values at full level below zh[k]
            ex[k-1]  = exner(pref[k-1]);
            ql       = sat_adjust_g(thlmean[k-1],qtmean[k-1],pref[k-1],ex[k-1]).ql;
            thv[k-1] = (thlmean[k-1] + Lv<TF>*ql/(cp<TF>*ex[k-1])) * (1. - (1. - Rv<TF>/Rd<TF>)*qtmean[k-1] - Rv<TF>/Rd<TF>*ql);
            rho[k-1] = pref[k-1] / (Rd<TF> * ex[k-1] * thv[k-1]);

            // 2. Calculate half level pressure at zh[k] using values at z[k-1]
            prefh[k] = pow((pow(prefh[k-1],rdcp) - grav<TF> * pow(p0<TF>,rdcp) * dz[k-1] / (cp<TF> * thv[k-1])),(1./rdcp));

            // 3. Interpolate conserved variables to zh[k] and calculate virtual temp and ql
            si     = interp2(thlmean[k-1],thlmean[k]);
            qti    = interp2(qtmean[k-1],qtmean[k]);

            exh[k]   = exner(prefh[k]);
            qli      = sat_adjust_g(si,qti,prefh[k],exh[k]).ql;
            thvh[k]  = (si + Lv<TF>*qli/(cp<TF>*exh[k])) * (1. - (1. - Rv<TF>/Rd<TF>)*qti - Rv<TF>/Rd<TF>*qli);
            rhoh[k]  = prefh[k] / (Rd<TF> * exh[k] * thvh[k]);

            // 4. Calculate full level pressure at z[k]
            pref[k]  = pow((pow(pref[k-1],rdcp) - grav<TF> * pow(p0<TF>,rdcp) * dzh[k] / (cp<TF> * thvh[k])),(1./rdcp));
        }

        // Fill bottom and top full level ghost cells
        pref[kstart-1] = static_cast<TF>(2.)*prefh[kstart] - pref[kstart];
        pref[kend]     = static_cast<TF>(2.)*prefh[kend]   - pref[kend-1];
    }


    // BvS: no longer used, base state is calculated at the host
    template <typename TF> __global__
    void calc_hydrostatic_pressure_g(TF* __restrict__ pref,     TF* __restrict__ prefh,
                                     TF* __restrict__ ex,       TF* __restrict__ exh,
                                     TF* __restrict__ thlmean,  TF* __restrict__ qtmean,
                                     const TF* const __restrict__ z,        const TF* const __restrict__ dz,
                                     const TF* const __restrict__ dzh,
                                     const TF pbot, int kstart, int kend)
    {
        TF ql, si, qti, qli, thvh, thv;
        TF rdcp = Rd<TF>/cp<TF>;

        const TF ssurf  = interp2(thlmean[kstart-1], thlmean[kstart]);
        const TF qtsurf = interp2(qtmean[kstart-1],  qtmean[kstart]);

        // Calculate surface (half=kstart) values
        ql            = sat_adjust_g(ssurf,qtsurf,pbot,exh[kstart]).ql;
        thvh          = (ssurf + Lv<TF>*ql/(cp<TF>*exh[kstart])) * (1. - (1. - Rv<TF>/Rd<TF>)*qtsurf - Rv<TF>/Rd<TF>*ql);
        prefh[kstart] = pbot;

        // First full grid level pressure
        pref[kstart] = pow((pow(pbot,rdcp) - grav<TF> * pow(p0<TF>,rdcp) * z[kstart] / (cp<TF> * thvh)),(1./rdcp));

        for (int k=kstart+1; k<kend+1; k++)
        {
            // 1. Calculate values at full level below zh[k]
            ex[k-1]  = exner(pref[k-1]);
            ql       = sat_adjust_g(thlmean[k-1],qtmean[k-1],pref[k-1],ex[k-1]).ql;
            thv      = (thlmean[k-1] + Lv<TF>*ql/(cp<TF>*ex[k-1])) * (1. - (1. - Rv<TF>/Rd<TF>)*qtmean[k-1] - Rv<TF>/Rd<TF>*ql);

            // 2. Calculate half level pressure at zh[k] using values at z[k-1]
            prefh[k] = pow((pow(prefh[k-1],rdcp) - grav<TF> * pow(p0<TF>,rdcp) * dz[k-1] / (cp<TF> * thv)),(1./rdcp));

            // 3. Interpolate conserved variables to zh[k] and calculate virtual temp and ql
            si     = interp2(thlmean[k-1],thlmean[k]);
            qti    = interp2(qtmean[k-1],qtmean[k]);

            exh[k]   = exner(prefh[k]);
            qli      = sat_adjust_g(si,qti,prefh[k],exh[k]).ql;
            thvh     = (si + Lv<TF>*qli/(cp<TF>*exh[k])) * (1. - (1. - Rv<TF>/Rd<TF>)*qti - Rv<TF>/Rd<TF>*qli);

            // 4. Calculate full level pressure at z[k]
            pref[k]  = pow((pow(pref[k-1],rdcp) - grav<TF> * pow(p0<TF>,rdcp) * dzh[k] / (cp<TF> * thvh)),(1./rdcp));
        }

        // Fill bottom and top full level ghost cells
        pref[kstart-1] = static_cast<TF>(2.)*prefh[kstart] - pref[kstart];
        pref[kend]     = static_cast<TF>(2.)*prefh[kend]   - pref[kend-1];
    }
    */
} // end name    space

template<typename TF>
void Thermo_moist<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();
    const int nmemsize = gd.kcells*sizeof(TF);

    // Allocate fields for Boussinesq and anelastic solver
    cuda_safe_call(cudaMalloc(&bs.thvref_g,  nmemsize));
    cuda_safe_call(cudaMalloc(&bs.thvrefh_g, nmemsize));
    cuda_safe_call(cudaMalloc(&bs.pref_g,    nmemsize));
    cuda_safe_call(cudaMalloc(&bs.prefh_g,   nmemsize));
    cuda_safe_call(cudaMalloc(&bs.exnref_g,  nmemsize));
    cuda_safe_call(cudaMalloc(&bs.exnrefh_g, nmemsize));

    // Copy fields to device
    cuda_safe_call(cudaMemcpy(bs.thvref_g,  bs.thvref.data(),  nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.thvrefh_g, bs.thvrefh.data(), nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.pref_g,    bs.pref.data(),    nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.prefh_g,   bs.prefh.data(),   nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.exnref_g,  bs.exnref.data(),  nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.exnrefh_g, bs.exnrefh.data(), nmemsize, cudaMemcpyHostToDevice));
}

template<typename TF>
void Thermo_moist<TF>::clear_device()
{
    cuda_safe_call(cudaFree(bs.thvref_g ));
    cuda_safe_call(cudaFree(bs.thvrefh_g));
    cuda_safe_call(cudaFree(bs.pref_g   ));
    cuda_safe_call(cudaFree(bs.prefh_g  ));
    cuda_safe_call(cudaFree(bs.exnref_g ));
    cuda_safe_call(cudaFree(bs.exnrefh_g));
    tdep_pbot->clear_device();
}

template<typename TF>
void Thermo_moist<TF>::forward_device()
{
    // Copy fields to device
    auto& gd = grid.get_grid_data();
    const int nmemsize = gd.kcells*sizeof(TF);
    cuda_safe_call(cudaMemcpy(bs.pref_g,    bs.pref.data(),    nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.prefh_g,   bs.prefh.data(),   nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.exnref_g,  bs.exnref.data(),  nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.exnrefh_g, bs.exnrefh.data(), nmemsize, cudaMemcpyHostToDevice));
}

template<typename TF>
void Thermo_moist<TF>::backward_device()
{
    auto& gd = grid.get_grid_data();
    const int nmemsize = gd.kcells*sizeof(TF);
    cudaMemcpy(bs.pref_g,    bs.pref.data(),    nmemsize, cudaMemcpyHostToDevice);
    cudaMemcpy(bs.prefh_g,   bs.prefh.data(),   nmemsize, cudaMemcpyHostToDevice);
    cudaMemcpy(bs.exnref_g,  bs.exnref.data(),  nmemsize, cudaMemcpyHostToDevice);
    cudaMemcpy(bs.exnrefh_g, bs.exnrefh.data(), nmemsize, cudaMemcpyHostToDevice);

    bs_stats = bs;
}

#ifdef USECUDA
template<typename TF>
void Thermo_moist<TF>::exec(const double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    // Re-calculate hydrostatic pressure and exner
    if (bs.swupdatebasestate)
    {
        //calc_hydrostatic_pressure<TF><<<1, 1>>>(bs.pref_g, bs.prefh_g, bs.exnref_g, bs.exnrefh_g,
        //                                        fields.sp.at("thl")->fld_mean_g, fields.sp.at("qt")->fld_mean_g,
        //                                        gd.z_g, gd.dz_g, gd.dzh_g, bs.pbot, gd.kstart, gd.kend);
        //cuda_check_error();

        // BvS: Calculating hydrostatic pressure on GPU is extremely slow. As temporary solution, copy back mean profiles to host,
        //      calculate pressure there and copy back the required profiles.
        cudaMemcpy(fields.sp.at("thl")->fld_mean.data(), fields.sp.at("thl")->fld_mean_g, gd.kcells*sizeof(TF), cudaMemcpyDeviceToHost);
        cudaMemcpy(fields.sp.at("qt")->fld_mean.data(),  fields.sp.at("qt")->fld_mean_g,  gd.kcells*sizeof(TF), cudaMemcpyDeviceToHost);

        auto tmp = fields.get_tmp();

        calc_base_state(
                bs.pref.data(), bs.prefh.data(),
                bs.rhoref.data(), bs.rhorefh.data(), bs.thvref.data(), bs.thvrefh.data(),
                bs.exnref.data(), bs.exnrefh.data(),
                fields.sp.at("thl")->fld_mean.data(), fields.sp.at("qt")->fld_mean.data(),
                bs.pbot, gd.kstart, gd.kend, gd.z.data(), gd.dz.data(), gd.dzh.data());

        fields.release_tmp(tmp);

        // Only half level pressure, exner and virtual potential temperature are needed for calc_buoyancy_tend
        cudaMemcpy(bs.prefh_g,   bs.prefh.data(),   gd.kcells*sizeof(TF), cudaMemcpyHostToDevice);
        cudaMemcpy(bs.exnrefh_g, bs.exnrefh.data(), gd.kcells*sizeof(TF), cudaMemcpyHostToDevice);
        cudaMemcpy(bs.thvrefh_g, bs.thvrefh.data(), gd.kcells*sizeof(TF), cudaMemcpyHostToDevice);
    }

    calc_buoyancy_tend_2nd_g<<<gridGPU, blockGPU>>>(
            fields.mt.at("w")->fld_g, fields.sp.at("thl")->fld_g,
            fields.sp.at("qt")->fld_g, bs.thvrefh_g, bs.exnrefh_g, bs.prefh_g,
            gd.istart, gd.jstart, gd.kstart+1,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
    cuda_check_error();

    cudaDeviceSynchronize();
    stats.calc_tend(*fields.mt.at("w"), tend_name);

}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_moist<TF>::get_thermo_field_g(
        Field3d<TF>& fld, const std::string& name, const bool cyclic )
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    dim3 gridGPU2 (gridi, gridj, gd.kmax);
    dim3 blockGPU2(blocki, blockj, 1);

    // BvS: getthermofield() is called from subgrid-model, before thermo(), so re-calculate the hydrostatic pressure
    if (bs.swupdatebasestate && (name == "b" || name == "ql" || name == "bh" || name == "qlh"))
    {
        //calc_hydrostatic_pressure_g<TF><<<1, 1>>>(bs.pref_g, bs.prefh_g, bs.exnref_g, bs.exnrefh_g,
        //                                          fields.sp.at("thl")->fld_mean_g, fields.sp.at("qt")->fld_mean_g,
        //                                          gd.z_g, gd.dz_g, gd.dzh_g, bs.pbot, gd.kstart, gd.kend);
        //cuda_check_error();

        // BvS: Calculating hydrostatic pressure on GPU is extremely slow. As temporary solution, copy back mean profiles to host,
        //      calculate pressure there and copy back the required profiles.
        cudaMemcpy(fields.sp.at("thl")->fld_mean.data(), fields.sp.at("thl")->fld_mean_g, gd.kcells*sizeof(TF), cudaMemcpyDeviceToHost);
        cudaMemcpy(fields.sp.at("qt")->fld_mean.data(),  fields.sp.at("qt")->fld_mean_g,  gd.kcells*sizeof(TF), cudaMemcpyDeviceToHost);

        auto tmp = fields.get_tmp();

        calc_base_state(bs.pref.data(), bs.prefh.data(),
                        &tmp->fld[0*gd.kcells], &tmp->fld[1*gd.kcells], &tmp->fld[2*gd.kcells], &tmp->fld[3*gd.kcells],
                        bs.exnref.data(), bs.exnrefh.data(), fields.sp.at("thl")->fld_mean.data(), fields.sp.at("qt")->fld_mean.data(),
                        bs.pbot, gd.kstart, gd.kend, gd.z.data(), gd.dz.data(), gd.dzh.data());

        fields.release_tmp(tmp);

        // Only full level pressure and bs.exner needed for calculating buoyancy of ql
        cudaMemcpy(bs.pref_g,   bs.pref.data(),   gd.kcells*sizeof(TF), cudaMemcpyHostToDevice);
        cudaMemcpy(bs.prefh_g,  bs.prefh.data(),  gd.kcells*sizeof(TF), cudaMemcpyHostToDevice);
        cudaMemcpy(bs.exnref_g, bs.exnref.data(), gd.kcells*sizeof(TF), cudaMemcpyHostToDevice);
        cudaMemcpy(bs.exnref_g, bs.exnref.data(), gd.kcells*sizeof(TF), cudaMemcpyHostToDevice);
    }

    if (name == "b")
    {
        calc_buoyancy_g<<<gridGPU, blockGPU>>>(
            fld.fld_g, fields.sp.at("thl")->fld_g, fields.sp.at("qt")->fld_g,
            bs.thvref_g, bs.pref_g, bs.exnref_g,
            gd.istart,  gd.jstart, gd.kstart,
            gd.iend, gd.jend, gd.kcells,
            gd.icells, gd.ijcells);
        cuda_check_error();
    }
    else if (name == "b_h")
    {
        calc_buoyancy_g<<<gridGPU, blockGPU>>>(
            fld.fld_g, fields.sp.at("thl")->fld_g, fields.sp.at("qt")->fld_g,
            bs.thvrefh_g, bs.prefh_g, bs.exnrefh_g,
            gd.istart,  gd.jstart, gd.kstart,
            gd.iend, gd.jend, gd.kcells,
            gd.icells, gd.ijcells);
        cuda_check_error();
    }
    else if (name == "ql")
    {
        calc_liquid_water_g<<<gridGPU2, blockGPU2>>>(
            fld.fld_g, fields.sp.at("thl")->fld_g, fields.sp.at("qt")->fld_g,
            bs.exnref_g, bs.pref_g,
            gd.istart,  gd.jstart,  gd.kstart,
            gd.iend,    gd.jend,    gd.kend,
            gd.icells, gd.ijcells);
        cuda_check_error();
    }
    else if (name == "ql_h")
    {
        calc_liquid_water_h_g<<<gridGPU2, blockGPU2>>>(
            fld.fld_g, fields.sp.at("thl")->fld_g, fields.sp.at("qt")->fld_g,
            bs.exnrefh_g, bs.prefh_g,
            gd.istart,  gd.jstart,  gd.kstart,
            gd.iend,    gd.jend,    gd.kend,
            gd.icells, gd.ijcells);
        cuda_check_error();
    }
    else if (name == "N2")
    {
        calc_N2_g<<<gridGPU2, blockGPU2>>>(
            fld.fld_g, fields.sp.at("thl")->fld_g, bs.thvref_g, gd.dzi_g,
            gd.istart,  gd.jstart, gd.kstart,
            gd.iend,    gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
        cuda_check_error();
    }
    else
    {
        std::string msg = "get_thermo_field_g \"" + name + "\" not supported";
        throw std::runtime_error(msg);
    }

    if (cyclic)
        boundary_cyclic.exec_g(fld.fld_g);
}
#endif

#ifdef USECUDA
template<typename TF>
TF* Thermo_moist<TF>::get_basestate_fld_g(std::string name)
{
    // BvS TO-DO: change std::string to enum
    if (name == "pref")
        return bs.pref_g;
    else if (name == "prefh")
        return bs.prefh_g;
    else if (name == "exner")
        return bs.exnref_g;
    else if (name == "exnerh")
        return bs.exnrefh_g;
    else
    {
        std::string error_message = "Can not get basestate field \"" + name + "\" from thermo_moist";
        throw std::runtime_error(error_message);
    }
}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_moist<TF>::get_buoyancy_fluxbot_g(Field3d<TF>& bfield)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU (gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);

    calc_buoyancy_flux_bot_g<<<gridGPU, blockGPU>>>(
        bfield.flux_bot_g,
        fields.sp.at("thl")->fld_g, fields.sp.at("thl")->flux_bot_g,
        fields.sp.at("qt")->fld_g, fields.sp.at("qt")->flux_bot_g,
        bs.thvrefh_g, gd.kstart, gd.icells, gd.jcells,
        gd.icells, gd.ijcells);
    cuda_check_error();
}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_moist<TF>::get_buoyancy_surf_g(Field3d<TF>& bfield)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU (gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);

    calc_buoyancy_bot_g<<<gridGPU, blockGPU>>>(
        bfield.fld_g, bfield.fld_bot_g,
        fields.sp.at("thl")->fld_g, fields.sp.at("thl")->fld_bot_g,
        fields.sp.at("qt")->fld_g, fields.sp.at("qt")->fld_bot_g,
        bs.thvref_g, bs.thvrefh_g, gd.kstart, gd.icells, gd.jcells,
        gd.icells, gd.ijcells);
    cuda_check_error();

    calc_buoyancy_flux_bot_g<<<gridGPU, blockGPU>>>(
        bfield.flux_bot_g,
        fields.sp.at("thl")->fld_g, fields.sp.at("thl")->flux_bot_g,
        fields.sp.at("qt")->fld_g, fields.sp.at("qt")->flux_bot_g,
        bs.thvrefh_g, gd.kstart, gd.icells, gd.jcells,
        gd.icells, gd.ijcells);
    cuda_check_error();
}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_moist<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;
    auto output = fields.get_tmp_g();

    get_thermo_field_g(*output, "b", false);
    column.calc_column("b", output->fld_g, no_offset);

    get_thermo_field_g(*output, "ql", false);
    column.calc_column("ql", output->fld_g, no_offset);

    fields.release_tmp_g(output);
}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_moist<TF>::get_mask(Stats<TF>& stats, std::string mask_name)
{
    if (mask_name == "ql")
    {
        auto ql  = fields.get_tmp();
        auto qlh = fields.get_tmp();
        auto ql_g = fields.get_tmp_g();

        get_thermo_field_g(*ql_g, "ql", true);

        fields.backward_field_device_3d(ql->fld.data(), ql_g->fld_g);
        get_thermo_field_g(*ql_g, "ql_h", true);
        fields.backward_field_device_3d(qlh->fld.data(), ql_g->fld_g);

        stats.set_mask_thres(mask_name, *ql, *qlh, 0., Stats_mask_type::Plus);

        fields.release_tmp_g(ql_g);
        fields.release_tmp(ql);
        fields.release_tmp(qlh);
    }
    else if (mask_name == "qlcore")
    {
        auto ql  = fields.get_tmp();
        auto qlh = fields.get_tmp();
        auto tmp_g = fields.get_tmp_g();

        get_thermo_field_g(*tmp_g, "ql", true);
        fields.backward_field_device_3d(ql->fld.data(), tmp_g->fld_g);
        get_thermo_field_g(*tmp_g, "ql_h", true);
        fields.backward_field_device_3d(qlh->fld.data(), tmp_g->fld_g);

        stats.set_mask_thres(mask_name, *ql, *qlh, 0., Stats_mask_type::Plus);

        fields.release_tmp(ql);
        fields.release_tmp(qlh);

        auto b = fields.get_tmp();
        auto bh = fields.get_tmp();

        get_thermo_field_g(*tmp_g, "b", true);
        fields.backward_field_device_3d(b->fld.data(), tmp_g->fld_g);
        get_thermo_field_g(*tmp_g, "b_h", true);
        fields.backward_field_device_3d(bh->fld.data(), tmp_g->fld_g);

        field3d_operators.calc_mean_profile(b->fld_mean.data(), b->fld.data());
        field3d_operators.subtract_mean_profile(b->fld.data(), b->fld_mean.data());

        field3d_operators.calc_mean_profile(bh->fld_mean.data(), bh->fld.data());
        field3d_operators.subtract_mean_profile(bh->fld.data(), bh->fld_mean.data());

        stats.set_mask_thres(mask_name, *b, *bh, 0., Stats_mask_type::Plus);

        fields.release_tmp(b);
        fields.release_tmp(bh);
        fields.release_tmp_g(tmp_g);
    }
    else
    {
        std::string message = "Moist thermodynamics can not provide mask: \"" + mask_name +"\"";
        throw std::runtime_error(message);
    }
}
#endif

template class Thermo_moist<double>;
template class Thermo_moist<float>;
