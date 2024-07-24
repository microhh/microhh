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

#ifndef BOUNDARY_SURFACE_KERNELS_GPU_H
#define BOUNDARY_SURFACE_KERNELS_GPU_H

#include "monin_obukhov.h"
#include "boundary.h"

namespace Boundary_surface_kernels_g
{
    namespace most = Monin_obukhov;

    /* Calculate absolute wind speed */
    template<typename TF> __global__
    void calc_dutot_g(
            TF* const __restrict__ dutot,
            const TF* const __restrict__ u,
            const TF* const __restrict__ v,
            const TF* const __restrict__ ubot,
            const TF* const __restrict__ vbot,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ii  = 1;
            const int ii2 = 2;
            const int jj  = icells;
            const int jj2 = 2*icells;

            const int ij  = i + j*icells;
            const int ijk = i + j*icells + kstart*ijcells;
            const TF minval = 1.e-1;

            const TF u_filtered = TF(1./9) *
                ( TF(0.5)*u[ijk-ii-jj] + u[ijk-jj] + u[ijk+ii-jj] + TF(0.5)*u[ijk+ii2-jj]
                + TF(0.5)*u[ijk-ii   ] + u[ijk   ] + u[ijk+ii   ] + TF(0.5)*u[ijk+ii2   ]
                + TF(0.5)*u[ijk-ii+jj] + u[ijk+jj] + u[ijk+ii+jj] + TF(0.5)*u[ijk+ii2+jj] );

            const TF v_filtered = TF(1./9) *
                ( TF(0.5)*v[ijk-ii-jj] + v[ijk-ii] + v[ijk-ii+jj] + TF(0.5)*v[ijk-ii+jj2]
                + TF(0.5)*v[ijk   -jj] + v[ijk   ] + v[ijk   +jj] + TF(0.5)*v[ijk   +jj2]
                + TF(0.5)*v[ijk+ii-jj] + v[ijk+ii] + v[ijk+ii+jj] + TF(0.5)*v[ijk+ii+jj2] );

            const TF du2 = fm::pow2(u_filtered - TF(0.5)*(ubot[ij] + ubot[ij+ii]))
                         + fm::pow2(v_filtered - TF(0.5)*(vbot[ij] + vbot[ij+jj]));

            // Prevent the absolute wind gradient from reaching values less than 0.01 m/s,
            // otherwise evisc at k = kstart blows up
            dutot[ij] = fmax(std::pow(du2, TF(0.5)), minval);
        }
    }

    template<typename TF> __global__
    void calc_duvdz_mo_g(
            TF* const __restrict__ dudz,
            TF* const __restrict__ dvdz,
            const TF* const __restrict__ u,
            const TF* const __restrict__ v,
            const TF* const __restrict__ ubot,
            const TF* const __restrict__ vbot,
            const TF* const __restrict__ ufluxbot,
            const TF* const __restrict__ vfluxbot,
            const TF* const __restrict__ ustar,
            const TF* const __restrict__ obuk,
            const TF* const __restrict__ z0m,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        const int ii = 1;
        const int jj = icells;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*icells;
            const int ijk = i + j*icells + kstart*ijcells;

            const TF du_c = TF(0.5)*((u[ijk] - ubot[ij]) + (u[ijk+ii] - ubot[ij+ii]));
            const TF dv_c = TF(0.5)*((v[ijk] - vbot[ij]) + (v[ijk+jj] - vbot[ij+jj]));

            const TF ufluxbot = -du_c * ustar[ij] * most::fm(zsl, z0m[ij], obuk[ij]);
            const TF vfluxbot = -dv_c * ustar[ij] * most::fm(zsl, z0m[ij], obuk[ij]);

            dudz[ij] = -ufluxbot / (Constants::kappa<TF> * zsl * ustar[ij]) * most::phim(zsl/obuk[ij]);
            dvdz[ij] = -vfluxbot / (Constants::kappa<TF> * zsl * ustar[ij]) * most::phim(zsl/obuk[ij]);
        }
    }

    template<typename TF> __global__
    void calc_dbdz_mo_g(
            TF* const __restrict__ dbdz,
            const TF* const __restrict__ bfluxbot,
            const TF* const __restrict__ ustar,
            const TF* const __restrict__ obuk,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij = i + j*icells;
            dbdz[ij] = -bfluxbot[ij] / (Constants::kappa<TF> * zsl * ustar[ij]) * most::phih(zsl/obuk[ij]);
        }
    }

    template<typename TF> __device__
    TF find_zL_g(
            const float* const __restrict__ zL,
            const float* const __restrict__ f,
            int &n,
            const float Ri,
            const TF zsl)
    {
        // Determine search direction.
        if ((f[n]-Ri) > 0.f)
            while ( (f[n-1]-Ri) > 0.f && n > 0) { --n; }
        else
            while ( (f[n]-Ri) < 0.f && n < (nzL_lut-1) ) { ++n; }

        const TF zL0 = (n == 0 || n == nzL_lut-1) ? zL[n] : zL[n-1] + (Ri-f[n-1]) / (f[n]-f[n-1]) * (zL[n]-zL[n-1]);
        return zL0;
    }

    template<typename TF> __device__
    TF calc_obuk_noslip_dirichlet_lookup_g(
            const float* const __restrict__ zL_lut,
            const float* const __restrict__ f_lut,
            int& n,
            const TF du,
            const TF db,
            const TF zsl)
    {
        // Calculate the appropriate Richardson number.
        const float Ri = Constants::kappa<TF> * db * zsl / fm::pow2(du);
        const TF zL = find_zL_g(zL_lut, f_lut, n, Ri, zsl);
        return zsl/zL;
    }


    template<typename TF> __device__
    TF calc_obuk_noslip_flux_lookup_g(
            const float* const __restrict__ zL_lut,
            const float* const __restrict__ f_lut,
            int& n,
            const TF du,
            const TF bfluxbot,
            const TF zsl)
    {
        // Calculate the appropriate Richardson number.
        const TF Ri = -Constants::kappa<TF> * bfluxbot * zsl / fm::pow3(du);
        const TF zL = find_zL_g(zL_lut, f_lut, n, Ri, zsl);
        return zsl/zL;
    }


    template<typename TF> __device__
    TF calc_obuk_noslip_flux_iterative_g(
            TF L, const TF du, TF bfluxbot, const TF zsl, const TF z0m)
    {
        TF L0;
        TF Lstart, Lend;
        TF fx, fxdif;

        int m = 0;
        int nlim = 10;

        const TF Lmax = 1.e16;

        // Limiter max z/L stable conditions:
        const TF L_min_stable = zsl/Constants::zL_max<TF>;

        // Avoid bfluxbot to be zero
        if (bfluxbot >= 0.)
            bfluxbot = fmax(TF(Constants::dsmall), bfluxbot);
        else
            bfluxbot = fmin(-TF(Constants::dsmall), bfluxbot);

        // Allow for one restart
        while (m <= 1)
        {
            // If L and bfluxbot are of the same sign, or the last calculation did not converge,
            // the stability has changed and the procedure needs to be reset
            if (m == 1 || L*bfluxbot >= 0.)
            {
                nlim = 200;
                if (bfluxbot >= 0.)
                    L = -Constants::dsmall;
                else
                    L = L_min_stable;
            }

            if (bfluxbot >= 0.)
                L0 = -Constants::dhuge;
            else
                L0 = Constants::dhuge;

            int n = 0;

            // Exit on convergence or on iteration count
            while (abs((L - L0)/L0) > 0.001 && n < nlim && abs(L) < Lmax)
            {
                L0     = L;
                fx     = zsl/L + Constants::kappa<TF>*zsl*bfluxbot / fm::pow3(du * most::fm(zsl, z0m, L));
                Lstart = L - 0.001*L;
                Lend   = L + 0.001*L;
                fxdif  = ( (zsl/Lend   + Constants::kappa<TF>*zsl*bfluxbot / fm::pow3(du * most::fm(zsl, z0m, Lend  )))
                         - (zsl/Lstart + Constants::kappa<TF>*zsl*bfluxbot / fm::pow3(du * most::fm(zsl, z0m, Lstart))) )
                         / (Lend - Lstart);

                if (abs(fxdif) < TF(1e-16))
                    break;

                L = L - fx/fxdif;

                // Limit L for stable conditions (similar to the limits of the lookup table)
                if (L >= TF(0) && L < L_min_stable)
                    L = L_min_stable;
                ++n;
            }

            if (n < nlim && abs(L) < Lmax)
                break; // Convergence has been reached
            else
            {
                // Convergence has not been reached, procedure restarted once
                ++m;
                nlim = 200;
            }
        }

        if (m > 1)
        {
            if (abs(L) > Lmax)
                L = -copysign(TF(1), bfluxbot) * Lmax;  // Non-converged in neutral regime; return Lmax with correct sign.
            else
                L = L_min_stable;  // Non-converged in stable regime.

            #ifdef PRINT_RIB_L_ERRORS
            printf("ERROR: no convergence Rib->L: du=%f, B0=%f, z0m=%f | returning L=%f\n", du, bfluxbot, z0m, L);
            #endif
        }

        return zsl/fmin(fmax(zsl/L, Constants::zL_min<TF>), Constants::zL_max<TF>);
    }


    template<typename TF> __device__
    TF calc_obuk_noslip_dirichlet_iterative_g(
            TF L, const TF du, TF db, const TF zsl, const TF z0m, const TF z0h)
    {
        TF L0;
        TF Lstart, Lend;
        TF fx, fxdif;

        int m = 0;
        int nlim = 10;
        const TF Lmax = 1.e16;

        const TF L_min_stable = zsl / Constants::zL_max<TF>;

        // The solver does not have a solution for large Ri numbers,
        // i.e. the `fx` equation below has no zero crossing.
        // The limit of 0.13 typically results in a minimum (positive) L of ~1
        const TF Ri = Constants::kappa<TF> * db * zsl / fm::pow2(du);
        if (Ri > TF(0.13))
            return L_min_stable;

        // Avoid db to be zero
        if (db >= 0.)
            db = fmax(TF(Constants::dsmall), db);
        else
            db = fmin(-TF(Constants::dsmall), db);

        // Allow for one restart
        while (m <= 1)
        {
            // If L and db are of different sign, or the last calculation did not converge,
            // the stability has changed and the procedure needs to be reset
            if (m == 1 || L*db <= 0.)
            {
                nlim = 200;
                if(db >= 0.)
                    L = L_min_stable;
                else
                    L = -Constants::dsmall;
            }

            if (db >= 0.)
                L0 = Constants::dhuge;
            else
                L0 = -Constants::dhuge;

            int n = 0;

            // Exit on convergence or on iteration count
            while (abs((L - L0)/L0) > 0.001 && n < nlim && abs(L) < Lmax)
            {
                L0     = L;
                fx     = zsl/L - Constants::kappa<TF>*zsl*db*most::fh(zsl, z0h, L) / fm::pow2(du * most::fm(zsl, z0m, L));
                Lstart = L - 0.001*L;
                Lend   = L + 0.001*L;
                fxdif  = ( (zsl/Lend   - Constants::kappa<TF>*zsl*db*most::fh(zsl, z0h, Lend)   / fm::pow2(du * most::fm(zsl, z0m, Lend  )))
                         - (zsl/Lstart - Constants::kappa<TF>*zsl*db*most::fh(zsl, z0h, Lstart) / fm::pow2(du * most::fm(zsl, z0m, Lstart))) )
                       / (Lend - Lstart);

                if (abs(fxdif) < TF(1e-16))
                    break;

                L = L - fx/fxdif;

                if (L >= TF(0) && L < L_min_stable)
                    L = L_min_stable;

                ++n;
            }

            if (n < nlim && abs(L) < Lmax)
                // Convergence has been reached
                break;
            else
            {
                // Convergence has not been reached, procedure restarted once
                ++m;
                nlim = 200;
            }
        }

        if (m > 1)
        {
            if (abs(L) > Lmax)
                L = copysign(TF(1), db) * Lmax; // Non-converged in neutral regime; return Lmax with correct sign.
            else
                L = L_min_stable; // Non-converged in stable regime.

            #ifdef PRINT_RIB_L_ERRORS
            printf("ERROR: no convergence Rib->L: du=%f, db=%f, z0m=%f, z0h=%f | returning L=%f\n", du, db, z0m, z0h, L);
            #endif
        }

        // Limits same as LUT solver:
        return zsl/fmin(fmax(zsl/L, Constants::zL_min<TF>), Constants::zL_max<TF>);
    }
}
#endif
