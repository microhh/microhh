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

#ifndef BOUNDARY_SURFACE_KERNELS_H
#define BOUNDARY_SURFACE_KERNELS_H

#include <iostream>
#include "constants.h"
#include "monin_obukhov.h"

namespace Boundary_surface_kernels
{
    namespace fm = Fast_math;
    namespace most = Monin_obukhov;

    // Limits on Obukhov length:
    template<typename TF> constexpr TF zL_max = 10.;
    template<typename TF> constexpr TF zL_min = -1.e4;

    template<typename TF>
    void set_bc(
            TF* const restrict a, TF* const restrict agrad, TF* const restrict aflux,
            const Boundary_type sw, const TF aval, const TF visc, const TF offset,
            const int icells, const int jcells)
    {
        const int jj = icells;

        if (sw == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    a[ij] = aval - offset;
                }
        }
        else if (sw == Boundary_type::Neumann_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    agrad[ij] = aval;
                    aflux[ij] = -aval*visc;
                }
        }
        else if (sw == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    aflux[ij] = aval;
                    agrad[ij] = -aval/visc;
                }
        }
    }

    template<typename TF>
    void prepare_lut(
        float* const restrict zL_sl,
        float* const restrict f_sl,
        const TF z0m,
        const TF z0h,
        const TF zsl,
        const int nzL_lut,
        const Boundary_type mbcbot,
        const Boundary_type thermobc)
    {
        std::vector<TF> zL_tmp(nzL_lut);

        // Calculate the non-streched part between -5 to 10 z/L with 9/10 of the points,
        // and stretch up to -1e4 in the negative limit.
        // Alter next three values in case the range need to be changed.
        const TF zLrange_min = -5.;

        TF dzL = (zL_max<TF> - zLrange_min) / (9.*nzL_lut/10.-1.);
        zL_tmp[0] = -zL_max<TF>;
        for (int n=1; n<9*nzL_lut/10; ++n)
            zL_tmp[n] = zL_tmp[n-1] + dzL;

        // Stretch the remainder of the z/L values far down for free convection.
        const TF zLend = -(zL_min<TF> - zLrange_min);

        // Find stretching that ends up at the correct value using geometric progression.
        TF r  = 1.01;
        TF r0 = Constants::dhuge;
        while (std::abs( (r-r0)/r0 ) > 1.e-10)
        {
            r0 = r;
            r  = std::pow( 1. - (zLend/dzL)*(1.-r), (1./ (nzL_lut/10.) ) );
        }

        for (int n=9*nzL_lut/10; n<nzL_lut; ++n)
        {
            zL_tmp[n] = zL_tmp[n-1] + dzL;
            dzL *= r;
        }

        // Calculate the final array and delete the temporary array.
        for (int n=0; n<nzL_lut; ++n)
            zL_sl[n] = -zL_tmp[nzL_lut-n-1];

        // Calculate the evaluation function.
        if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Flux_type)
        {
            for (int n=0; n<nzL_lut; ++n)
                f_sl[n] = zL_sl[n] * std::pow(most::fm(zsl, z0m, zsl/zL_sl[n]), 3);
        }
        else if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Dirichlet_type)
        {
            for (int n=0; n<nzL_lut; ++n)
                f_sl[n] = zL_sl[n] * std::pow(most::fm(zsl, z0m, zsl/zL_sl[n]), 2) / most::fh(zsl, z0h, zsl/zL_sl[n]);
        }
    }

    template<typename TF>
    void calc_duvdz(
            TF* const restrict dudz,
            TF* const restrict dvdz,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict ubot,
            const TF* const restrict vbot,
            const TF* const restrict ufluxbot,
            const TF* const restrict vfluxbot,
            const TF* const restrict ustar,
            const TF* const restrict obuk,
            const TF* const restrict z0m,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int ijcells)
    {
        const int ii=1;
        const int jj=icells;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = ij + kstart*ijcells;

                const TF du_c = TF(0.5)*((u[ijk] - ubot[ij]) + (u[ijk+ii] - ubot[ij+ii]));
                const TF dv_c = TF(0.5)*((v[ijk] - vbot[ij]) + (v[ijk+jj] - vbot[ij+jj]));

                const TF ufluxbot = -du_c * ustar[ij] * most::fm(zsl, z0m[ij], obuk[ij]);
                const TF vfluxbot = -dv_c * ustar[ij] * most::fm(zsl, z0m[ij], obuk[ij]);

                dudz[ij] = -ufluxbot / (Constants::kappa<TF> * zsl * ustar[ij]) * most::phim(zsl/obuk[ij]);
                dvdz[ij] = -vfluxbot / (Constants::kappa<TF> * zsl * ustar[ij]) * most::phim(zsl/obuk[ij]);
            }
    }

    template<typename TF>
    void calc_dbdz(
            TF* const restrict dbdz,
            const TF* const restrict bfluxbot,
            const TF* const restrict ustar,
            const TF* const restrict obuk,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                dbdz[ij] = -bfluxbot[ij] / (Constants::kappa<TF> * zsl * ustar[ij]) * most::phih(zsl/obuk[ij]);
            }
    }

    template<typename TF>
    TF find_zL(
            const float* const restrict zL,
            const float* const restrict f,
            int& n, const float Ri)
    {
        // Determine search direction. All checks are at float accuracy.
        if ( (f[n]-Ri) > 0.f )
            while ( (f[n-1]-Ri) > 0.f && n > 0.f) { --n; }
        else
            while ( (f[n]-Ri) < 0.f && n < (nzL_lut-1) ) { ++n; }

        const TF zL0 = (n == 0 || n == nzL_lut-1) ? zL[n] : zL[n-1] + (Ri-f[n-1]) / (f[n]-f[n-1]) * (zL[n]-zL[n-1]);

        return zL0;
    }

    template<typename TF>
    TF calc_obuk_noslip_flux_lookup(
            const float* restrict zL, const float* restrict f,
            int& n,
            const TF du, const TF bfluxbot, const TF zsl)
    {
        // Calculate the appropriate Richardson number and reduce precision.
        const float Ri = -Constants::kappa<TF> * bfluxbot * zsl / fm::pow3(du);
        return zsl/find_zL<TF>(zL, f, n, Ri);
    }

    template<typename TF>
    TF calc_obuk_noslip_dirichlet_lookup(
            const float* restrict zL, const float* restrict f,
            int& n,
            const TF du, const TF db, const TF zsl)
    {
        // Calculate the appropriate Richardson number and reduce precision.
        const float Ri = Constants::kappa<TF> * db * zsl / fm::pow2(du);
        return zsl/find_zL<TF>(zL, f, n, Ri);
    }

    template<typename TF>
    TF calc_obuk_noslip_flux_iterative(
            TF L, const TF du, TF bfluxbot, const TF zsl, const TF z0m)
    {
        TF L0;
        TF Lstart, Lend;
        TF fx, fxdif;

        int m = 0;
        int nlim = 10;

        const TF Lmax = 1.e20;

        // Limiter max z/L stable conditions:
        const TF L_min_stable = zsl/zL_max<TF>;

        // Avoid bfluxbot to be zero
        if (bfluxbot >= 0.)
            bfluxbot = std::max(TF(Constants::dsmall), bfluxbot);
        else
            bfluxbot = std::min(-TF(Constants::dsmall), bfluxbot);

        // Allow for one restart
        while (m <= 1)
        {
            // If L and bfluxbot are of the same sign, or the last calculation did not converge,
            // the stability has changed and the procedure needs to be reset
            if (L*bfluxbot >= 0.)
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
            while (std::abs((L - L0)/L0) > 0.001 && n < nlim && std::abs(L) < Lmax)
            {
                L0     = L;
                fx     = zsl/L + Constants::kappa<TF>*zsl*bfluxbot / fm::pow3(du * most::fm(zsl, z0m, L));
                Lstart = L - 0.001*L;
                Lend   = L + 0.001*L;
                fxdif  = ( (zsl/Lend   + Constants::kappa<TF>*zsl*bfluxbot / fm::pow3(du * most::fm(zsl, z0m, Lend  )))
                         - (zsl/Lstart + Constants::kappa<TF>*zsl*bfluxbot / fm::pow3(du * most::fm(zsl, z0m, Lstart))) )
                       / (Lend - Lstart);
                L      = L - fx/fxdif;

                // Limit L for stable conditions (similar to the limits of the lookup table)
                if (L >= TF(0) && L < L_min_stable)
                    L = L_min_stable;
                ++n;
            }

            if (n < nlim && std::abs(L) < Lmax)
                // Convergence has been reached
                break;
            else
            {
                // Convergence has not been reached, procedure restarted once
                if (bfluxbot >= 0.)
                    L = -Constants::dsmall;
                else
                    L = L_min_stable;

                ++m;
                nlim = 200;
            }
        }

        if (m > 1)
            std::cout << "ERROR: convergence has not been reached in Obukhov length calculation" << std::endl;

        return zsl/std::min(std::max(zsl/L, zL_min<TF>), zL_max<TF>);
    }

    template<typename TF>
    TF calc_obuk_noslip_dirichlet_iterative(
            TF L, const TF du, TF db, const TF zsl, const TF z0m, const TF z0h)
    {
        TF L0;
        TF Lstart, Lend;
        TF fx, fxdif;

        int m = 0;
        int nlim = 10;

        const TF Lmax = 1.e20;

        // Avoid db to be zero
        if (db >= 0.)
            db = std::max(TF(Constants::dsmall), db);
        else
            db = std::min(-TF(Constants::dsmall), db);

        // Allow for one restart
        while (m <= 1)
        {
            // If L and db are of different sign, or the last calculation did not converge,
            // the stability has changed and the procedure needs to be reset
            if (L*db <= 0.)
            {
                nlim = 200;
                if(db >= 0.)
                    L = Constants::dsmall;
                else
                    L = -Constants::dsmall;
            }

            if (db >= 0.)
                L0 = Constants::dhuge;
            else
                L0 = -Constants::dhuge;

            int n = 0;

            // exit on convergence or on iteration count
            while (std::abs((L - L0)/L0) > 0.001 && n < nlim && std::abs(L) < Lmax)
            {
                L0     = L;
                fx     = zsl/L - Constants::kappa<TF>*zsl*db*most::fh(zsl, z0h, L) / fm::pow2(du * most::fm(zsl, z0m, L));
                Lstart = L - 0.001*L;
                Lend   = L + 0.001*L;
                fxdif  = ( (zsl/Lend   - Constants::kappa<TF>*zsl*db*most::fh(zsl, z0h, Lend)   / fm::pow2(du * most::fm(zsl, z0m, Lend  )))
                         - (zsl/Lstart - Constants::kappa<TF>*zsl*db*most::fh(zsl, z0h, Lstart) / fm::pow2(du * most::fm(zsl, z0m, Lstart))) )
                       / (Lend - Lstart);
                L      = L - fx/fxdif;
                ++n;
            }

            if (n < nlim && std::abs(L) < Lmax)
                // Convergence has been reached
                break;
            else
            {
                // Convergence has not been reached, procedure restarted once
                L = Constants::dsmall;
                ++m;
                nlim = 200;
            }
        }

        if (m > 1)
        {
            std::cout << "ERROR: no convergence obukhov iteration!" << std::endl;
            std::cout << "INPUT: du=" << du << ", db=" << db << ", z0m=" << z0m << ", z0h=" << z0h << ", OUTPUT: L=" << L <<  std::endl;
        }

        return zsl/std::min(std::max(zsl/L, zL_min<TF>), zL_max<TF>);
    }

    template<typename TF>
    void calc_ra(
            TF* const restrict ra,
            const TF* const restrict ustar,
            const TF* const restrict obuk,
            const TF* const restrict z0h,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                ra[ij]  = TF(1) / (ustar[ij] * most::fh(zsl, z0h[ij], obuk[ij]));
            }
    }
}
#endif
