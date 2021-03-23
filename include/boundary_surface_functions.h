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

#ifndef BOUNDARY_SURFACE_FUNCTIONS_H
#define BOUNDARY_SURFACE_FUNCTIONS_H

#include "constants.h"
#include "monin_obukhov.h"

namespace Boundary_surface_functions
{
    namespace fm = Fast_math;
    namespace most = Monin_obukhov;

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
                dbdz[ij] = -bfluxbot[ij]/(Constants::kappa<TF>*zsl*ustar[ij])*most::phih(zsl/obuk[ij]);
            }
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
        const TF zL_max_stable = 10.;
        const TF L_min_stable = zsl/zL_max_stable;

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

        return L;
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

        return L;
    }
}
#endif
