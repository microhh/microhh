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
}
#endif
