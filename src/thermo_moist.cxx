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
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "thermo_moist.h"
#include "diff.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "netcdf_interface.h"
#include "model.h"
#include "stats.h"
#include "master.h"
#include "cross.h"
#include "dump.h"
#include "column.h"
#include "thermo_moist_functions.h"
#include "timeloop.h"
#include "field3d_operators.h"

using Finite_difference::O2::interp2;
using namespace Constants;
using namespace Thermo_moist_functions;

namespace
{
    // Help function(s) to switch between the different NetCDF data types
    template<typename TF> TF netcdf_fp_fillvalue();
    template<> double netcdf_fp_fillvalue<double>() { return NC_FILL_DOUBLE; }
    template<> float  netcdf_fp_fillvalue<float>()  { return NC_FILL_FLOAT; }

    template<typename TF>
    void calc_top_and_bot(TF* restrict thl0, TF* restrict qt0,
                          const TF* const z, const TF* const zh,
                          const TF* const dzhi,
                          const int kstart, const int kend)
    {
        // Calculate surface and model top values thl and qt
        TF thl0s, qt0s, thl0t, qt0t;
        thl0s = thl0[kstart] - z[kstart]*(thl0[kstart+1]-thl0[kstart])*dzhi[kstart+1];
        qt0s  = qt0[kstart]  - z[kstart]*(qt0[kstart+1] -qt0[kstart] )*dzhi[kstart+1];
        thl0t = thl0[kend-1] + (zh[kend]-z[kend-1])*(thl0[kend-1]-thl0[kend-2])*dzhi[kend-1];
        qt0t  = qt0[kend-1]  + (zh[kend]-z[kend-1])*(qt0[kend-1]- qt0[kend-2] )*dzhi[kend-1];

        // Set the ghost cells for the reference temperature and moisture
        thl0[kstart-1]  = TF(2.)*thl0s - thl0[kstart];
        thl0[kend]      = TF(2.)*thl0t - thl0[kend-1];
        qt0[kstart-1]   = TF(2.)*qt0s  - qt0[kstart];
        qt0[kend]       = TF(2.)*qt0t  - qt0[kend-1];
    }

    template<typename TF>
    void calc_buoyancy_tend_2nd(
            TF* restrict wt, TF* restrict thl, TF* restrict qt,
            TF* restrict ph, TF* restrict thlh, TF* restrict qth,
            TF* restrict ql, TF* restrict qi, TF* restrict thvrefh,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        #pragma omp parallel for
        for (int k=kstart+1; k<kend; k++)
        {
            const TF exnh = exner(ph[k]);
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ij  = i + j*jj;
                    thlh[ij] = interp2(thl[ijk-kk], thl[ijk]);
                    qth[ij]  = interp2(qt[ijk-kk], qt[ijk]);
                }

            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ij  = i + j*jj;
                    Struct_sat_adjust<TF> ssa = sat_adjust(thlh[ij], qth[ij], ph[k], exnh);
                    ql[ij] = ssa.ql;
                    qi[ij] = ssa.qi;
                }

            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ij  = i + j*jj;
                    wt[ijk] += buoyancy(exnh, thlh[ij], qth[ij], ql[ij], qi[ij], thvrefh[k]);
                }
        }
    }

    template<typename TF>
    void calc_buoyancy(
            TF* restrict b, TF* restrict thl, TF* restrict qt,
            TF* restrict p, TF* restrict ql, TF* restrict qi, TF* restrict thvref,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int kcells, const int jj, const int kk)
    {
        #pragma omp parallel for
        for (int k=0; k<kcells; k++)
        {
            const TF ex = exner(p[k]);
            if (k >= kstart && k < kend)
            {
                for (int j=jstart; j<jend; j++)
                    #pragma ivdep
                    for (int i=istart; i<iend; i++)
                    {
                        const int ijk = i + j*jj + k*kk;
                        Struct_sat_adjust<TF> ssa = sat_adjust(thl[ijk], qt[ijk], p[k], ex);
                        ql[ijk] = ssa.ql;
                        qi[ijk] = ssa.qi;
                    }
            }
            else
            {
                for (int j=jstart; j<jend; j++)
                    #pragma ivdep
                    for (int i=istart; i<iend; i++)
                    {
                        const int ijk  = i + j*jj+k*kk;
                        ql[ijk] = 0.;
                        qi[ijk] = 0.;
                    }

            }
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    b[ijk] = buoyancy(ex, thl[ijk], qt[ijk], ql[ijk], qi[ijk], thvref[k]);
                }
        }
    }

    template<typename TF>
    void calc_buoyancy_h(TF* restrict bh, TF* restrict thl, TF* restrict qt,
                         TF* restrict ph, TF* restrict thvrefh, TF* restrict thlh, TF* restrict qth,
                         TF* restrict ql, TF* restrict qi,
                         const int istart, const int iend,
                         const int jstart, const int jend,
                         const int kstart, const int kend,
                         const int jj, const int kk)
    {
        using Finite_difference::O2::interp2;

        for (int k=kstart; k<kend+1; k++)
        {
            const TF exnh = exner(ph[k]);

            if (k>=kstart)
            {
                for (int j=jstart; j<jend; j++)
                    #pragma ivdep
                    for (int i=istart; i<iend; i++)
                    {
                        const int ijk = i + j*jj + k*kk;
                        const int ij  = i + j*jj;

                        thlh[ij] = interp2(thl[ijk-kk], thl[ijk]);
                        qth[ij]  = interp2(qt[ijk-kk], qt[ijk]);
                    }

                for (int j=jstart; j<jend; j++)
                    #pragma ivdep
                    for (int i=istart; i<iend; i++)
                    {
                        const int ij  = i + j*jj;
                        Struct_sat_adjust<TF> ssa = sat_adjust(thlh[ij], qth[ij], ph[k], exnh);
                        ql[ij] = ssa.ql;
                        qi[ij] = ssa.qi;
                    }
            }
            else
            {
                for (int j=jstart; j<jend; j++)
                    #pragma ivdep
                    for (int i=istart; i<iend; i++)
                    {
                        const int ij  = i + j*jj;
                        ql[ij] = 0.;
                        qi[ij] = 0.;
                    }
            }
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ij  = i + j*jj;

                    bh[ijk] = buoyancy(exnh, thlh[ij], qth[ij], ql[ij], qi[ij], thvrefh[k]);
                }
        }
    }

    template<typename TF>
    void calc_liquid_water(TF* restrict ql, TF* restrict thl, TF* restrict qt, TF* restrict p,
                           const int istart, const int iend,
                           const int jstart, const int jend,
                           const int kstart, const int kend,
                           const int jj, const int kk)
    {
        // Calculate the ql field
        #pragma omp parallel for
        for (int k=kstart; k<kend; k++)
        {
            const TF ex = exner(p[k]);
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    ql[ijk] = sat_adjust(thl[ijk], qt[ijk], p[k], ex).ql;
                }
        }
    }

    template<typename TF>
    void calc_saturated_water_vapor(
            TF* restrict qsat, TF* restrict thl, TF* restrict qt, TF* restrict p,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        // Calculate the ql field
        #pragma omp parallel for
        for (int k=kstart; k<kend; k++)
        {
            const TF ex = exner(p[k]);
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    qsat[ijk] = sat_adjust(thl[ijk], qt[ijk], p[k], ex).qs;
                }
        }
    }

    template<typename TF>
    void calc_relative_humidity(
            TF* restrict rh, TF* restrict thl, TF* restrict qt, TF* restrict p,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        // Calculate the ql field
        #pragma omp parallel for
        for (int k=kstart; k<kend; k++)
        {
            const TF ex = exner(p[k]);
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    rh[ijk] = std::min( qt[ijk] / sat_adjust(thl[ijk], qt[ijk], p[k], ex).qs, TF(1.));
                }
        }
    }

    template<typename TF>
    void calc_w500hpa(
            TF* restrict w500, const TF* restrict w, const TF* restrict ph,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        // Find the first index that is above 500 hPa.
        int index500 = 0;
        for (int k=kstart; k<kend+1; ++k)
        {
            if (ph[k] <= 5.e4)
            {
                index500 = k;
                break;
            }
        }
        if (index500 == 0 || index500 == kend)
            throw std::runtime_error("calc_w500hPa cannot find pressure level inside of domain");

        #pragma omp parallel for
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + index500*kk;
                w500[ij] = w[ijk];
            }
    }

    template<typename TF>
    void calc_qlqicore_max_thv_prime(
            TF* const restrict thv_prime,
            const TF* restrict qlqi,
            const TF* const restrict thv,
            const TF* const restrict thv_mean,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*icells;
                thv_prime[ij] = -Constants::dbig;

                bool has_cloud = false;
                for (int k=kstart; k<kend; k++)
                {
                    const int ijk = ij + k*ijcells;

                    if (qlqi[ijk] > 0)
                    {
                        has_cloud = true;

                        if (thv[ijk]-thv_mean[k] > thv_prime[ij])
                            thv_prime[ij] = thv[ijk]-thv_mean[k];
                    }
                }

                if (!has_cloud)
                    thv_prime[ij] = netcdf_fp_fillvalue<TF>();
            }
    }

    template<typename TF>
    void calc_liquid_water_h(TF* restrict qlh, TF* restrict thl,  TF* restrict qt,
                             TF* restrict ph, TF* restrict thlh, TF* restrict qth,
                             const int istart, const int iend,
                             const int jstart, const int jend,
                             const int kstart, const int kend,
                             const int jj, const int kk)
    {
        using Finite_difference::O2::interp2;

        for (int k=kstart+1; k<kend+1; k++)
        {
            const TF exnh = exner(ph[k]);

            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ij  = i + j*jj;

                    thlh[ij] = interp2(thl[ijk-kk], thl[ijk]);
                    qth[ij]  = interp2(qt[ijk-kk], qt[ijk]);
                }

            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ij  = i + j*jj;
                    const int ijk  = i + j*jj+k*kk;

                    qlh[ijk] = sat_adjust(thlh[ij], qth[ij], ph[k], exnh).ql;
                }
        }

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk  = i + j*jj+kstart*kk;
                qlh[ijk] = 0.;
            }
    }

    template<typename TF>
    void calc_ice(
            TF* restrict qi, TF* restrict thl, TF* restrict qt, TF* restrict p,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        // Calculate the ql field
        #pragma omp parallel for
        for (int k=kstart; k<kend; k++)
        {
            const TF ex = exner(p[k]);
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    qi[ijk] = sat_adjust(thl[ijk], qt[ijk], p[k], ex).qi;
                }
        }
    }

    template<typename TF>
    void calc_condensate(
            TF* restrict qc, TF* restrict thl, TF* restrict qt, TF* restrict p,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        // Calculate the ql field
        #pragma omp parallel for
        for (int k=kstart; k<kend; k++)
        {
            const TF ex = exner(p[k]);
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    qc[ijk] = std::max(qt[ijk] - sat_adjust(thl[ijk], qt[ijk], p[k], ex).qs, TF(0.));
                }
        }
    }

    template<typename TF>
    void calc_N2(TF* restrict N2, const TF* const restrict thl, const TF* const restrict dzi, TF* restrict thvref,
                 const int istart, const int iend,
                 const int jstart, const int jend,
                 const int kstart, const int kend,
                 const int jj, const int kk)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    N2[ijk] = grav<TF>/thvref[k]*TF(0.5)*(thl[ijk+kk] - thl[ijk-kk])*dzi[k];
                }
    }

    template<typename TF>
    void calc_T(TF* const restrict T, const TF* const restrict thl, const TF* const restrict qt,
                const TF* const restrict pref, const TF* const restrict exnref,
                const int istart, const int iend,
                const int jstart, const int jend,
                const int kstart, const int kend,
                const int jj, const int kk)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj+ k*kk;
                    T[ijk] = sat_adjust(thl[ijk], qt[ijk], pref[k], exnref[k]).t;
                }
        }
    }

    template<typename TF>
    void calc_path(
        TF* const restrict path,
        const TF* const restrict fld,
        const TF* const restrict rhoref,
        const TF* const restrict dz,
        const int istart, const int iend,
        const int jstart, const int jend,
        const int kstart, const int kend,
        const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*icells;
                path[ij] = TF(0);
            }

        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*icells;
                    const int ijk = ij + k*ijcells;

                    path[ij] += rhoref[k] * fld[ijk] * dz[k];
                }
        }
    }

    template<typename TF>
    void calc_T_h(TF* restrict Th, TF* restrict thl,  TF* restrict qt,
                  TF* restrict ph, TF* restrict thlh, TF* restrict qth, TF* restrict ql,
                  const int istart, const int iend,
                  const int jstart, const int jend,
                  const int kstart, const int kend,
                  const int jj, const int kk)
    {
        using Finite_difference::O2::interp2;

        for (int k=kstart+1; k<kend+1; k++)
        {
            const TF exnh = exner(ph[k]);
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ij  = i + j*jj;

                    thlh[ij] = interp2(thl[ijk-kk], thl[ijk]);
                    qth[ij]  = interp2(qt[ijk-kk], qt[ijk]);
                }

            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ij  = i + j*jj;
                    const int ijk  = i + j*jj+k*kk;

                    Th[ijk] = sat_adjust(thlh[ij], qth[ij], ph[k], exnh).t;
                }
        }
    }

    template<typename TF>
    void calc_thv(
            TF* const restrict thv, const TF* const restrict thl,
            const TF* const restrict qt, const TF* const restrict p,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart_minus_one, const int kend_plus_one,
            const int jj, const int kk)
    {
        // Note: this function is only called from statistics,
        // using `kstart-1` and `kend+1` as bounds, since the ghost
        // cells are needed to calculate the surface gradients.
        // Make sure that we don't `sat_adjust()` the ghost cells...

        // Calculate the thv field: ghost cells
        #pragma omp parallel for
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj + kstart_minus_one*kk;
                thv[ijk] = virtual_temperature_no_ql(thl[ijk], qt[ijk]);
            }

        #pragma omp parallel for
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj + (kend_plus_one-1)*kk;
                thv[ijk] = virtual_temperature_no_ql(thl[ijk], qt[ijk]);
            }

        // Calculate the thv field: interior
        #pragma omp parallel for
        for (int k=kstart_minus_one+1; k<kend_plus_one-1; k++)
        {
            const TF ex = exner(p[k]);
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;

                    Struct_sat_adjust<TF> ssa = sat_adjust(thl[ijk], qt[ijk], p[k], ex);
                    thv[ijk] = virtual_temperature(ex, thl[ijk], qt[ijk], ssa.ql, ssa.qi);
                }
        }
    }

    template<typename TF>
    void calc_T_bot(TF* const restrict T_bot, const TF* const restrict th,
                    const TF* const restrict exnrefh, const TF* const restrict threfh,
                    const int istart, const int iend, const int jstart, const int jend, const int kstart,
                    const int jj, const int kk)
    {
        using Finite_difference::O2::interp2;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                T_bot[ij] = exnrefh[kstart]*threfh[kstart] + (interp2(th[ijk-kk], th[ijk]) - threfh[kstart]);
            }
    }

    template<typename TF>
    void calc_buoyancy_bot(TF* restrict b,      TF* restrict bbot,
                           TF* restrict thl,    TF* restrict thlbot,
                           TF* restrict qt,     TF* restrict qtbot,
                           TF* restrict thvref, TF* restrict thvrefh,
                           const int icells, const int jcells,
                           const int ijcells, const int kstart)
    {
        // assume no liquid water at the lowest model level
        for (int j=0; j<jcells; j++)
            #pragma ivdep
            for (int i=0; i<icells; i++)
            {
                const int ij  = i + j*icells;
                const int ijk = i + j*icells + kstart*ijcells;
                bbot[ij ] = buoyancy_no_ql(thlbot[ij], qtbot[ij], thvrefh[kstart]);
                b   [ijk] = buoyancy_no_ql(thl[ijk], qt[ijk], thvref[kstart]);
            }
    }

    template<typename TF>
    void calc_buoyancy_bot(
            TF* const restrict b_bot,
            const TF* const restrict thl_bot,
            const TF* const restrict qt_bot,
            const TF* const restrict thvrefh,
            const int icells, const int jcells, const int kstart)
    {
        // assume no liquid water at the lowest model level
        for (int j=0; j<jcells; j++)
            #pragma ivdep
            for (int i=0; i<icells; i++)
            {
                const int ij  = i + j*icells;
                b_bot[ij] = buoyancy_no_ql(thl_bot[ij], qt_bot[ij], thvrefh[kstart]);
            }
    }

    template<typename TF>
    void calc_buoyancy_fluxbot(TF* restrict bfluxbot, TF* restrict thl, TF* restrict thlfluxbot,
                               TF* restrict qt, TF* restrict qtfluxbot, TF* restrict thvrefh,
                               const int icells, const int jcells, const int kstart,
                               const int ijcells)
    {

        // Assume no liquid water at the lowest model level.
        // Pass the temperature and moisture of the first model level, because the surface values are
        // unknown before the surface layer solver.
        for (int j=0; j<jcells; j++)
            #pragma ivdep
            for (int i=0; i<icells; i++)
            {
                const int ij  = i + j*icells;
                const int ijk = i + j*icells + kstart*ijcells;
                bfluxbot[ij] = buoyancy_flux_no_ql(thl[ijk], thlfluxbot[ij], qt[ijk], qtfluxbot[ij], thvrefh[kstart]);
            }
    }

    template<typename TF>
    void calc_thv_fluxbot(
            TF* const restrict thv_fluxbot,
            const TF* const restrict thl, const TF* const restrict thl_fluxbot,
            const TF* const restrict qt,  const TF* const restrict qt_fluxbot,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int ijcells)
    {
        // Assume no liquid water at the lowest model level.
        // Pass the temperature and moisture of the first model level, because the surface values are
        // unknown before the surface layer solver.
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij  = i + j*icells;
                const int ijk = ij + kstart*ijcells;
                thv_fluxbot[ij] = virtual_temperature_flux_no_ql(
                        thl[ijk], thl_fluxbot[ij], qt[ijk], qt_fluxbot[ij]);
            }
    }

    template<typename TF>
    void calc_buoyancy_tend_4th(TF* restrict wt, TF* restrict thl,  TF* restrict qt,
                                TF* restrict ph, TF* restrict thlh, TF* restrict qth,
                                TF* restrict ql, TF* restrict thvrefh,
                                const int istart, const int iend,
                                const int jstart, const int jend,
                                const int kstart, const int kend,
                                const int icells, const int ijcells)
    {
        const int jj  = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        for (int k=kstart+1; k<kend; k++)
        {
            const TF exnh = exner(ph[k]);
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk1;
                    const int ij  = i + j*jj;

                    thlh[ij]    = interp4c(thl[ijk-kk2], thl[ijk-kk1], thl[ijk], thl[ijk+kk1]);
                    qth[ij]     = interp4c(qt[ijk-kk2],  qt[ijk-kk1],  qt[ijk],  qt[ijk+kk1]);
                    const TF tl = thlh[ij] * exnh;

                    // Calculate first estimate of ql using Tl
                    // if ql(Tl)>0, saturation adjustment routine needed
                    ql[ij]  = qth[ij]-qsat(ph[k], tl);
                }

            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ij = i + j*jj;

                    if (ql[ij] > 0)   // already doesn't vectorize because of iteration in sat_adjust()
                        ql[ij] = sat_adjust(thlh[ij], qth[ij], ph[k], exnh).ql;
                    else
                        ql[ij] = 0.;
                }

            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk1;
                    const int ij  = i + j*jj;

                    wt[ijk] += buoyancy(exnh, thlh[ij], qth[ij], ql[ij], thvrefh[k]);
                }
        }
    }

    template<typename TF>
    int calc_zi(const TF* const restrict fldmean, const int kstart, const int kend, const int plusminus)
    {
        TF maxgrad = 0.;
        TF grad = 0.;
        int kinv = kstart;
        for (int k=kstart+1; k<kend; ++k)
        {
            grad = plusminus * (fldmean[k] - fldmean[k-1]);
            if (grad > maxgrad)
            {
                maxgrad = grad;
                kinv = k;
            }
        }
        return kinv;
    }

    template<typename TF>
    void calc_radiation_fields(
            TF* restrict T, TF* restrict T_h, TF* restrict vmr_h2o,
            TF* restrict clwp, TF* restrict ciwp, TF* restrict T_sfc,
            TF* restrict thlh, TF* restrict qth,
            const TF* restrict thl, const TF* restrict qt, const TF* restrict thl_bot,
            const TF* restrict p, const TF* restrict ph,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int igc, const int jgc, const int kgc,
            const int jj, const int kk,
            const int jj_nogc, const int kk_nogc)
    {
        // This routine strips off the ghost cells, because of the data handling in radiation.
        using Finite_difference::O2::interp2;

        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
        {
            const TF ex = exner(p[k]);
            const TF dpg = (ph[k] - ph[k+1]) / Constants::grav<TF>;
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ijk_nogc = (i-igc) + (j-jgc)*jj_nogc + (k-kgc)*kk_nogc;
                    const Struct_sat_adjust<TF> ssa = sat_adjust(thl[ijk], qt[ijk], p[k], ex);

                    clwp[ijk_nogc] = ssa.ql * dpg;
                    ciwp[ijk_nogc] = ssa.qi * dpg;

                    const TF qv = qt[ijk] - ssa.ql - ssa.qi;
                    vmr_h2o[ijk_nogc] = qv / (ep<TF> - ep<TF>*qv);

                    T[ijk_nogc] = ssa.t;
                }
        }

        for (int k=kstart; k<kend+1; ++k)
        {
            const TF exnh = exner(ph[k]);
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + k*kk;
                    thlh[ij] = interp2(thl[ijk-kk], thl[ijk]);
                    qth [ij] = interp2(qt [ijk-kk], qt [ijk]);
                }

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    const int ijk_nogc = (i-igc) + (j-jgc)*jj_nogc + (k-kgc)*kk_nogc;
                    T_h[ijk_nogc] = sat_adjust(thlh[ij], qth[ij], ph[k], exnh).t;
                }
        }

        // Calculate surface temperature (assuming no liquid water)
        const TF exn_bot = exner(ph[kstart]);
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*jj;
                const int ij_nogc = (i-igc) + (j-jgc)*jj_nogc;

                T_sfc[ij_nogc] = thl_bot[ij] * exn_bot;
            }
    }

    template<typename TF>
    void calc_radiation_columns(
            TF* const restrict T, TF* const restrict T_h, TF* const restrict vmr_h2o,
            TF* const restrict clwp, TF* const restrict ciwp, TF* const restrict T_sfc,
            const TF* const restrict thl, const TF* const restrict qt, const TF* const restrict thl_bot,
            const TF* const restrict p, const TF* const restrict ph,
            const int* const col_i, const int* const col_j,
            const int n_cols,
            const int kgc, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        // This routine strips off the ghost cells, because of the data handling in radiation.
        using Finite_difference::O2::interp2;

        const int ktot = kend-kstart;

        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
        {
            const TF ex = exner(p[k]);
            const TF dpg = (ph[k] - ph[k+1]) / Constants::grav<TF>;

            #pragma ivdep
            for (int n=0; n<n_cols; ++n)
            {
                const int i = col_i[n];
                const int j = col_j[n];

                const int ijk = i + j*icells + k*ijcells;
                const int ijk_out = n + (k-kgc)*n_cols;

                const Struct_sat_adjust<TF> ssa = sat_adjust(thl[ijk], qt[ijk], p[k], ex);

                clwp[ijk_out] = ssa.ql * dpg;
                ciwp[ijk_out] = ssa.qi * dpg;

                const TF qv = qt[ijk] - ssa.ql - ssa.qi;
                vmr_h2o[ijk_out] = qv / (ep<TF> - ep<TF>*qv);
                T[ijk_out] = ssa.t;
            }
        }

        for (int k=kstart; k<kend+1; ++k)
        {
            const TF exnh = exner(ph[k]);

            #pragma ivdep
            for (int n=0; n<n_cols; ++n)
            {
                const int i = col_i[n];
                const int j = col_j[n];

                const int ijk = i + j*icells + k*ijcells;
                const int ijk_out = n + (k-kgc)*n_cols;

                const TF thlh = interp2(thl[ijk-ijcells], thl[ijk]);
                const TF qth  = interp2(qt [ijk-ijcells], qt [ijk]);

                T_h[ijk_out] = sat_adjust(thlh, qth, ph[k], exnh).t;
            }
        }

        // Calculate surface temperature (assuming no liquid water)
        const TF exn_bot = exner(ph[kstart]);

        #pragma ivdep
        for (int n=0; n<n_cols; ++n)
        {
            const int i = col_i[n];
            const int j = col_j[n];

            const int ij = i + j*icells;
            const int ij_out = n;

            T_sfc[ij_out] = thl_bot[ij] * exn_bot;
        }
    }

    template<typename TF>
    void calc_land_surface_fields(
            TF* const restrict T_bot,
            TF* const restrict T_a,
            TF* const restrict vpd,
            TF* const restrict qs,
            TF* const restrict dqsdT,
            const TF* const restrict thl_bot,
            const TF* const restrict thl,
            const TF* const restrict qt,
            const TF* const restrict exner,
            const TF* const restrict exnerh,
            const TF* const restrict p,
            const TF* const restrict ph,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int ijcells)
    {
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij  = i + j*icells;
                const int ijk = i + j*icells + kstart*ijcells;

                // Saturation adjustment for first model level
                Struct_sat_adjust<TF> sa = sat_adjust(
                        thl[ijk], qt[ijk], p[kstart], exner[kstart]);

                T_bot[ij] = exnerh[kstart] * thl_bot[ij];   // Assuming no ql
                T_a[ij]   = sa.t;

                // Vapor pressure deficit first model level
                const TF es = esat(sa.t);
                const TF e = qt[ijk]/sa.qs * es;
                vpd[ij] = es-e;

                // qsat(T_bot) + dqsatdT(T_bot)
                qs[ij] = qsat(ph[kstart], T_bot[ij]);
                dqsdT[ij] = dqsatdT(ph[kstart], T_bot[ij]);
            }
    }

}


template<typename TF>
Thermo_moist<TF>::Thermo_moist(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Thermo<TF>(masterin, gridin, fieldsin, inputin),
    boundary_cyclic(masterin, gridin),
    field3d_operators(master, grid, fieldsin),
    field3d_io(master, grid)
{
    auto& gd = grid.get_grid_data();

    swthermo = "moist";

    // 4th order code is not implemented in thermo_moist
    if (grid.get_spatial_order() == Grid_order::Fourth)
        throw std::runtime_error("swthermo=moist is not supported for swspatialorder=4\n");

    // Initialize the prognostic fields
    const std::string group_name = "thermo";

    fields.init_prognostic_field("thl", "Liquid water potential temperature", "K", group_name, gd.sloc);
    fields.init_prognostic_field("qt", "Total water mixing ratio", "kg kg-1", group_name, gd.sloc);

    // Get the diffusivities of temperature and moisture
    fields.sp.at("thl")->visc = inputin.get_item<TF>("fields", "svisc", "thl");
    fields.sp.at("qt")->visc = inputin.get_item<TF>("fields", "svisc", "qt");

    // Test if the diffusivities of theta and qt are equal, else throw error
    if (fields.sp.at("thl")->visc != fields.sp.at("qt")->visc)
        throw std::runtime_error("The diffusivities of temperature and moisture must be equal\n");

    bs.pbot = inputin.get_item<TF>("thermo", "pbot", "");

    // Get base state option (boussinesq or anelastic)
    std::string swbasestate_in = inputin.get_item<std::string>("thermo", "swbasestate", "", "");
    if (swbasestate_in == "boussinesq")
        bs.swbasestate = Basestate_type::boussinesq;
    else if (swbasestate_in == "anelastic")
        bs.swbasestate = Basestate_type::anelastic;
    else
        throw std::runtime_error("Invalid option for \"swbasestate\"");


    // BvS test for updating hydrostatic prssure during run
    // swupdate..=0 -> initial base state pressure used in saturation calculation
    // swupdate..=1 -> base state pressure updated before saturation calculation
    bs.swupdatebasestate = inputin.get_item<bool>("thermo", "swupdatebasestate", "", true);

    // Time variable surface pressure
    tdep_pbot = std::make_unique<Timedep<TF>>(master, grid, "p_sbot", inputin.get_item<bool>("thermo", "swtimedep_pbot", "", false));

    available_masks.insert(available_masks.end(), {"ql", "qlcore", "bplus", "bmin"});
}

template<typename TF>
Thermo_moist<TF>::~Thermo_moist()
{
}

template<typename TF>
void Thermo_moist<TF>::init()
{
    auto& gd = grid.get_grid_data();

    bs.thl0.resize(gd.kcells);
    bs.qt0.resize(gd.kcells);
    bs.thvref.resize(gd.kcells);
    bs.thvrefh.resize(gd.kcells);
    bs.exnref.resize(gd.kcells);
    bs.exnrefh.resize(gd.kcells);
    bs.pref.resize(gd.kcells);
    bs.prefh.resize(gd.kcells);
    bs.rhoref.resize(gd.kcells);
    bs.rhorefh.resize(gd.kcells);
}

template<typename TF>
void Thermo_moist<TF>::save(const int iotime)
{
    auto& gd = grid.get_grid_data();

    int nerror = 0;

    if ( (master.get_mpiid() == 0) && bs.swupdatebasestate)
    {
        // Save the base state to disk
        FILE *pFile;
        char filename[256];
        std::sprintf(filename, "%s.%07d", "thermo_basestate", iotime);
        pFile = fopen(filename, "wbx");
        master.print_message("Saving \"%s\" ... ", filename);

        if (pFile == NULL)
        {
            master.print_message("FAILED\n");
            nerror++;
        }
        else
            master.print_message("OK\n");

        fwrite(&bs.thvref [gd.kstart], sizeof(TF), gd.ktot  , pFile);
        fwrite(&bs.thvrefh[gd.kstart], sizeof(TF), gd.ktot+1, pFile);
        fclose(pFile);
    }

    auto tmp1 = fields.get_tmp();

    // Save surface values thl + qt, which are needed for bitwise identical restarts
    auto save_2d_field = [&](
            TF* const restrict field, const std::string& name)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", name.c_str(), iotime);
        master.print_message("Saving \"%s\" ... ", filename);

        const int kslice = 0;
        if (field3d_io.save_xy_slice(
                field, tmp1->fld.data(), filename, kslice))
        {
            master.print_message("FAILED\n");
            nerror += 1;
        }
        else
            master.print_message("OK\n");
    };

    save_2d_field(fields.sp.at("thl")->fld_bot.data(), "thl_bot");
    save_2d_field(fields.sp.at("qt")->fld_bot.data(), "qt_bot");

    fields.release_tmp(tmp1);

    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("Error in writing thermo_basestate");
}

template<typename TF>
void Thermo_moist<TF>::load(const int iotime)
{
    auto& gd = grid.get_grid_data();

    int nerror = 0;

    if ( (master.get_mpiid() == 0) && bs.swupdatebasestate)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", "thermo_basestate", iotime);

        std::printf("Loading \"%s\" ... ", filename);

        FILE* pFile;
        pFile = fopen(filename, "rb");
        if (pFile == NULL)
        {
            master.print_message("FAILED\n");
            ++nerror;
        }
        else
        {
            master.print_message("OK\n");

            fread(&bs.thvref [gd.kstart], sizeof(TF), gd.ktot  , pFile);
            fread(&bs.thvrefh[gd.kstart], sizeof(TF), gd.ktot+1, pFile);
            fclose(pFile);
        }
    }

    auto tmp1 = fields.get_tmp();

    // Lambda function to load 2D fields.
    auto load_2d_field = [&](
            TF* const restrict field, const std::string& name)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", name.c_str(), iotime);
        master.print_message("Loading \"%s\" ... ", filename);

        if (field3d_io.load_xy_slice(
                field, tmp1->fld.data(),
                filename))
        {
            master.print_message("FAILED\n");
            nerror += 1;
        }
        else
            master.print_message("OK\n");

        boundary_cyclic.exec_2d(field);
    };

    load_2d_field(fields.sp.at("thl")->fld_bot.data(), "thl_bot");
    load_2d_field(fields.sp.at("qt")->fld_bot.data(), "qt_bot");

    fields.release_tmp(tmp1);

    // Communicate the file read error over all procs.
    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("Error in loading thermo_moist basestate");

    master.broadcast(&bs.thvref [gd.kstart], gd.ktot  );
    master.broadcast(&bs.thvrefh[gd.kstart], gd.ktot+1);
}

template<typename TF>
void Thermo_moist<TF>::create_basestate(Input& inputin, Netcdf_handle& input_nc)
{
    auto& gd = grid.get_grid_data();

    // Enable automated calculation of horizontally averaged fields
    fields.set_calc_mean_profs(true);

    // Calculate the base state profiles. With swupdatebasestate=1, these profiles are updated on every iteration.
    // 1. Take the initial profile as the reference
    const std::vector<int> start = {0};
    const std::vector<int> count = {gd.ktot};

    Netcdf_group& group_nc = input_nc.get_group("init");
    group_nc.get_variable(bs.thl0, "thl", start, count);
    group_nc.get_variable(bs.qt0, "qt", start, count);

    // Shift the vector
    std::rotate(bs.thl0.rbegin(), bs.thl0.rbegin() + gd.kstart, bs.thl0.rend());
    std::rotate(bs.qt0.rbegin(), bs.qt0.rbegin() + gd.kstart, bs.qt0.rend());

    calc_top_and_bot(
            bs.thl0.data(), bs.qt0.data(), gd.z.data(), gd.zh.data(), gd.dzhi.data(), gd.kstart, gd.kend);

    // 4. Calculate the initial/reference base state
    calc_base_state(
            bs.pref.data(), bs.prefh.data(), bs.rhoref.data(), bs.rhorefh.data(), bs.thvref.data(),
            bs.thvrefh.data(), bs.exnref.data(), bs.exnrefh.data(), bs.thl0.data(), bs.qt0.data(), bs.pbot,
            gd.kstart, gd.kend, gd.z.data(), gd.dz.data(), gd.dzh.data());

    // 5. In Boussinesq mode, overwrite reference temperature and density
    if (bs.swbasestate == Basestate_type::boussinesq)
    {
        bs.thvref0 = inputin.get_item<TF>("thermo", "thvref0", "");

        for (int k=0; k<gd.kcells; ++k)
        {
            bs.rhoref[k]  = 1.;
            bs.rhorefh[k] = 1.;
            bs.thvref[k]  = bs.thvref0;
            bs.thvrefh[k] = bs.thvref0;
        }
    }

    // 6. Copy the initial reference to the fields. This is the reference used in the dynamics.
    //    This one is not updated throughout the simulation to be consistent with the anelastic approximation.
    fields.rhoref = bs.rhoref;
    fields.rhorefh = bs.rhorefh;
}

template<typename TF>
void Thermo_moist<TF>::create(Input& inputin, Netcdf_handle& input_nc, Stats<TF>& stats, Column<TF>& column, Cross<TF>& cross, Dump<TF>& dump)
{
    create_basestate(inputin, input_nc);

    // 7. Process the time dependent surface pressure
    std::string timedep_dim = "time_surface";
    tdep_pbot->create_timedep(input_nc, timedep_dim);

    // Init the toolbox classes.
    boundary_cyclic.init();

    // Set up output classes
    create_stats(stats);
    create_column(column);
    create_dump(dump);
    create_cross(cross);
}

#ifndef USECUDA
template<typename TF>
void Thermo_moist<TF>::exec(const double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // Re-calculate hydrostatic pressure and exner, pass dummy as thvref to prevent overwriting base state
    auto tmp = fields.get_tmp();
    if (bs.swupdatebasestate)
    {
        calc_base_state(
                bs.pref.data(), bs.prefh.data(),
                bs.rhoref.data(), bs.rhorefh.data(),
                bs.thvref.data(), bs.thvrefh.data(),
                bs.exnref.data(), bs.exnrefh.data(),
                fields.sp.at("thl")->fld_mean.data(),
                fields.sp.at("qt")->fld_mean.data(),
                bs.pbot, gd.kstart, gd.kend,
                gd.z.data(), gd.dz.data(), gd.dzh.data());
    }

    // extend later for gravity vector not normal to surface
    calc_buoyancy_tend_2nd(
            fields.mt.at("w")->fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), bs.prefh.data(),
            &tmp->fld[0*gd.ijcells], &tmp->fld[1*gd.ijcells],
            &tmp->fld[2*gd.ijcells], &tmp->fld[3*gd.ijcells], bs.thvrefh.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    fields.release_tmp(tmp);

    stats.calc_tend(*fields.mt.at("w"), tend_name);
}
#endif

template<typename TF>
unsigned long Thermo_moist<TF>::get_time_limit(unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}

//#ifndef USECUDA
template<typename TF>
void Thermo_moist<TF>::get_mask(Stats<TF>& stats, std::string mask_name)
{

    if (mask_name == "ql")
    {
        auto ql = fields.get_tmp();
        auto qlh = fields.get_tmp();

        get_thermo_field(*ql, "ql", true, false);
        get_thermo_field(*qlh, "ql_h", true, false);

        stats.set_mask_thres(mask_name, *ql, *qlh, 0., Stats_mask_type::Plus);

        fields.release_tmp(ql);
        fields.release_tmp(qlh);
    }
    else if (mask_name == "qlcore")
    {
        auto ql = fields.get_tmp();
        auto qlh = fields.get_tmp();

        get_thermo_field(*ql, "ql", true, false);
        get_thermo_field(*qlh, "ql_h", true, false);

        stats.set_mask_thres(mask_name, *ql, *qlh, 0., Stats_mask_type::Plus);

        fields.release_tmp(ql);
        fields.release_tmp(qlh);

        auto b = fields.get_tmp();
        auto bh = fields.get_tmp();

        get_thermo_field(*b, "b", true, true);
        get_thermo_field(*bh, "b_h", true, true);

        field3d_operators.calc_mean_profile(b->fld_mean.data(), b->fld.data());
        field3d_operators.calc_mean_profile(bh->fld_mean.data(), bh->fld.data());
        field3d_operators.subtract_mean_profile(b->fld.data(), b->fld_mean.data());
        field3d_operators.subtract_mean_profile(bh->fld.data(), bh->fld_mean.data());

        stats.set_mask_thres(mask_name, *b, *bh, 0., Stats_mask_type::Plus);

        fields.release_tmp(b);
        fields.release_tmp(bh);
    }
    else if (mask_name == "bplus" || mask_name == "bmin")
    {
        auto b = fields.get_tmp();
        auto bh = fields.get_tmp();

        get_thermo_field(*b, "b", true, true);
        get_thermo_field(*bh, "b_h", true, true);

        field3d_operators.calc_mean_profile(b->fld_mean.data(), b->fld.data());
        field3d_operators.calc_mean_profile(bh->fld_mean.data(), bh->fld.data());
        field3d_operators.subtract_mean_profile(b->fld.data(), b->fld_mean.data());
        field3d_operators.subtract_mean_profile(bh->fld.data(), bh->fld_mean.data());

        if (mask_name == "bplus")
            stats.set_mask_thres(mask_name, *b, *bh, 0., Stats_mask_type::Plus);
        else
            stats.set_mask_thres(mask_name, *b, *bh, 0., Stats_mask_type::Min);

        fields.release_tmp(b);
        fields.release_tmp(bh);
    }
    else
    {
        std::string message = "Moist thermodynamics can not provide mask: \"" + mask_name +"\"";
        throw std::runtime_error(message);
    }
}
//#endif

template<typename TF>
bool Thermo_moist<TF>::has_mask(std::string mask_name)
{
    if (std::find(available_masks.begin(), available_masks.end(), mask_name) != available_masks.end())
        return true;
    else
        return false;
}

template<typename TF>
bool Thermo_moist<TF>::check_field_exists(const std::string name)
{
    if (name == "b" || name == "ql" || name == "T" || name == "qi")
        return true;
    else
        return false;
}

template<typename TF>
void Thermo_moist<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    tdep_pbot->update_time_dependent(bs.pbot, timeloop);
}

template<typename TF>
void Thermo_moist<TF>::get_thermo_field(
        Field3d<TF>& fld, const std::string& name, const bool cyclic, const bool is_stat)
{
    auto& gd = grid.get_grid_data();

    Background_state base;
    if (is_stat)
        base = bs_stats;
    else
        base = bs;

    // BvS: get_thermo_field() is called from subgrid-model, before thermo(), so re-calculate the hydrostatic pressure
    // Pass dummy as rhoref,bs.thvref to prevent overwriting base state
    if (bs.swupdatebasestate)
    {
        auto tmp = fields.get_tmp();
        calc_base_state(base.pref.data(), base.prefh.data(), &tmp->fld[0*gd.kcells], &tmp->fld[1*gd.kcells], &tmp->fld[2*gd.kcells],
                        &tmp->fld[3*gd.kcells], base.exnref.data(), base.exnrefh.data(), fields.sp.at("thl")->fld_mean.data(),
                        fields.sp.at("qt")->fld_mean.data(), base.pbot, gd.kstart, gd.kend, gd.z.data(), gd.dz.data(), gd.dzh.data());
        fields.release_tmp(tmp);
    }

    if (name == "b")
    {
        auto tmp  = fields.get_tmp();
        auto tmp2 = fields.get_tmp();
        calc_buoyancy(
                fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.pref.data(),
                tmp->fld.data(), tmp2->fld.data(), base.thvref.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.kcells, gd.icells, gd.ijcells);
        fields.release_tmp(tmp );
        fields.release_tmp(tmp2);
    }
    else if (name == "b_h")
    {
        auto tmp = fields.get_tmp();
        calc_buoyancy_h(
                fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.prefh.data(), base.thvrefh.data(),
                &tmp->fld[0*gd.ijcells], &tmp->fld[1*gd.ijcells], &tmp->fld[2*gd.ijcells], &tmp->fld[3*gd.ijcells],
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        fields.release_tmp(tmp);
    }
    else if (name == "ql")
    {
        calc_liquid_water(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.pref.data(),
                          gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
    }
    else if (name == "ql_h")
    {
        auto tmp = fields.get_tmp();
        calc_liquid_water_h(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.prefh.data(), &tmp->fld[0*gd.ijcells], &tmp->fld[1*gd.ijcells],
                            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        fields.release_tmp(tmp);
    }
    else if (name == "qi")
    {
        calc_ice(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.pref.data(),
                 gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
    }
    else if (name == "ql_qi")
    {
        calc_condensate(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.pref.data(),
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
    }
    else if (name == "qsat")
    {
        calc_saturated_water_vapor(
                fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.pref.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
    }
    else if (name == "rh")
    {
        calc_relative_humidity(
                fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.pref.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
    }
    else if (name == "N2")
    {
        calc_N2(fld.fld.data(), fields.sp.at("thl")->fld.data(), gd.dzi.data(), base.thvref.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
    }
    else if (name == "T")
    {
        calc_T(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.pref.data(), base.exnref.data(),
               gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
    }
    else if (name == "T_h")
    {
        auto tmp = fields.get_tmp();
        calc_T_h(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.prefh.data(), &tmp->fld[0*gd.ijcells], &tmp->fld[1*gd.ijcells],
                 &tmp->fld[2*gd.ijcells], gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        fields.release_tmp(tmp);
    }
    else if (name == "thv")
    {
        // Calculate thv including one ghost cell, for the calculation of gradients in the statistics
        calc_thv(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.pref.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart-1, gd.kend+1, gd.icells, gd.ijcells);
    }
    else if (name == "thv_fluxbot")
    {
        calc_thv_fluxbot(
                fld.flux_bot.data(),
                fields.sp.at("thl")->fld.data(), fields.sp.at("thl")->flux_bot.data(),
                fields.sp.at("qt") ->fld.data(), fields.sp.at("qt") ->flux_bot.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.icells, gd.ijcells);
    }
    else
    {
        std::string error_message = "Can not get thermo field: \"" + name + "\"";
        throw std::runtime_error(error_message);
    }

    if (cyclic)
        boundary_cyclic.exec(fld.fld.data());
}

template<typename TF>
void Thermo_moist<TF>::get_radiation_fields(
        Field3d<TF>& T, Field3d<TF>& T_h, Field3d<TF>& qv, Field3d<TF>& clwp, Field3d<TF>& ciwp) const
{
    auto& gd = grid.get_grid_data();

    calc_radiation_fields(
            T.fld.data(), T_h.fld.data(), qv.fld.data(),
            clwp.fld.data(), ciwp.fld.data(), T_h.fld_bot.data(),
            T.fld_bot.data(), T.fld_top.data(),  // These 2d fields are used as tmp arrays.
            fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(),
            fields.sp.at("thl")->fld_bot.data(),
            bs.pref.data(), bs.prefh.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.igc, gd.jgc, gd.kgc,
            gd.icells, gd.ijcells,
            gd.imax, gd.imax*gd.jmax);
}

template<typename TF>
void Thermo_moist<TF>::get_radiation_columns(
    Field3d<TF>& tmp, std::vector<int>& col_i, std::vector<int>& col_j) const
{
    auto& gd = grid.get_grid_data();

    // Get slices from tmp field.
    const int n_cols = col_i.size();
    int offset = 0;
    TF* t_lay_a = &tmp.fld.data()[offset]; offset += n_cols * gd.ktot;
    TF* t_lev_a = &tmp.fld.data()[offset]; offset += n_cols * (gd.ktot+1);
    TF* t_sfc_a = &tmp.fld.data()[offset]; offset += n_cols;
    TF* h2o_a   = &tmp.fld.data()[offset]; offset += n_cols * (gd.ktot);
    TF* clwp_a  = &tmp.fld.data()[offset]; offset += n_cols * (gd.ktot);
    TF* ciwp_a  = &tmp.fld.data()[offset];

    calc_radiation_columns(
            t_lay_a, t_lev_a, h2o_a, clwp_a, ciwp_a, t_sfc_a,
            fields.sp.at("thl")->fld.data(),
            fields.sp.at("qt")->fld.data(),
            fields.sp.at("thl")->fld_bot.data(),
            bs.pref.data(),
            bs.prefh.data(),
            col_i.data(),
            col_j.data(),
            n_cols,
            gd.kgc, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);
}

template<typename TF>
void Thermo_moist<TF>::get_land_surface_fields(
        std::vector<TF>& T_bot, std::vector<TF>& T_a, std::vector<TF>& vpd,
        std::vector<TF>& qsat, std::vector<TF>& dqsatdT)
{
    /* Calculate the thermo fields required by the LSM in
       2D slices in the 3D tmp field */
    auto& gd = grid.get_grid_data();

    calc_land_surface_fields(
            T_bot.data(),
            T_a.data(),
            vpd.data(),
            qsat.data(),
            dqsatdT.data(),
            fields.sp.at("thl")->fld_bot.data(),
            fields.sp.at("thl")->fld.data(),
            fields.sp.at("qt")->fld.data(),
            bs.exnref.data(),
            bs.exnrefh.data(),
            bs.pref.data(),
            bs.prefh.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart,
            gd.icells, gd.ijcells);
}

template<typename TF>
void Thermo_moist<TF>::get_buoyancy_surf(
        std::vector<TF>& b, std::vector<TF>& bbot, bool is_stat)
{
    auto& gd = grid.get_grid_data();

    Background_state base;
    if (is_stat)
        base = bs_stats;
    else
        base = bs;

    calc_buoyancy_bot(
            b.data(), bbot.data(),
            fields.sp.at("thl")->fld.data(),
            fields.sp.at("thl")->fld_bot.data(),
            fields.sp.at("qt")->fld.data(),
            fields.sp.at("qt")->fld_bot.data(),
            base.thvref.data(),
            base.thvrefh.data(),
            gd.icells, gd.jcells,
            gd.ijcells, gd.kstart);

    //calc_buoyancy_fluxbot(
    //        b.flux_bot.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("thl")->flux_bot.data(),
    //        fields.sp.at("qt")->fld.data(), fields.sp.at("qt")->flux_bot.data(), base.thvrefh.data(),
    //        gd.icells, gd.jcells, gd.kstart, gd.ijcells);
}

template<typename TF>
void Thermo_moist<TF>::get_buoyancy_surf(
        std::vector<TF>& b_bot, std::vector<TF>& thl_bot, std::vector<TF>& qt_bot)
{
    auto& gd = grid.get_grid_data();

    calc_buoyancy_bot(
            b_bot.data(), thl_bot.data(), qt_bot.data(),
            bs.thvrefh.data(), gd.icells, gd.jcells, gd.kstart);
}

template<typename TF>
void Thermo_moist<TF>::get_buoyancy_fluxbot(
        std::vector<TF>& bfluxbot, bool is_stat)
{
    auto& gd = grid.get_grid_data();

    Background_state base;
    if (is_stat)
        base = bs_stats;
    else
        base = bs;

    calc_buoyancy_fluxbot(
            bfluxbot.data(),
            fields.sp.at("thl")->fld.data(),
            fields.sp.at("thl")->flux_bot.data(),
            fields.sp.at("qt")->fld.data(),
            fields.sp.at("qt")->flux_bot.data(),
            base.thvrefh.data(),
            gd.icells, gd.jcells,
            gd.kstart, gd.ijcells);
}

template<typename TF>
void Thermo_moist<TF>::get_temperature_bot(Field3d<TF>& T_bot, bool is_stat)
{
    auto& gd = grid.get_grid_data();
    Background_state base;
    if (is_stat)
        base = bs_stats;
    else
        base = bs;

    calc_T_bot(T_bot.fld_bot.data(), fields.sp.at("thl")->fld.data(), base.exnrefh.data(), base.thl0.data(),
               gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.icells, gd.ijcells);
}

template<typename TF>
const std::vector<TF>& Thermo_moist<TF>::get_basestate_vector(std::string name) const
{
    if (name == "p")
        return bs.pref;
    else if (name == "ph")
        return bs.prefh;
    else if (name == "exner")
        return bs.exnref;
    else if (name == "exnerh")
        return bs.exnrefh;
    else if (name == "rho")
        return bs.rhoref;
    else if (name == "rhoh")
        return bs.rhorefh;
    else if (name == "thv")
        return bs.thvref;
    else if (name == "thvh")
        return bs.thvrefh;
    else
    {
        std::string error = "Thermo_moist::get_basestate_vector() can't return \"" + name + "\"";
        throw std::runtime_error(error);
    }
}

template<typename TF>
TF Thermo_moist<TF>::get_db_ref() const
{
    auto& gd = grid.get_grid_data();
    return Constants::grav<TF>/bs.thvref[gd.kstart]*(bs.thvref[gd.kstart] - bs.thvrefh[gd.kstart]);
}



template<typename TF>
void Thermo_moist<TF>::get_prog_vars(std::vector<std::string>& list)
{
    list.push_back("thl");
    list.push_back("qt");
}

template<typename TF>
TF Thermo_moist<TF>::get_buoyancy_diffusivity()
{
    // Use the diffusivity from the liquid water potential temperature
    return fields.sp.at("thl")->visc;
}

template<typename TF>
int Thermo_moist<TF>::get_bl_depth()
{
    // Use the liquid water potential temperature gradient to find the BL depth
    auto& gd = grid.get_grid_data();
    return calc_zi(fields.sp.at("thl")->fld_mean.data(), gd.kstart, gd.kend, 1);
}

template<typename TF>
void Thermo_moist<TF>::create_stats(Stats<TF>& stats)
{
    const std::string group_name = "thermo";

    bs_stats = bs;

    // Add variables to the statistics
    if (stats.get_switch())
    {
        /* Add fixed base-state density and temperature profiles. Density should probably be in fields (?), but
           there the statistics are initialized before thermo->create() is called */
        stats.add_fixed_prof("rhoref",  "Full level basic state density", "kg m-3", "z" , group_name, bs.rhoref );
        stats.add_fixed_prof("rhorefh", "Half level basic state density", "kg m-3", "zh", group_name, bs.rhorefh);
        stats.add_fixed_prof("thvref", "Full level basic state virtual potential temperature", "K", "z" , group_name, bs.thvref);
        stats.add_fixed_prof("thvrefh", "Half level basic state virtual potential temperature", "K", "zh", group_name, bs.thvrefh);

        if (bs_stats.swupdatebasestate)
        {
            stats.add_prof("phydro", "Full level hydrostatic pressure", "Pa", "z" , group_name);
            stats.add_prof("phydroh", "Half level hydrostatic pressure", "Pa", "zh", group_name);
            stats.add_prof("rho",  "Full level density", "kg m-3", "z" , group_name);
            stats.add_prof("rhoh", "Half level density", "kg m-3", "zh", group_name);
        }
        else
        {
            stats.add_fixed_prof("pydroh",  "Full level hydrostatic pressure", "Pa", "z" , group_name, bs.pref);
            stats.add_fixed_prof("phydroh", "Half level hydrostatic pressure", "Pa", "zh", group_name, bs.prefh);
        }

        auto thv = fields.get_tmp();
        thv->name = "thv";
        thv->longname = "Virtual potential temperature";
        thv->unit = "K";
        stats.add_profs(*thv, "z", {"mean", "2", "w", "grad", "diff", "flux"}, group_name);
        fields.release_tmp(thv);

        auto T = fields.get_tmp();
        T->name = "T";
        T->longname = "Absolute temperature";
        T->unit = "K";
        stats.add_profs(*T, "z", {"mean", "2"}, group_name);
        fields.release_tmp(T);

        auto ql = fields.get_tmp();
        ql->name = "ql";
        ql->longname = "Liquid water";
        ql->unit = "kg kg-1";
        stats.add_profs(*ql, "z", {"mean", "frac", "path", "cover"}, group_name);
        fields.release_tmp(ql);

        auto qi = fields.get_tmp();
        qi->name = "qi";
        qi->longname = "Ice";
        qi->unit = "kg kg-1";
        stats.add_profs(*qi, "z", {"mean", "frac", "path", "cover"}, group_name);
        fields.release_tmp(qi);

        auto qlqi = fields.get_tmp();
        qlqi->name = "qlqi";
        qlqi->longname = "Liquid water and ice";
        qlqi->unit = "kg kg-1";
        stats.add_profs(*qlqi, "z", {"mean", "frac", "path", "cover"}, group_name);
        fields.release_tmp(qlqi);

        auto qsat = fields.get_tmp();
        qsat->name = "qsat";
        qsat->longname = "Saturated water vapor";
        qsat->unit = "kg kg-1";
        stats.add_profs(*qsat, "z", {"mean", "path"}, group_name);
        fields.release_tmp(qsat);

        auto rh = fields.get_tmp();
        rh->name = "rh";
        rh->longname = "Relative humidity";
        rh->unit = "-";
        stats.add_profs(*rh, "z", {"mean"}, group_name);
        fields.release_tmp(rh);

        stats.add_time_series("zi", "Boundary Layer Depth", "m", group_name);

        stats.add_time_series("thl_bot", "Surface liquid water potential temperature", "K", group_name);
        stats.add_time_series("qt_bot", "Surface specific humidity", "kg kg-1", group_name);

        stats.add_tendency(*fields.mt.at("w"), "zh", tend_name, tend_longname, group_name);
    }
}

template<typename TF>
void Thermo_moist<TF>::create_column(Column<TF>& column)
{
    // add the profiles to the columns
    if (column.get_switch())
    {
        // Vertical profiles
        column.add_prof("thv", "Virtual potential temperature", "K", "z");
        column.add_prof("ql", "Liquid water mixing ratio", "kg kg-1", "z");
        column.add_prof("qi", "Ice mixing ratio", "kg kg-1", "z");

        column.add_time_series("ql_path", "Liquid water path", "kg m-2");
        column.add_time_series("qi_path", "Ice path", "kg m-2");
    }
}

template<typename TF>
void Thermo_moist<TF>::create_cross(Cross<TF>& cross)
{
    if (cross.get_switch())
    {
        swcross_b = false;
        swcross_ql = false;
        swcross_qi = false;
        swcross_qlqi = false;
        swcross_qsat = false;
        swcross_qlqithv = false;

        // Vectors with allowed cross variables for buoyancy and liquid water.
        const std::vector<std::string> allowed_crossvars_b = {"b", "bbot", "bfluxbot"};
        const std::vector<std::string> allowed_crossvars_ql = {"ql", "qlpath", "qlbase", "qltop"};
        const std::vector<std::string> allowed_crossvars_qi = {"qi", "qipath"};
        const std::vector<std::string> allowed_crossvars_qlqi = {"qlqipath", "qlqibase", "qlqitop"};
        const std::vector<std::string> allowed_crossvars_qsat = {"qsatpath"};
        const std::vector<std::string> allowed_crossvars_misc = {"w500hpa"};
        const std::vector<std::string> allowed_crossvars_qlqithv = {"qlqicore_max_thv_prime"};

        std::vector<std::string> bvars  = cross.get_enabled_variables(allowed_crossvars_b);
        std::vector<std::string> qlvars = cross.get_enabled_variables(allowed_crossvars_ql);
        std::vector<std::string> qivars = cross.get_enabled_variables(allowed_crossvars_qi);
        std::vector<std::string> qlqivars = cross.get_enabled_variables(allowed_crossvars_qlqi);
        std::vector<std::string> qsatvars = cross.get_enabled_variables(allowed_crossvars_qsat);
        std::vector<std::string> miscvars = cross.get_enabled_variables(allowed_crossvars_misc);
        std::vector<std::string> qlqithvvars = cross.get_enabled_variables(allowed_crossvars_qlqithv);

        if (bvars.size() > 0)
            swcross_b  = true;

        if (qlvars.size() > 0)
            swcross_ql = true;

        if (qivars.size() > 0)
            swcross_qi = true;

        if (qlqivars.size() > 0)
            swcross_qlqi = true;

        if (qsatvars.size() > 0)
            swcross_qsat = true;

        if (qlqithvvars.size() > 0)
            swcross_qlqithv = true;

        // Merge into one vector
        crosslist = bvars;
        crosslist.insert(crosslist.end(), qlvars.begin(), qlvars.end());
        crosslist.insert(crosslist.end(), qivars.begin(), qivars.end());
        crosslist.insert(crosslist.end(), qlqivars.begin(), qlqivars.end());
        crosslist.insert(crosslist.end(), qsatvars.begin(), qsatvars.end());
        crosslist.insert(crosslist.end(), miscvars.begin(), miscvars.end());
        crosslist.insert(crosslist.end(), qlqithvvars.begin(), qlqithvvars.end());
    }
}

template<typename TF>
void Thermo_moist<TF>::create_dump(Dump<TF>& dump)
{
    if (dump.get_switch())
    {
        // Get global cross-list from cross.cxx
        std::vector<std::string>& dumplist_global = dump.get_dumplist();

        // Check if fields in dumplist are retrievable thermo fields
        std::vector<std::string>::iterator dumpvar = dumplist_global.begin();
        while (dumpvar != dumplist_global.end())
        {
            if (check_field_exists(*dumpvar))
            {
                // Remove variable from global list, put in local list
                dumplist.push_back(*dumpvar);
                dumplist_global.erase(dumpvar); // erase() returns iterator of next element..
            }
            else
                ++dumpvar;
        }
    }
}

template<typename TF>
void Thermo_moist<TF>::exec_stats(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    #ifndef USECUDA
    bs_stats = bs;
    #endif

    const TF no_offset = 0.;
    const TF no_threshold = 0.;

    // Calculate the virtual temperature stats.
    auto thv = fields.get_tmp();
    thv->loc = gd.sloc;
    get_thermo_field(*thv, "thv", true, true);
    get_thermo_field(*thv, "thv_fluxbot", true, true);

    stats.calc_stats("thv", *thv, no_offset, no_threshold);

    fields.release_tmp(thv);

    // Calculate the absolute temperature stats.
    auto T = fields.get_tmp();
    T->loc = gd.sloc;

    get_thermo_field(*T, "T", true, true);
    stats.calc_stats("T", *T, no_offset, no_threshold);

    fields.release_tmp(T);

    // Calculate the liquid water stats
    auto ql = fields.get_tmp();
    ql->loc = gd.sloc;

    get_thermo_field(*ql, "ql", true, true);
    stats.calc_stats("ql", *ql, no_offset, no_threshold);

    fields.release_tmp(ql);

    // Calculate the ice stats
    auto qi = fields.get_tmp();
    qi->loc = gd.sloc;

    get_thermo_field(*qi, "qi", true, true);
    stats.calc_stats("qi", *qi, no_offset, no_threshold);

    fields.release_tmp(qi);

    // Calculate the combined liquid water and ice stats
    auto qlqi = fields.get_tmp();
    qlqi->loc = gd.sloc;

    get_thermo_field(*qlqi, "ql_qi", true, true);
    stats.calc_stats("qlqi", *qlqi, no_offset, no_threshold);

    fields.release_tmp(qlqi);

    // Calculate the saturated water vapor stats
    auto qsat = fields.get_tmp();
    qsat->loc = gd.sloc;

    get_thermo_field(*qsat, "qsat", true, true);
    stats.calc_stats("qsat", *qsat, no_offset, no_threshold);

    fields.release_tmp(qsat);

    // Calculate the relative humidity
    auto rh = fields.get_tmp();
    rh->loc = gd.sloc;

    get_thermo_field(*rh, "rh", true, true);
    stats.calc_stats("rh", *rh, no_offset, no_threshold);

    // Surface values
    stats.calc_stats_2d("thl_bot", fields.ap.at("thl")->fld_bot, no_offset);
    stats.calc_stats_2d("qt_bot", fields.ap.at("qt")->fld_bot, no_offset);

    fields.release_tmp(rh);

    if (bs_stats.swupdatebasestate)
    {
        stats.set_prof("phydro" , bs_stats.pref);
        stats.set_prof("phydroh", bs_stats.prefh);
        stats.set_prof("rho"    , bs_stats.rhoref);
        stats.set_prof("rhoh"   , bs_stats.rhorefh);
    }

    stats.set_time_series("zi", gd.z[get_bl_depth()]);
}

#ifndef USECUDA
template<typename TF>
void Thermo_moist<TF>::exec_column(Column<TF>& column)
{
    auto& gd = grid.get_grid_data();

    #ifndef USECUDA
    bs_stats = bs;
    #endif

    const TF no_offset = 0.;
    auto output = fields.get_tmp();

    // Vertical profiles
    get_thermo_field(*output, "thv", false, true);
    column.calc_column("thv", output->fld.data(), no_offset);

    get_thermo_field(*output, "ql", false, true);

    calc_path(
        output->fld_bot.data(),
        output->fld.data(),
        bs_stats.rhoref.data(),
        gd.dz.data(),
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);

    column.calc_column("ql", output->fld.data(), no_offset);
    column.calc_time_series("ql_path", output->fld_bot.data(), no_offset);

    get_thermo_field(*output, "qi", false, true);

    calc_path(
        output->fld_bot.data(),
        output->fld.data(),
        bs_stats.rhoref.data(),
        gd.dz.data(),
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);

    column.calc_column("qi", output->fld.data(), no_offset);
    column.calc_time_series("qi_path", output->fld_bot.data(), no_offset);

    fields.release_tmp(output);
}
#endif


template<typename TF>
void Thermo_moist<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    auto& gd = grid.get_grid_data();

    #ifndef USECUDA
    bs_stats = bs;
    #endif

    auto output = fields.get_tmp();

    if (swcross_b)
    {
        get_thermo_field(*output, "b", false, true);
        get_buoyancy_fluxbot(output->flux_bot, true);
    }

    for (auto& it : crosslist)
    {
        if (it == "b")
            cross.cross_simple(output->fld.data(), "b", iotime, gd.sloc);
        else if (it == "blngrad")
            cross.cross_lngrad(output->fld.data(), "blngrad", iotime);
        else if (it == "bbot")
            cross.cross_plane(output->fld_bot.data(), "bbot", iotime);
        else if (it == "bfluxbot")
            cross.cross_plane(output->flux_bot.data(), "bfluxbot", iotime);
    }

    if (swcross_ql)
        get_thermo_field(*output, "ql", false, true);

    for (auto& it : crosslist)
    {
        if (it == "ql")
            cross.cross_simple(output->fld.data(), "ql", iotime, gd.sloc);
        if (it == "qlpath")
            cross.cross_path(output->fld.data(), "qlpath", iotime);
        if (it == "qlbase")
            cross.cross_height_threshold(output->fld.data(), 0., Cross_direction::Bottom_to_top, "qlbase", iotime);
        if (it == "qltop")
            cross.cross_height_threshold(output->fld.data(), 0., Cross_direction::Top_to_bottom, "qltop", iotime);
    }

    if (swcross_qi)
        get_thermo_field(*output, "qi", false, true);

    for (auto& it : crosslist)
    {
        if (it == "qi")
            cross.cross_simple(output->fld.data(), "qi", iotime, gd.sloc);
        if (it == "qipath")
            cross.cross_path(output->fld.data(), "qipath", iotime);
    }

    if (swcross_qlqi)
        get_thermo_field(*output, "ql_qi", false, true);

    for (auto& it : crosslist)
    {
        if (it == "qlqipath")
            cross.cross_path(output->fld.data(), "qlqipath", iotime);
        if (it == "qlqibase")
            cross.cross_height_threshold(output->fld.data(), 0., Cross_direction::Bottom_to_top, "qlqibase", iotime);
        if (it == "qlqitop")
            cross.cross_height_threshold(output->fld.data(), 0., Cross_direction::Top_to_bottom, "qlqitop", iotime);
    }

    if (swcross_qsat)
        get_thermo_field(*output, "qsat", false, true);

    for (auto& it : crosslist)
    {
        if (it == "qsatpath")
            cross.cross_path(output->fld.data(), "qsatpath", iotime);
    }

    for (auto& it : crosslist)
    {
        if (it == "w500hpa")
        {
            calc_w500hpa(
                    output->fld_bot.data(), fields.mp.at("w")->fld.data(), bs_stats.prefh.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

            cross.cross_plane(output->fld_bot.data(), "w500hpa", iotime);
        }
    }

    fields.release_tmp(output);

    if (swcross_qlqithv)
    {
        auto qlqi = fields.get_tmp();
        auto thv  = fields.get_tmp();

        get_thermo_field(*qlqi, "ql_qi", false, true);
        get_thermo_field(*thv,  "thv", false, true);

        field3d_operators.calc_mean_profile(thv->fld_mean.data(), thv->fld.data());

        for (auto& it : crosslist)
        {
            if (it == "qlqicore_max_thv_prime")
            {
                calc_qlqicore_max_thv_prime(
                        thv->fld_bot.data(),
                        qlqi->fld.data(), thv->fld.data(), thv->fld_mean.data(),
                        gd.istart, gd.iend,
                        gd.jstart, gd.jend,
                        gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);

                cross.cross_plane(thv->fld_bot.data(), "qlqicore_max_thv_prime", iotime);
            }
        }

        fields.release_tmp(qlqi);
        fields.release_tmp(thv);
    }
}

template<typename TF>
void Thermo_moist<TF>::exec_dump(Dump<TF>& dump, unsigned long iotime)
{
    #ifndef USECUDA
    bs_stats = bs;
    #endif

    auto output = fields.get_tmp();

    for (auto& it : dumplist)
    {
        if (check_field_exists(it))
            get_thermo_field(*output, it, false, true);
        else
        {
            std::string msg = "Thermo dump of field \"" + it + "\" not supported";
            throw std::runtime_error(msg);
        }

        dump.save_dump(output->fld.data(), it, iotime);
    }

    fields.release_tmp(output);
}

template class Thermo_moist<double>;
template class Thermo_moist<float>;
