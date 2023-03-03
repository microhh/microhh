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

#include "constants.h"

namespace Sb_common
{
    using namespace Constants;

    template<typename TF>
    void convert_unit(
            TF* const restrict a,
            const TF* const restrict rho,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride,
            bool to_kgm3)
    {
        for (int k=kstart; k<kend; k++)
        {
            const TF fac = (to_kgm3) ? rho[k] : TF(1.)/rho[k];

            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jstride + k*kstride;
                    a[ijk] *= fac;
                }
        }
    }


    template<typename TF>
    void copy_slice(
            TF* const restrict fld_2d,
            const TF* const restrict fld_3d,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        for (int j = jstart; j < jend; j++)
            #pragma ivdep
            for (int i = istart; i < iend; i++)
            {
                const int ij = i + j * jstride;
                const int ijk = i + j * jstride + k * kstride;

                fld_2d[ij] = fld_3d[ijk];
            }
    }


    template<typename TF>
    void copy_slice_and_integrate(
            TF* const restrict fld_2d,
            const TF* const restrict fld_3d,
            const TF* const restrict fld_3d_tend,
            const TF* const restrict rho,
            const TF dt,
            bool do_integration,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        const TF fac = do_integration ? 1 : 0;

        for (int j = jstart; j < jend; j++)
            #pragma ivdep
            for (int i = istart; i < iend; i++)
            {
                const int ij = i + j * jstride;
                const int ijk = i + j * jstride + k * kstride;

                // fld_3d_tend is still per kg, while fld_2d and fld_3d per m-3.
                fld_2d[ij] = fld_3d[ijk] + fac*dt*rho[k]*fld_3d_tend[ijk];
            }
    }


    template<typename TF>
    void implicit_core(
            TF* const restrict q_val,
            TF* const restrict q_sum,
            TF* const restrict q_impl,
            TF* const restrict vsed_new,
            TF* const restrict vsed_now,
            TF* const restrict flux_new,
            TF* const restrict flux_now,
            const TF rdzdt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        for (int j = jstart; j < jend; j++)
            #pragma ivdep
            for (int i = istart; i < iend; i++)
            {
                const int ij = i + j * jstride;

                // `new` on r.h.s. is `now` value from level above
                vsed_new[ij] = TF(0.5) * (vsed_now[ij] + vsed_new[ij]);

                // `flux_new` are the updated flux values from the level above
                // `flux_now` are here the old (current time step) flux values from the level above
                const TF flux_sum = flux_new[ij] + flux_now[ij];

                // `flux_now` are here overwritten with the current level
                flux_now[ij] = std::min(vsed_now[ij] * q_val[ij], flux_sum);   // loop dependency
                flux_now[ij] = std::max(flux_now[ij], TF(0));                  // Maybe not necessary

                // Time integrated value without implicit weight
                q_sum[ij] = q_val[ij] + rdzdt * (flux_sum - flux_now[ij]);

                // Implicit weight
                q_impl[ij] = TF(1) / (TF(1) + vsed_new[ij] * rdzdt);

                // prepare for source term calculation
                const TF q_star = q_impl[ij] * q_sum[ij];
                q_val[ij]  = q_star;       // source/sinks work on star-values
                q_sum[ij]  = q_sum[ij] - q_star;
            }
    }


    template<typename TF>
    void integrate_process(
            TF* const restrict val,
            const TF* const restrict tend,
            const double dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        for (int j = jstart; j < jend; j++)
            #pragma ivdep
            for (int i = istart; i < iend; i++)
            {
                const int ij = i + j * jstride;

                // Time integration
                val[ij] += tend[ij] * TF(dt);
            }
    }


    template<typename TF>
    void implicit_time(
            TF* const restrict q_val,
            TF* const restrict q_sum,
            TF* const restrict q_impl,
            TF* const restrict vsed_new,
            TF* const restrict vsed_now,
            TF* const restrict flux_new,
            const TF rdzdt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        for (int j = jstart; j < jend; j++)
            #pragma ivdep
            for (int i = istart; i < iend; i++)
            {
                const int ij = i + j * jstride;

                // Time integration
                q_val[ij] = std::max(TF(0), q_impl[ij] * (q_sum[ij] + q_val[ij]));

                // Prepare for next level
                flux_new[ij] = q_val[ij] * vsed_new[ij];
                vsed_new[ij] = vsed_now[ij];
            }
    }


    template<typename TF>
    void calc_thermo_tendencies_cloud(
            TF* const restrict thlt,
            TF* const restrict qtt,
            const TF* const restrict qr_tend_conversion,
            const TF* const restrict rho,
            const TF* const restrict exner,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        const TF rho_i = TF(1) / rho[k];

        for (int j = jstart; j < jend; j++)
            #pragma ivdep
            for (int i = istart; i < iend; i++)
            {
                const int ij  = i + j * jstride;
                const int ijk = i + j * jstride + k*kstride;

                // Tendencies `qt` and `thl` only include tendencies from microphysics conversions.
                qtt[ijk]  -= rho_i * qr_tend_conversion[ij];
                thlt[ijk] += rho_i * Lv<TF> / (cp<TF> * exner[k]) * qr_tend_conversion[ij];
            }
    }


    template<typename TF, bool sw_prognostic_ice>
    void calc_thermo_tendencies_cloud_ice(
            TF* const restrict thlt,
            TF* const restrict qtt,
            const TF* const restrict qv_to_ql,
            const TF* const restrict qv_to_qf,
            const TF* const restrict ql_to_qf,
            const TF* const restrict rho,
            const TF* const restrict exner,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        const TF rho_i = TF(1) / rho[k];

        for (int j = jstart; j < jend; j++)
                #pragma ivdep
                for (int i = istart; i < iend; i++)
                {
                    const int ij  = i + j * jstride;
                    const int ijk = i + j * jstride + k*kstride;

                    if (sw_prognostic_ice)
                        qtt[ijk] -= rho_i * qv_to_ql[ij];
                    else
                        qtt[ijk] -= rho_i * (qv_to_ql[ij] + qv_to_qf[ij]);

                    thlt[ijk] +=
                            ((rho_i * Lv<TF> / (cp<TF> * exner[k]) * qv_to_ql[ij]) +
                             (rho_i * Lf<TF> / (cp<TF> * exner[k]) * ql_to_qf[ij]) +
                             (rho_i * Ls<TF> / (cp<TF> * exner[k]) * qv_to_qf[ij]));
                }
    }

    template<typename TF>
    void diagnose_tendency(
            TF* const restrict tend,
            const TF* const restrict fld_old,
            const TF* const restrict fld_new,
            const TF* const restrict rho,
            const double dt,
            bool do_integration,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        const TF dt_i = TF(1) / dt;
        const TF rho_i = TF(1) / rho[k];
        const TF fac = do_integration ? 1 : 0;

        for (int j = jstart; j < jend; j++)
                #pragma ivdep
                for (int i = istart; i < iend; i++)
                {
                    const int ij = i + j * jstride;
                    const int ijk= i + j * jstride + k*kstride;

                    // Evaluate tendencies. This includes the tendencies from both conversions and implicit sedimentation.
                    // `Old` versions are integrated first with only the dynamics tendencies to avoid double counting.
                    tend[ijk] += rho_i * (fld_new[ij] - (fld_old[ijk] + fac*dt*rho[k]*tend[ijk])) * dt_i;
                }
    }
}