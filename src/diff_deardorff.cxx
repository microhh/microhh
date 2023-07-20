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

#include <algorithm>
#include <cmath>
#include <iostream>

#include "grid.h"
#include "fields.h"
#include "master.h"
#include "defines.h"
#include "constants.h"
#include "monin_obukhov.h"
#include "thermo.h"
#include "boundary.h"
#include "stats.h"
#include "fast_math.h"

#include "diff_deardorff.h"
#include "diff_kernels.h"

namespace
{
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;
    namespace dk = Diff_kernels;

    template <typename TF>
    void check_for_minval(
            TF* const restrict a,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells, const int ijcells,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int jj = icells;
        const int kk = ijcells;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    a[ijk] = std::max(a[ijk], Constants::sgstke_min<TF>);
                }

        boundary_cyclic.exec(a);
    }

    template <typename TF, Surface_model surface_model, bool sw_mason>
    void calc_evisc_neutral(
            TF* const restrict evisc,
            const TF* const restrict a,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* z0m,
            const TF dx, const TF dy,
            const TF cn, const TF cm,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells, const int ijcells,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int jj = icells;
        const int kk = ijcells;

        // Wall damping constant.
        constexpr TF n_mason = TF(2.);
        constexpr TF A_vandriest = TF(26.);

        if (surface_model == Surface_model::Disabled)
            throw std::runtime_error("Resolved wall not supported in Deardorff SGSm.");
        else
        {
            for (int k=kstart; k<kend; ++k) // Counter starts at kstart (as sgstke is defined here)
            {
                // Calculate geometric filter width, based on Deardorff (1980)
                const TF mlen0 = std::pow(dx*dy*dz[k], TF(1./3.));

                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ij = i + j*jj;
                        const int ijk = i + j*jj + k*kk;
                        TF fac;

                        if (sw_mason) // Apply Mason's wall correction
                            fac = std::pow(TF(1.)/(TF(1.)/std::pow(mlen0, n_mason) + TF(1.)/
                                        (std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);
                        else
                            fac = mlen0;

                        // Calculate eddy diffusivity for momentum.
                        evisc[ijk] = cm * fac * std::sqrt(a[ijk]);
                    }
            }
        }

        boundary_cyclic.exec(evisc);
    }

    template<typename TF, Surface_model surface_model, bool sw_mason>
    void calc_evisc(
            TF* const restrict evisc,
            const TF* const restrict a,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict N2,
            const TF* const restrict bgradbot,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* const restrict z0m,
            const TF dx, const TF dy,
            const TF cn, const TF cm,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells, const int ijcells,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int jj = icells;
        const int kk = ijcells;

       if (surface_model == Surface_model::Disabled)
            throw std::runtime_error("Resolved wall not supported in Deardorff SGSm.");
       else
       {
            // Variables for the wall damping and length scales
            const TF n_mason = TF(2.);
            TF mlen;
            TF fac;

            // Calculate geometric filter width, based on Deardorff (1980)
            const TF mlen0 = std::pow(dx*dy*dz[kstart], TF(1./3.));

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;

                    if ( bgradbot[ij] > 0 ) // Only if stably stratified, adapt length scale
                        mlen = cn * std::sqrt(a[ijk]) / std::sqrt(bgradbot[ij]);
                    else
                        mlen = mlen0;

                    fac  = std::min(mlen0, mlen);

                    if (sw_mason) // Apply Mason's wall correction here
                        fac = std::pow(TF(1.)/(TF(1.)/std::pow(fac, n_mason) + TF(1.)/
                                    (std::pow(Constants::kappa<TF>*(z[kstart]+z0m[ij]), n_mason))), TF(1.)/n_mason);

                    // Calculate eddy diffusivity for momentum.
                    evisc[ijk] = cm * fac * std::sqrt(a[ijk]);
                }

            for (int k=kstart+1; k<kend; ++k) // Counter starts at kstart (as sgstke is defined here)
            {
                // Calculate geometric filter width, based on Deardorff (1980)
                const TF mlen0 = std::pow(dx*dy*dz[k], TF(1./3.));

                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ij = i + j*jj;
                        const int ijk = i + j*jj + k*kk;

                        if (N2[ijk] > 0) // Only if stably stratified, adapt length scale
                            mlen = cn * std::sqrt(a[ijk]) / std::sqrt(N2[ijk]);
                        else
                            mlen = mlen0;

                        fac  = std::min(mlen0, mlen);

                        if (sw_mason) // Apply Mason's wall correction here
                            fac = std::pow(TF(1.)/(TF(1.)/std::pow(fac, n_mason) + TF(1.)/
                                        (std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);

                        // Calculate eddy diffusivity for momentum.
                        evisc[ijk] = cm * fac * std::sqrt(a[ijk]);
                    }
            }
        }

        boundary_cyclic.exec(evisc);
    }

    template<typename TF, Surface_model surface_model, bool sw_mason>
    void calc_evisc_heat(
            TF* const restrict evisch,
            const TF* const restrict evisc,
            const TF* const restrict a,
            const TF* const restrict N2,
            const TF* const restrict bgradbot,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* const restrict z0m,
            const TF dx, const TF dy,
            const TF cn, const TF ch1, const TF ch2,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells, const int ijcells,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int jj = icells;
        const int kk = ijcells;

        if (surface_model == Surface_model::Disabled)
             throw std::runtime_error("Resolved wall not supported in Deardorff SGSm.");
        else
        {
            // Variables for the wall damping and length scales
            const TF n_mason = TF(2.);
            TF mlen;
            TF fac;

            // Calculate geometric filter width, based on Deardorff (1980)
            const TF mlen0 = std::pow(dx*dy*dz[kstart], TF(1./3.));

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;

                    if ( bgradbot[ij] > 0 ) // Only if stably stratified, adapt length scale
                        mlen = cn * std::sqrt(a[ijk]) / std::sqrt(bgradbot[ij]);
                    else
                        mlen = mlen0;

                    fac  = std::min(mlen0, mlen);

                    if (sw_mason) // Apply Mason's wall correction here
                        fac = std::pow(TF(1.)/(TF(1.)/std::pow(fac, n_mason) + TF(1.)/
                                    (std::pow(Constants::kappa<TF>*(z[kstart]+z0m[ij]), n_mason))), TF(1.)/n_mason);

                    // Calculate eddy diffusivity for momentum.
                    evisch[ijk] = (ch1 + ch2 * fac / mlen0 ) * evisc[ijk];
                }

            for (int k=kstart+1; k<kend; ++k) // Counter starts at kstart (as sgstke is defined here)
            {
                // Calculate geometric filter width, based on Deardorff (1980)
                const TF mlen0 = std::pow(dx*dy*dz[k], TF(1./3.));

                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ij = i + j*jj;
                        const int ijk = i + j*jj + k*kk;

                        if ( N2[ijk] > 0 ) // Only if stably stratified, adapt length scale
                            mlen = cn * std::sqrt(a[ijk]) / std::sqrt(N2[ijk]);
                        else
                            mlen = mlen0;

                        fac  = std::min(mlen0, mlen);

                        if (sw_mason) // Apply Mason's wall correction here
                            fac = std::pow(TF(1.)/(TF(1.)/std::pow(fac, n_mason) + TF(1.)/
                                        (std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);

                        // Calculate eddy diffusivity for momentum.
                        evisch[ijk] = (ch1 + ch2 * fac / mlen0 ) * evisc[ijk];
                    }
            }
        }

        boundary_cyclic.exec(evisch);
    }

    template <typename TF>
    void sgstke_shear_tend(
            TF* const restrict at,
            const TF* const restrict a,
            const TF* const restrict evisc,
            const TF* const restrict strain2,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // Calculate shear production of SGS TKE based on Deardorff (1980)
                    // NOTE: `strain2` is defined/calculated as:
                    // S^2 = 0.5 * (dui/dxj + duj/dxi)^2 = dui/dxj * (dui/dxj + duj/dxi)
                    at[ijk] += evisc[ijk] * strain2[ijk];
                }
    }

    template <typename TF>
    void sgstke_buoy_tend(
            TF* const restrict at,
            const TF* const restrict a,
            const TF* const restrict evisch,
            const TF* const restrict N2,
            const TF* const restrict bgradbot,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;

                // Calculate buoyancy destruction of SGS TKE based on Deardorff (1980)
                at[ijk] -= evisch[ijk] * bgradbot[ij];
            }

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // Calculate buoyancy destruction of SGS TKE based on Deardorff (1980)
                    at[ijk] -= evisch[ijk] * N2[ijk];
                }
    }

    template <typename TF>
    void sgstke_diss_tend(
            TF* const restrict at,
            const TF* const restrict a,
            const TF* const restrict N2,
            const TF* const restrict bgradbot,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* const restrict z0m,
            const TF dx, const TF dy,
            const TF cn,
            const TF ce1, const TF ce2,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk,
            const bool sw_mason)
    {
        const TF n_mason = TF(2.);
        TF mlen ;
        TF fac  ;

        // Calculate geometric filter width, based on Deardorff (1980)
        const TF mlen0 = std::pow(dx*dy*dz[kstart], TF(1./3.));

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;

                if (bgradbot[ij] > 0) // Only if stably stratified, adapt length scale
                    mlen = cn * std::sqrt(a[ijk]) / std::sqrt(bgradbot[ij]);
                else
                    mlen = mlen0;

                fac  = std::min(mlen0, mlen);

                if (sw_mason) // Apply Mason's wall correction here
                    fac = std::pow(TF(1.)/(TF(1.)/std::pow(fac, n_mason) + TF(1.)/
                                (std::pow(Constants::kappa<TF>*(z[kstart]+z0m[ij]), n_mason))), TF(1.)/n_mason);

                // SvdL, 15-04-2023: quite strange (so check later), because why would (fac) be
                // altered by the Mason correction but mlen0 not? Why not both?
                // Calculate dissipation of SGS TKE based on Deardorff (1980)
                at[ijk] -= (ce1 + ce2 * fac / mlen0 ) * std::pow(a[ijk], TF(3./2.)) / fac;
            }

        for (int k=kstart+1; k<kend; ++k)
        {
            // Calculate geometric filter width, based on Deardorff (1980)
            const TF mlen0 = std::pow(dx*dy*dz[k], TF(1./3.));

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    const int ijk = i + j*jj + k*kk;

                    if (N2[ijk] > 0) // Only if stably stratified, adapt length scale
                        mlen = cn * std::sqrt(a[ijk]) / std::sqrt(N2[ijk]);
                    else
                        mlen = mlen0;

                    fac  = std::min(mlen0, mlen);

                    if (sw_mason) // Apply Mason's wall correction here
                        fac = std::pow(TF(1.)/(TF(1.)/std::pow(fac, n_mason) + TF(1.)/
                                    (std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);

                    // SvdL, 15-04-2023: quite strange (so check later), because why would (fac) be
                    // altered by the Mason correction but mlen0 not? Why not both?
                    // Calculate dissipation of SGS TKE based on Deardorff (1980)
                    at[ijk] -= (ce1 + ce2 * fac / mlen0 ) * std::pow(a[ijk], TF(3./2.)) / fac;
                }
        }
    }

    template <typename TF>
    void sgstke_diss_tend_neutral(
            TF* const restrict at,
            const TF* const restrict a,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* const restrict z0m,
            const TF dx, const TF dy,
            const TF ce1, const TF ce2,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk,
            const bool sw_mason)
    {
        const TF n_mason = TF(2.);
        TF fac;

        for (int k=kstart; k<kend; ++k)
        {
            // Calculate geometric filter width, based on Deardorff (1980)
            const TF mlen0 = std::pow(dx*dy*dz[k], TF(1./3.));

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + k*kk;

                    if (sw_mason) // Apply Mason's wall correction here
                        fac = std::pow(TF(1.)/(TF(1.)/std::pow(mlen0, n_mason) + TF(1.)/
                                    (std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);
                    else
                        fac = mlen0;

                    // SvdL, 15-04-2023: quite strange (so check later), because why would (fac) be
                    // altered by the Mason correction but mlen0 not? Why not both?
                    // Calculate dissipation of SGS TKE based on Deardorff (1980)
                    at[ijk] -= (ce1 + ce2 * fac / mlen0 ) * std::pow(a[ijk], TF(3./2.)) / fac ;
                }
        }
    }
}

template<typename TF>
Diff_deardorff<TF>::Diff_deardorff(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Boundary<TF>& boundaryin, Input& inputin) :
    Diff<TF>(masterin, gridin, fieldsin, boundaryin, inputin),
    boundary_cyclic(master, grid),
    field3d_operators(master, grid, fields)
{
    auto& gd = grid.get_grid_data();
    dnmax = inputin.get_item<TF>("diff", "dnmax", "", 0.4  );

    // Read constants of the Deardorff subgrid tke scheme
    ap    = inputin.get_item<TF>("diff", "ap"   , "", 1.5  );
    cf    = inputin.get_item<TF>("diff", "cf"   , "", 2.5  );
    ce1   = inputin.get_item<TF>("diff", "ce1"  , "", 0.19 );
    ce2   = inputin.get_item<TF>("diff", "ce2"  , "", 0.51 );
    cm    = inputin.get_item<TF>("diff", "cm"   , "", 0.12 );
    ch1   = inputin.get_item<TF>("diff", "ch1"  , "", 1.   );
    ch2   = inputin.get_item<TF>("diff", "ch2"  , "", 2.   );
    cn    = inputin.get_item<TF>("diff", "cn"   , "", 0.76 );

    const std::string group_name = "sgstke";

    // Set the switch between buoy/no buoy once
    const std::string sw_thermo = inputin.get_item<std::string>("thermo", "swthermo", "");
    sw_buoy = (sw_thermo == "0") ? false : true;

    // Set the switch for use of Mason's wall correction
    sw_mason = inputin.get_item<bool>("diff", "swmason", "", true);

    // Initialize field of SGS TKE
    fields.init_prognostic_field("sgstke", "SGS TKE", "m2 s-2", group_name, gd.sloc);

    // SvdL, 15-04-2023: If I remember correctly, exactly this was needed to avoid zero divisions somewhere? ...
    // Maybe it was just because of a call to a non-existing variable?
    fields.sp.at("sgstke")->visc = inputin.get_item<TF>("fields", "svisc", "sgstke");

    fields.init_diagnostic_field("evisc",  "Eddy viscosity for momentum", "m2 s-1", group_name, gd.sloc);

    // Add additional eddy viscosity for heat/scalars, if there is buoyancy
    if (sw_buoy)
        fields.init_diagnostic_field("eviscs", "Eddy viscosity for scalars", "m2 s-1",  group_name, gd.sloc);

    // Checks on input
    if (grid.get_spatial_order() != Grid_order::Second)
        throw std::runtime_error("Diff_deardorff only runs with second order grids.");
    if (boundary.get_switch() == "default")
        throw std::runtime_error("Diff_deardorff does not support resolved walls.");
}

template<typename TF>
Diff_deardorff<TF>::~Diff_deardorff()
{
}

template<typename TF>
void Diff_deardorff<TF>::init()
{
    boundary_cyclic.init();
}

template<typename TF>
Diffusion_type Diff_deardorff<TF>::get_switch() const
{
    return swdiff;
}

#ifndef USECUDA
template<typename TF>
unsigned long Diff_deardorff<TF>::get_time_limit(const unsigned long idt, const double dt)
{
    auto& gd = grid.get_grid_data();

    const TF tPr_dummy = 1;

    // When no buoyancy, use eddy viscosity for momentum.
    TF* evisc = !sw_buoy
        ? fields.sd.at("evisc")->fld.data()
        : fields.sd.at("eviscs")->fld.data();

    dnmul = dk::calc_dnmul<TF>(
            evisc,
            gd.dzi.data(),
            1./(gd.dx*gd.dx),
            1./(gd.dy*gd.dy),
            tPr_dummy,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    master.max(&dnmul, 1);

    // Avoid zero division.
    dnmul = std::max(Constants::dsmall, dnmul);

    return idt * dnmax / (dt * dnmul);
}
#endif

#ifndef USECUDA
template<typename TF>
double Diff_deardorff<TF>::get_dn(const double dt)
{
    auto& gd = grid.get_grid_data();

    const TF tPr_dummy = 1;

    // When no buoyancy, use eddy viscosity for momentum.
    TF* evisc = !sw_buoy
        ? fields.sd.at("evisc")->fld.data()
        : fields.sd.at("eviscs")->fld.data();

    dnmul = dk::calc_dnmul<TF>(
            evisc,
            gd.dzi.data(),
            1./(gd.dx*gd.dx),
            1./(gd.dy*gd.dy),
            tPr_dummy,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    master.max(&dnmul, 1);

    return dnmul*dt;
}
#endif

template<typename TF>
void Diff_deardorff<TF>::create(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // Get the maximum viscosity
    TF viscmax = fields.visc;
    for (auto& it : fields.sp)
        viscmax = std::max(it.second->visc, viscmax);

    // Calculate time step multiplier for diffusion number
    dnmul = 0;
    for (int k=gd.kstart; k<gd.kend; ++k)
        dnmul = std::max(dnmul, std::abs(viscmax * (1./(gd.dx*gd.dx) + 1./(gd.dy*gd.dy) + 1./(gd.dz[k]*gd.dz[k]))));

    create_stats(stats);

    // SvdL, 15-04-2023: If sgs tke from input >= 0, this function shouldn't be necessary.
    // However, sgs tke blows up, when this check is removed. Still don't fully understand why...
    // update 14-04-2023: probably related to the same type of crash Chiel experiences with sgs_diss
    check_for_minval<TF>(
            fields.sp.at("sgstke")->fld.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.jcells, gd.ijcells,
            boundary_cyclic);
}

#ifndef USECUDA
template<typename TF>
void Diff_deardorff<TF>::exec(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // Dummy tPr value for `diff_c`.
    const TF tPr_dummy = 1;

    dk::diff_u<TF, Surface_model::Enabled>(
            fields.mt.at("u")->fld.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dzhi.data(),
            1./gd.dx, 1./gd.dy,
            fields.sd.at("evisc")->fld.data(),
            fields.mp.at("u")->flux_bot.data(),
            fields.mp.at("u")->flux_top.data(),
            fields.rhoref.data(), fields.rhorefh.data(),
            fields.visc,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    dk::diff_v<TF, Surface_model::Enabled>(
            fields.mt.at("v")->fld.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dzhi.data(),
            1./gd.dx, 1./gd.dy,
            fields.sd.at("evisc")->fld.data(),
            fields.mp.at("v")->flux_bot.data(),
            fields.mp.at("v")->flux_top.data(),
            fields.rhoref.data(), fields.rhorefh.data(),
            fields.visc,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    dk::diff_w<TF>(
            fields.mt.at("w")->fld.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dzhi.data(),
            1./gd.dx, 1./gd.dy,
            fields.sd.at("evisc")->fld.data(),
            fields.rhoref.data(), fields.rhorefh.data(),
            fields.visc,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    for (auto it : fields.st)
    {
        if( it.first == "sgstke" ) // sgstke diffuses with eddy viscosity for momentum
        {
            dk::diff_c<TF, Surface_model::Enabled>(
                    it.second->fld.data(),
                    fields.sp.at(it.first)->fld.data(),
                    gd.dzi.data(), gd.dzhi.data(),
                    1./(gd.dx*gd.dx), 1./(gd.dy*gd.dy),
                    fields.sd.at("evisc")->fld.data(),
                    fields.sp.at(it.first)->flux_bot.data(),
                    fields.sp.at(it.first)->flux_top.data(),
                    fields.rhoref.data(), fields.rhorefh.data(),
                    tPr_dummy, fields.sp.at(it.first)->visc,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
        }
        else // all other scalars, normally diffuse with eddy viscosity for heat/scalars
        {
            if(!sw_buoy) // but not if there is no buoyancy (then eviscs not defined)
            {
                dk::diff_c<TF, Surface_model::Enabled>(
                        it.second->fld.data(),
                        fields.sp.at(it.first)->fld.data(),
                        gd.dzi.data(), gd.dzhi.data(),
                        1./(gd.dx*gd.dx), 1./(gd.dy*gd.dy),
                        fields.sd.at("evisc")->fld.data(),
                        fields.sp.at(it.first)->flux_bot.data(),
                        fields.sp.at(it.first)->flux_top.data(),
                        fields.rhoref.data(), fields.rhorefh.data(),
                        tPr_dummy, fields.sp.at(it.first)->visc,
                        gd.istart, gd.iend,
                        gd.jstart, gd.jend,
                        gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);
            }
            else // assume buoyancy calculation is needed
            {
                dk::diff_c<TF, Surface_model::Enabled>(
                        it.second->fld.data(),
                        fields.sp.at(it.first)->fld.data(),
                        gd.dzi.data(), gd.dzhi.data(),
                        1./(gd.dx*gd.dx), 1./(gd.dy*gd.dy),
                        fields.sd.at("eviscs")->fld.data(),
                        fields.sp.at(it.first)->flux_bot.data(),
                        fields.sp.at(it.first)->flux_top.data(),
                        fields.rhoref.data(), fields.rhorefh.data(),
                        tPr_dummy, fields.sp.at(it.first)->visc,
                        gd.istart, gd.iend,
                        gd.jstart, gd.jend,
                        gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);
            }
        }
    }

    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);

    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);
}

template<typename TF>
void Diff_deardorff<TF>::exec_viscosity(Stats<TF>& stats, Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();
    auto str2_tmp = fields.get_tmp();

    // Calculate strain rate using MO for velocity gradients lowest level.
    const std::vector<TF>& dudz = boundary.get_dudz();
    const std::vector<TF>& dvdz = boundary.get_dvdz();
    const std::vector<TF>& z0m = boundary.get_z0m();

    dk::calc_strain2<TF, Surface_model::Enabled>(
            str2_tmp->fld.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            dudz.data(),
            dvdz.data(),
            gd.z.data(),
            gd.dzi.data(),
            gd.dzhi.data(),
            1./gd.dx, 1./gd.dy,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    // Start with retrieving the stability information
    if (!sw_buoy)
    {
        auto evisc_neutral_wrapper = [&]<Surface_model surface_model, bool sw_mason>()
        {
            calc_evisc_neutral<TF, surface_model, sw_mason>(
                    fields.sd.at("evisc")->fld.data(),
                    fields.sp.at("sgstke")->fld.data(),
                    fields.mp.at("u")->fld.data(),
                    fields.mp.at("v")->fld.data(),
                    fields.mp.at("w")->fld.data(),
                    gd.z.data(), gd.dz.data(),
                    z0m.data(),
                    gd.dx, gd.dy,
                    this->cn, this->cm,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.jcells,
                    gd.ijcells,
                    boundary_cyclic);
        };

        // Calculate eddy viscosity using MO at lowest model level
        if (sw_mason)
            evisc_neutral_wrapper.template operator()<Surface_model::Enabled, true>();
        else
            evisc_neutral_wrapper.template operator()<Surface_model::Enabled, false>();

        sgstke_diss_tend_neutral<TF>(
                fields.st.at("sgstke")->fld.data(),
                fields.sp.at("sgstke")->fld.data(),
                gd.z.data(),
                gd.dz.data(),
                z0m.data(),
                gd.dx,
                gd.dy,
                this->ce1,
                this->ce2,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells,
                sw_mason);

                stats.calc_tend(*fields.st.at("sgstke"), tend_name_diss);
    }
    else
    {
        // Assume buoyancy calculation is needed
        auto buoy_tmp = fields.get_tmp();
        thermo.get_thermo_field(*buoy_tmp, "N2", false, false);
        const std::vector<TF>& dbdz = boundary.get_dbdz(); // SvdL, 19 April 2023: this should already be the "surface" N2

        auto evisc_wrapper = [&]<Surface_model surface_model, bool sw_mason>()
        {
            calc_evisc<TF, surface_model, sw_mason>(
                    fields.sd.at("evisc")->fld.data(),
                    fields.sp.at("sgstke")->fld.data(),
                    fields.mp.at("u")->fld.data(),
                    fields.mp.at("v")->fld.data(),
                    fields.mp.at("w")->fld.data(),
                    buoy_tmp->fld.data(),
                    dbdz.data(),
                    gd.z.data(), gd.dz.data(),
                    z0m.data(),
                    gd.dx, gd.dy,
                    this->cn, this->cm,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.jcells, gd.ijcells,
                    boundary_cyclic);
        };

        auto evisc_heat_wrapper = [&]<Surface_model surface_model, bool sw_mason>()
        {
            calc_evisc_heat<TF, surface_model, sw_mason>(
                    fields.sd.at("eviscs")->fld.data(),
                    fields.sd.at("evisc")->fld.data(),
                    fields.sp.at("sgstke")->fld.data(),
                    buoy_tmp->fld.data(),
                    dbdz.data(),
                    gd.z.data(), gd.dz.data(), z0m.data(),
                    gd.dx, gd.dy,
                    this->cn, this->ch1, this->ch2,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.jcells, gd.ijcells,
                    boundary_cyclic);
        };

        if (sw_mason)
        {
            evisc_wrapper.template operator()<Surface_model::Enabled, true>();
            evisc_heat_wrapper.template operator()<Surface_model::Enabled, true>();
        }
        else
        {
            evisc_wrapper.template operator()<Surface_model::Enabled, false>();
            evisc_heat_wrapper.template operator()<Surface_model::Enabled, false>();
        }

        // BvS: I left the tendency calculations of sgstke here; feels a bit strange
        // to calculate them in `exec_viscosity`, but otherwise strain^2 has to be
        // recalculated in diff->exec()...
        sgstke_buoy_tend<TF>(
                fields.st.at("sgstke")->fld.data(),
                fields.sp.at("sgstke")->fld.data(),
                fields.sd.at("eviscs")->fld.data(),
                buoy_tmp->fld.data(),
                dbdz.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_tend(*fields.st.at("sgstke"), tend_name_buoy);

        sgstke_diss_tend<TF>(
                fields.st.at("sgstke")->fld.data(),
                fields.sp.at("sgstke")->fld.data(),
                buoy_tmp->fld.data(),
                dbdz.data(),
                gd.z.data(),
                gd.dz.data(),
                z0m.data(),
                gd.dx, gd.dy,
                this->cn,
                this->ce1,
                this->ce2,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells, sw_mason); ///< SvdL, 14-04-2023: not the nicest, see above);

        stats.calc_tend(*fields.st.at("sgstke"), tend_name_diss);

        fields.release_tmp(buoy_tmp);
    }

    sgstke_shear_tend<TF>(
            fields.st.at("sgstke")->fld.data(),
            fields.sp.at("sgstke")->fld.data(),
            fields.sd.at("evisc")->fld.data(),
            str2_tmp->fld.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    stats.calc_tend(*fields.st.at("sgstke"), tend_name_shear);

    // Release temporary fields
    fields.release_tmp(str2_tmp);
}
#endif

template<typename TF>
void Diff_deardorff<TF>::create_stats(Stats<TF>& stats)
{
    const std::string group_name_tke = "sgstke";
    const std::string group_name_default = "default";

    // Add variables to the statistics
    if (stats.get_switch())
    {
        // Always add statistics of eddy viscosity for momentum (!)
        stats.add_profs(*fields.sd.at("evisc"), "z", {"mean", "2"}, group_name_default);

        // Add shear and dissipation of sgstke to the list of tendencies
        stats.add_tendency(*fields.st.at("sgstke"), "z", tend_name_shear, tend_longname_shear);
        stats.add_tendency(*fields.st.at("sgstke"), "z", tend_name_diss, tend_longname_diss);

        // Add additional profile of eddy viscosity for heat/scalars and tendency of buoyancy production of sgstke
        if (sw_buoy)
        {
            stats.add_profs(*fields.sd.at("eviscs"), "z", {"mean", "2"}, group_name_default);
            stats.add_tendency(*fields.st.at("sgstke"), "z", tend_name_buoy, tend_longname_buoy);
        }

        stats.add_tendency(*fields.mt.at("u"), "z",  tend_name, tend_longname);
        stats.add_tendency(*fields.mt.at("v"), "z",  tend_name, tend_longname);
        stats.add_tendency(*fields.mt.at("w"), "zh", tend_name, tend_longname);

        for (auto it : fields.st)
            stats.add_tendency(*it.second, "z", tend_name, tend_longname);
    }
}

template<typename TF>
void Diff_deardorff<TF>::exec_stats(Stats<TF>& stats, Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();

    const TF no_offset = 0.;
    const TF no_threshold = 0.;

    stats.calc_stats("evisc", *fields.sd.at("evisc"), no_offset, no_threshold);
    if (sw_buoy)
        stats.calc_stats("eviscs", *fields.sd.at("eviscs"), no_offset, no_threshold);
}

template<typename TF>
void Diff_deardorff<TF>::diff_flux(
        Field3d<TF>& restrict out, const Field3d<TF>& restrict fld_in)
{
    auto& gd = grid.get_grid_data();

    const TF tPr_dummy = 1;

    // SvdL. 15-04-2023: still check if these boundary fluxes for sgstke are correct?
    // Calculate the boundary fluxes.
    dk::calc_diff_flux_bc(
            out.fld.data(), fld_in.flux_bot.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.icells, gd.ijcells);

    dk::calc_diff_flux_bc(
            out.fld.data(), fld_in.flux_top.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kend, gd.icells, gd.ijcells);

    // Calculate the interior.
    if (fld_in.loc[0] == 1)
        dk::calc_diff_flux_u<TF, Surface_model::Enabled>(
                out.fld.data(), fld_in.fld.data(),
                fields.mp.at("w")->fld.data(),
                fields.sd.at("evisc")->fld.data(),
                gd.dxi, gd.dzhi.data(),
                fields.visc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

    else if (fld_in.loc[1] == 1)
        dk::calc_diff_flux_v<TF, Surface_model::Enabled>(
                out.fld.data(), fld_in.fld.data(),
                fields.mp.at("w")->fld.data(),
                fields.sd.at("evisc")->fld.data(),
                gd.dyi, gd.dzhi.data(),
                fields.visc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    else
    {
        // SvdL, 14-04-2023: if no buoyancy scalars diffuse with eddy viscosity for momentum,
        // sgstke and w always diffuse with this one.
        std::string varname = fld_in.name;
        if (!sw_buoy || varname == "sgstke" || varname == "w")
            dk::calc_diff_flux_c<TF, Surface_model::Enabled>(
                    out.fld.data(), fld_in.fld.data(),
                    fields.sd.at("evisc")->fld.data(),
                    gd.dzhi.data(),
                    tPr_dummy, fld_in.visc,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
        else
            dk::calc_diff_flux_c<TF, Surface_model::Enabled>(
                    out.fld.data(), fld_in.fld.data(),
                    fields.sd.at("eviscs")->fld.data(),
                    gd.dzhi.data(),
                    tPr_dummy, fld_in.visc,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
    }
}

template class Diff_deardorff<double>;
template class Diff_deardorff<float>;
