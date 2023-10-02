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

#include "diff_smag2.h"
#include "diff_kernels.h"

namespace
{
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;
    namespace dk = Diff_kernels;

    template <typename TF, Surface_model surface_model, bool sw_mason>
    void calc_evisc_neutral(
            TF* const restrict evisc,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict ufluxbot,
            const TF* const restrict vfluxbot,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* const restrict dzhi,
            const TF* const restrict z0m,
            const TF dx, const TF dy, const TF zsize,
            const TF cs, const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells, const int ijcells,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int jj = icells;
        const int kk = ijcells;

        TF mlen;

        // Wall damping constant.
        constexpr TF n_mason = TF(1.);
        constexpr TF A_vandriest = TF(26.);

        if (surface_model == Surface_model::Disabled)
        {
            for (int k=kstart; k<kend; ++k)
            {
                // const TF mlen_wall = Constants::kappa<TF>*std::min(z[k], zsize-z[k]);
                const TF mlen_smag = cs*std::pow(dx*dy*dz[k], TF(1./3.));

                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk_bot = i + j*jj + kstart*kk;
                        const int ijk_top = i + j*jj + kend*kk;

                        const TF u_tau_bot = std::pow(
                                fm::pow2( visc*(u[ijk_bot] - u[ijk_bot-kk] )*dzhi[kstart] )
                              + fm::pow2( visc*(v[ijk_bot] - v[ijk_bot-kk] )*dzhi[kstart] ), TF(0.25) );
                        const TF u_tau_top = std::pow(
                                fm::pow2( visc*(u[ijk_top] - u[ijk_top-kk] )*dzhi[kend] )
                              + fm::pow2( visc*(v[ijk_top] - v[ijk_top-kk] )*dzhi[kend] ), TF(0.25) );

                        const TF fac_bot = TF(1.) - std::exp( -(       z[k] *u_tau_bot) / (A_vandriest*visc) );
                        const TF fac_top = TF(1.) - std::exp( -((zsize-z[k])*u_tau_top) / (A_vandriest*visc) );
                        const TF fac = std::min( fac_bot, fac_top );

                        const int ijk = i + j*jj + k*kk;
                        evisc[ijk] = fm::pow2(fac * mlen_smag) * std::sqrt(evisc[ijk]);
                    }
            }

            // For a resolved wall the viscosity at the wall is needed. For now, assume that the eddy viscosity
            // is mirrored around the surface.
            const int kb = kstart;
            const int kt = kend-1;
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ijkb = i + j*jj + kb*kk;
                    const int ijkt = i + j*jj + kt*kk;
                    evisc[ijkb-kk] = evisc[ijkb];
                    evisc[ijkt+kk] = evisc[ijkt];
                }
        }
        else
        {
            for (int k=kstart; k<kend; ++k)
            {
                // Calculate smagorinsky constant times filter width squared, use wall damping according to Mason's paper.
                const TF mlen0 = cs*std::pow(dx*dy*dz[k], TF(1./3.));

                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ij  = i + j*jj;
                        const int ijk = i + j*jj + k*kk;

                        if (sw_mason) // Apply Mason's wall correction
                            mlen = std::pow(TF(1.)/(TF(1.)/std::pow(mlen0, n_mason) + TF(1.)/
                                        (std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);
                        else
                            mlen = mlen0;

                        evisc[ijk] = fm::pow2(mlen) * std::sqrt(evisc[ijk]);
                    }
            }
        }

        boundary_cyclic.exec(evisc);
    }

    template<typename TF, Surface_model surface_model, bool sw_mason>
    void calc_evisc(
            TF* const restrict evisc,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict N2,
            const TF* const restrict bgradbot,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* const restrict dzi,
            const TF* const restrict z0m,
            const TF dx, const TF dy,
            const TF cs, const TF tPr,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells, const int ijcells,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int jj = icells;
        const int kk = ijcells;

        TF mlen;

        if (surface_model == Surface_model::Disabled)
        {
            for (int k=kstart; k<kend; ++k)
            {
                // calculate smagorinsky constant times filter width squared, do not use wall damping with resolved walls.
                const TF mlen = cs*std::pow(dx*dy*dz[k], TF(1./3.));
                const TF fac = fm::pow2(mlen);

                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;

                        // Add the buoyancy production to the TKE
                        TF RitPrratio = N2[ijk] / evisc[ijk] / tPr;
                        RitPrratio = std::min(RitPrratio, TF(1.-Constants::dsmall));

                        evisc[ijk] = fac * std::sqrt(evisc[ijk]) * std::sqrt(TF(1.)-RitPrratio);
                    }
            }

            // For a resolved wall the viscosity at the wall is needed. For now, assume that the eddy viscosity
            // is mirrored over the surface.
            const int kb = kstart;
            const int kt = kend-1;
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ijkb = i + j*jj + kb*kk;
                    const int ijkt = i + j*jj + kt*kk;
                    evisc[ijkb-kk] = evisc[ijkb];
                    evisc[ijkt+kk] = evisc[ijkt];
                }
        }
        else
        {
            // Variables for the wall damping.
            const TF n_mason = 2.;

            // Bottom boundary, here strain is fully parametrized using MO.
            // Calculate smagorinsky constant times filter width squared, use wall damping according to Mason.
            const TF mlen0 = cs*std::pow(dx*dy*dz[kstart], TF(1./3.));

            for (int j=jstart; j<jend; ++j)
            {
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;

                    // TODO use the thermal expansion coefficient from the input later, what to do if there is no buoyancy?
                    // Add the buoyancy production to the TKE
                    TF RitPrratio = bgradbot[ij] / evisc[ijk] / tPr;
                    RitPrratio = std::min(RitPrratio, TF(1.-Constants::dsmall));

                    if (sw_mason) // Apply Mason's wall correction
                        mlen = std::pow(TF(1.)/(TF(1.)/std::pow(mlen0, n_mason) + TF(1.)/
                                    (std::pow(Constants::kappa<TF>*(z[kstart]+z0m[ij]), n_mason))), TF(1.)/n_mason);
                    else
                        mlen = mlen0;

                    evisc[ijk] = fm::pow2(mlen) * std::sqrt(evisc[ijk]) * std::sqrt(TF(1.)-RitPrratio);
                }
            }

            for (int k=kstart+1; k<kend; ++k)
            {
                // Calculate smagorinsky constant times filter width squared, use wall damping according to Mason
                const TF mlen0 = cs*std::pow(dx*dy*dz[k], TF(1./3.));

                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ij  = i + j*jj;
                        const int ijk = i + j*jj + k*kk;

                        // Add the buoyancy production to the TKE
                        TF RitPrratio = N2[ijk] / evisc[ijk] / tPr;
                        RitPrratio = std::min(RitPrratio, TF(1.-Constants::dsmall));

                        if (sw_mason) // Apply Mason's wall correction
                            mlen = std::pow(TF(1.)/(TF(1.)/std::pow(mlen0, n_mason) + TF(1.)/
                                        (std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);
                        else
                            mlen = mlen0;

                        evisc[ijk] = fm::pow2(mlen) * std::sqrt(evisc[ijk]) * std::sqrt(TF(1.)-RitPrratio);
                    }
            }
        }

        boundary_cyclic.exec(evisc);
    }
} // End namespace.

template<typename TF>
Diff_smag2<TF>::Diff_smag2(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Boundary<TF>& boundaryin, Input& inputin) :
    Diff<TF>(masterin, gridin, fieldsin, boundaryin, inputin),
    boundary_cyclic(master, grid),
    field3d_operators(master, grid, fields)
{
    auto& gd = grid.get_grid_data();
    dnmax = inputin.get_item<TF>("diff", "dnmax", "", 0.4  );
    cs    = inputin.get_item<TF>("diff", "cs"   , "", 0.23 );
    tPr   = inputin.get_item<TF>("diff", "tPr"  , "", 1./3.);

    const std::string group_name = "default";

    // Set the switch for use of Mason's wall correction
    sw_mason = inputin.get_item<bool>("diff", "swmason", "", true);

    fields.init_diagnostic_field("evisc", "Eddy viscosity", "m2 s-1", group_name, gd.sloc);

    if (grid.get_spatial_order() != Grid_order::Second)
        throw std::runtime_error("Diff_smag2 only runs with second order grids");
}

template<typename TF>
Diff_smag2<TF>::~Diff_smag2()
{
}

template<typename TF>
void Diff_smag2<TF>::init()
{
    boundary_cyclic.init();
}

template<typename TF>
Diffusion_type Diff_smag2<TF>::get_switch() const
{
    return swdiff;
}

#ifndef USECUDA
template<typename TF>
unsigned long Diff_smag2<TF>::get_time_limit(const unsigned long idt, const double dt)
{
    auto& gd = grid.get_grid_data();

    double dnmul = dk::calc_dnmul<TF>(
        fields.sd.at("evisc")->fld.data(),
        gd.dzi.data(),
        1./(gd.dx*gd.dx),
        1./(gd.dy*gd.dy),
        tPr,
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
double Diff_smag2<TF>::get_dn(const double dt)
{
    auto& gd = grid.get_grid_data();

    double dnmul = dk::calc_dnmul<TF>(
        fields.sd.at("evisc")->fld.data(),
        gd.dzi.data(),
        1./(gd.dx*gd.dx),
        1./(gd.dy*gd.dy),
        tPr,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    master.max(&dnmul, 1);

    return dnmul*dt;
}
#endif

template<typename TF>
void Diff_smag2<TF>::create(Stats<TF>& stats, const bool cold_start)
{
    if (cold_start)
        return;

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
}

#ifndef USECUDA
template<typename TF>
void Diff_smag2<TF>::exec(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    auto diff_wrapper = [&]<Surface_model surface_model>()
    {
        dk::diff_u<TF, surface_model>(
                fields.mt.at("u")->fld.data(),
                fields.mp.at("u")->fld.data(),
                fields.mp.at("v")->fld.data(),
                fields.mp.at("w")->fld.data(),
                gd.dzi.data(), gd.dzhi.data(),
                1./gd.dx, 1./gd.dy,
                fields.sd.at("evisc")->fld.data(),
                fields.mp.at("u")->flux_bot.data(),
                fields.mp.at("u")->flux_top.data(),
                fields.rhoref.data(),
                fields.rhorefh.data(),
                fields.visc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        dk::diff_v<TF, surface_model>(
                fields.mt.at("v")->fld.data(),
                fields.mp.at("u")->fld.data(),
                fields.mp.at("v")->fld.data(),
                fields.mp.at("w")->fld.data(),
                gd.dzi.data(), gd.dzhi.data(),
                1./gd.dx, 1./gd.dy,
                fields.sd.at("evisc")->fld.data(),
                fields.mp.at("v")->flux_bot.data(),
                fields.mp.at("v")->flux_top.data(),
                fields.rhoref.data(),
                fields.rhorefh.data(),
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
                fields.rhoref.data(),
                fields.rhorefh.data(),
                fields.visc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        for (auto it : fields.st)
        {
            dk::diff_c<TF, surface_model>(
                    it.second->fld.data(),
                    fields.sp.at(it.first)->fld.data(),
                    gd.dzi.data(), gd.dzhi.data(),
                    1./(gd.dx*gd.dx), 1./(gd.dy*gd.dy),
                    fields.sd.at("evisc")->fld.data(),
                    fields.sp.at(it.first)->flux_bot.data(),
                    fields.sp.at(it.first)->flux_top.data(),
                    fields.rhoref.data(),
                    fields.rhorefh.data(),
                    tPr,
                    fields.sp.at(it.first)->visc,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
        }
    };

    if (boundary.get_switch() != "default")
        diff_wrapper.template operator()<Surface_model::Enabled>();
    else
        diff_wrapper.template operator()<Surface_model::Disabled>();

    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);

    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);
}

template<typename TF>
void Diff_smag2<TF>::exec_viscosity(Stats<TF>&, Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();

    auto strain2_wrapper = [&]<Surface_model surface_model>(
            const TF* const restrict dudz,
            const TF* const restrict dvdz)
    {
        dk::calc_strain2<TF, Surface_model::Enabled>(
                fields.sd.at("evisc")->fld.data(),
                fields.mp.at("u")->fld.data(),
                fields.mp.at("v")->fld.data(),
                fields.mp.at("w")->fld.data(),
                dudz,
                dvdz,
                gd.z.data(),
                gd.dzi.data(),
                gd.dzhi.data(),
                1./gd.dx, 1./gd.dy,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    };

    if (boundary.get_switch() != "default")
    {
        // Calculate strain rate using MO for velocity gradients lowest level.
        const std::vector<TF>& dudz = boundary.get_dudz();
        const std::vector<TF>& dvdz = boundary.get_dvdz();

        strain2_wrapper.template operator()<Surface_model::Enabled>(dudz.data(), dvdz.data());
    }
    else
        strain2_wrapper.template operator()<Surface_model::Enabled>(nullptr, nullptr);

    // Start with retrieving the stability information
    if (thermo.get_switch() == Thermo_type::Disabled)
    {
        auto evisc_wrapper = [&]<Surface_model surface_model, bool sw_mason>(
                const TF* const restrict z0m)
        {
            calc_evisc_neutral<TF, surface_model, sw_mason>(
                    fields.sd.at("evisc")->fld.data(),
                    fields.mp.at("u")->fld.data(),
                    fields.mp.at("v")->fld.data(),
                    fields.mp.at("w")->fld.data(),
                    fields.mp.at("u")->flux_bot.data(),
                    fields.mp.at("v")->flux_bot.data(),
                    gd.z.data(), gd.dz.data(),
                    gd.dzhi.data(), z0m,
                    gd.dx, gd.dy, gd.zsize,
                    this->cs,
                    fields.visc,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.jcells, gd.ijcells,
                    boundary_cyclic);
        };

        if (boundary.get_switch() != "default")
        {
            const std::vector<TF>& z0m = boundary.get_z0m();

            if (sw_mason)
                evisc_wrapper.template operator()<Surface_model::Enabled, true>(z0m.data());
            else
                evisc_wrapper.template operator()<Surface_model::Enabled, false>(z0m.data());
        }
        else
        {
            if (sw_mason)
                evisc_wrapper.template operator()<Surface_model::Disabled, true>(nullptr);
            else
                evisc_wrapper.template operator()<Surface_model::Disabled, false>(nullptr);
        }
    }
    // assume buoyancy calculation is needed
    else
    {
        // Store the buoyancy flux in tmp1
        auto& gd = grid.get_grid_data();
        auto buoy_tmp = fields.get_tmp();
        auto tmp = fields.get_tmp();

        thermo.get_thermo_field(*buoy_tmp, "N2", false, false);

        auto evisc_wrapper = [&]<Surface_model surface_model, bool sw_mason>(
                const TF* const restrict dbdz,
                const TF* const restrict z0m)
        {
            calc_evisc<TF, surface_model, sw_mason>(
                    fields.sd.at("evisc")->fld.data(),
                    fields.mp.at("u")->fld.data(),
                    fields.mp.at("v")->fld.data(),
                    fields.mp.at("w")->fld.data(),
                    buoy_tmp->fld.data(),
                    dbdz,
                    gd.z.data(),
                    gd.dz.data(),
                    gd.dzi.data(),
                    z0m,
                    gd.dx, gd.dy,
                    this->cs, this->tPr,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.jcells, gd.ijcells,
                    boundary_cyclic);
        };

        if (boundary.get_switch() != "default")
        {
            const std::vector<TF>& z0m = boundary.get_z0m();
            const std::vector<TF>& dbdz = boundary.get_dbdz();

            if (sw_mason)
                evisc_wrapper.template operator()<Surface_model::Enabled, true>(dbdz.data(), z0m.data());
            else
                evisc_wrapper.template operator()<Surface_model::Enabled, false>(dbdz.data(), z0m.data());
        }
        else
        {
            if (sw_mason)
                evisc_wrapper.template operator()<Surface_model::Disabled, true>(nullptr, nullptr);
            else
                evisc_wrapper.template operator()<Surface_model::Disabled, false>(nullptr, nullptr);
        }

        fields.release_tmp(buoy_tmp);
        fields.release_tmp(tmp);
    }
}
#endif

template<typename TF>
void Diff_smag2<TF>::create_stats(Stats<TF>& stats)
{
    const std::string group_name = "default";

    // Add variables to the statistics
    if (stats.get_switch())
    {
        stats.add_profs(*fields.sd.at("evisc"), "z", {"mean", "2"}, group_name);
        stats.add_tendency(*fields.mt.at("u"), "z", tend_name, tend_longname);
        stats.add_tendency(*fields.mt.at("v"), "z", tend_name, tend_longname);
        stats.add_tendency(*fields.mt.at("w"), "zh", tend_name, tend_longname);

        for (auto it : fields.st)
            stats.add_tendency(*it.second, "z", tend_name, tend_longname);
    }
}

template<typename TF>
void Diff_smag2<TF>::exec_stats(Stats<TF>& stats, Thermo<TF>& thermo)
{
    const TF no_offset = 0.;
    const TF no_threshold = 0.;
    stats.calc_stats("evisc", *fields.sd.at("evisc"), no_offset, no_threshold);
}

template<typename TF>
void Diff_smag2<TF>::diff_flux(Field3d<TF>& restrict out, const Field3d<TF>& restrict fld_in)
{
    auto& gd = grid.get_grid_data();

    auto diff_flux_wrapper = [&]<Surface_model surface_model>()
    {
        // Calculate the interior.
        if (fld_in.loc[0] == 1)
            dk::calc_diff_flux_u<TF, surface_model>(
                    out.fld.data(),
                    fld_in.fld.data(),
                    fields.mp.at("w")->fld.data(),
                    fields.sd.at("evisc")->fld.data(),
                    gd.dxi, gd.dzhi.data(),
                    fields.visc,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

        else if (fld_in.loc[1] == 1)
            dk::calc_diff_flux_v<TF, surface_model>(
                    out.fld.data(),
                    fld_in.fld.data(),
                    fields.mp.at("w")->fld.data(),
                    fields.sd.at("evisc")->fld.data(),
                    gd.dyi, gd.dzhi.data(),
                    fields.visc,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

        else
            dk::calc_diff_flux_c<TF, surface_model>(
                    out.fld.data(),
                    fld_in.fld.data(),
                    fields.sd.at("evisc")->fld.data(),
                    gd.dzhi.data(),
                    tPr, fld_in.visc,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
    };

    if (boundary.get_switch() != "default")
    {
        // Calculate the boundary fluxes.
        dk::calc_diff_flux_bc(
                out.fld.data(),
                fld_in.flux_bot.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart,
                gd.icells, gd.ijcells);

        dk::calc_diff_flux_bc(
                out.fld.data(),
                fld_in.flux_top.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kend,
                gd.icells, gd.ijcells);

        // Calculate interior:
        diff_flux_wrapper.template operator()<Surface_model::Enabled>();
    }
    else
    {
        // Calculate boundary and interior:
        diff_flux_wrapper.template operator()<Surface_model::Disabled>();
    }
}


#ifdef FLOAT_SINGLE
template class Diff_smag2<float>;
#else
template class Diff_smag2<double>;
#endif
