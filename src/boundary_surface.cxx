/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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
#include <algorithm>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "diff.h"
#include "boundary_surface.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "master.h"
#include "cross.h"
#include "monin_obukhov.h"

namespace
{
    // Make a shortcut in the file scope.
    namespace most = Monin_obukhov;

    // Size of the lookup table.
    const int nzL = 10000; // Size of the lookup table for MO iterations.

    template<typename TF>
    TF find_zL(const float* const restrict zL, const float* const restrict f,
               int& n, const float Ri)
    {
        // Determine search direction. All checks are at float accuracy.
        if ( (f[n]-Ri) > 0.f )
            while ( (f[n-1]-Ri) > 0.f && n > 0.f) { --n; }
        else
            while ( (f[n]-Ri) < 0.f && n < (nzL-1) ) { ++n; }

        const TF zL0 = (n == 0 || n == nzL-1) ? zL[n] : zL[n-1] + (Ri-f[n-1]) / (f[n]-f[n-1]) * (zL[n]-zL[n-1]);

        return zL0;
    }

    template<typename TF>
    TF calc_obuk_noslip_flux(const float* restrict zL, const float* restrict f,
                             int& n,
                             const TF du, const TF bfluxbot, const TF zsl)
    {
        // Calculate the appropriate Richardson number and reduce precision.
        const float Ri = -Constants::kappa<TF> * bfluxbot * zsl / std::pow(du, 3);
        return zsl/find_zL<TF>(zL, f, n, Ri);
    }

    template<typename TF>
    TF calc_obuk_noslip_dirichlet(const float* restrict zL, const float* restrict f,
                                  int& n,
                                  const TF du, const TF db, const TF zsl)
    {
        // Calculate the appropriate Richardson number and reduce precision.
        const float Ri = Constants::kappa<TF> * db * zsl / std::pow(du, 2);
        return zsl/find_zL<TF>(zL, f, n, Ri);
    }

    template<typename TF>
    void set_bc(TF* const restrict a, TF* const restrict agrad, TF* const restrict aflux,
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
    void stability(TF* restrict ustar, TF* restrict obuk, TF* restrict bfluxbot,
                   TF* restrict u, TF* restrict v, TF* restrict b,
                   TF* restrict ubot , TF* restrict vbot, TF* restrict bbot,
                   TF* restrict dutot, const TF* restrict z,
                   const float* zL_sl, const float* f_sl, int* nobuk,
                   const TF z0m, const TF z0h,
                   const int istart, const int iend, const int jstart, const int jend, const int kstart,
                   const int icells, const int jcells, const int kk,
                   Boundary_type mbcbot, Boundary_type thermobc,
                   Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int ii = 1;
        const int jj = icells;

        // Calculate total wind.
        TF du2;
        const TF minval = 1.e-1;

        // First, interpolate the wind to the scalar location.
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                du2 = std::pow(static_cast<TF>(0.5)*(u[ijk] + u[ijk+ii]) - static_cast<TF>(0.5)*(ubot[ij] + ubot[ij+ii]), static_cast<TF>(2))
                    + std::pow(static_cast<TF>(0.5)*(v[ijk] + v[ijk+jj]) - static_cast<TF>(0.5)*(vbot[ij] + vbot[ij+jj]), static_cast<TF>(2));
                // prevent the absolute wind gradient from reaching values less than 0.01 m/s,
                // otherwise evisc at k = kstart blows up
                dutot[ij] = std::max(std::pow(du2, static_cast<TF>(0.5)), minval);
            }

        boundary_cyclic.exec_2d(dutot);

        // calculate Obukhov length
        // case 1: fixed buoyancy flux and fixed ustar
        if (mbcbot == Boundary_type::Ustar_type && thermobc == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    obuk[ij] = -std::pow(ustar[ij], static_cast<TF>(3)) / (Constants::kappa<TF>*bfluxbot[ij]);
                }
        }
        // case 2: fixed buoyancy surface value and free ustar
        else if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    obuk [ij] = calc_obuk_noslip_flux(zL_sl, f_sl, nobuk[ij], dutot[ij], bfluxbot[ij], z[kstart]);
                    ustar[ij] = dutot[ij] * most::fm(z[kstart], z0m, obuk[ij]);
                }
        }
        else if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    const TF db = b[ijk] - bbot[ij];
                    obuk [ij] = calc_obuk_noslip_dirichlet(zL_sl, f_sl, nobuk[ij], dutot[ij], db, z[kstart]);
                    ustar[ij] = dutot[ij] * most::fm(z[kstart], z0m, obuk[ij]);
                }
        }
    }

    template<typename TF>
    void stability_neutral(
            TF* restrict ustar, TF* restrict obuk,
            TF* restrict u, TF* restrict v,
            TF* restrict ubot , TF* restrict vbot,
            TF* restrict dutot, const TF* restrict z,
            const TF z0m,
            const int istart, const int iend, const int jstart, const int jend, const int kstart,
            const int icells, const int jcells, const int kk,
            Boundary_type mbcbot,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int ii = 1;
        const int jj = icells;

        // calculate total wind
        TF du2;
        const TF minval = 1.e-1;

        // first, interpolate the wind to the scalar location
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                du2 = std::pow(static_cast<TF>(0.5)*(u[ijk] + u[ijk+ii]) - static_cast<TF>(0.5)*(ubot[ij] + ubot[ij+ii]), static_cast<TF>(2.))
                    + std::pow(static_cast<TF>(0.5)*(v[ijk] + v[ijk+jj]) - static_cast<TF>(0.5)*(vbot[ij] + vbot[ij+jj]), static_cast<TF>(2.));
                // prevent the absolute wind gradient from reaching values less than 0.01 m/s,
                // otherwise evisc at k = kstart blows up
                dutot[ij] = std::max(std::pow(du2, static_cast<TF>(0.5)), minval);
            }

        boundary_cyclic.exec_2d(dutot);

        // set the Obukhov length to a very large negative number
        // case 1: fixed buoyancy flux and fixed ustar
        if (mbcbot == Boundary_type::Ustar_type)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    obuk[ij] = -Constants::dbig;
                }
        }
        // case 2: free ustar
        else if (mbcbot == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    obuk [ij] = -Constants::dbig;
                    ustar[ij] = dutot[ij] * most::fm(z[kstart], z0m, obuk[ij]);
                }
        }
    }

    template<typename TF>
    void surfm(TF* restrict ustar, TF* restrict obuk,
               TF* restrict u, TF* restrict ubot, TF* restrict ugradbot, TF* restrict ufluxbot,
               TF* restrict v, TF* restrict vbot, TF* restrict vgradbot, TF* restrict vfluxbot,
               const TF zsl, const TF z0m, const Boundary_type bcbot,
               const int istart, const int iend, const int jstart, const int jend, const int kstart,
               const int icells, const int jcells, const int kk,
               Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int ii = 1;
        const int jj = icells;

        // the surface value is known, calculate the flux and gradient
        if (bcbot == Boundary_type::Dirichlet_type)
        {
            // first calculate the surface value
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;

                    // interpolate the whole stability function rather than ustar or obuk
                    ufluxbot[ij] = -(u[ijk]-ubot[ij])*static_cast<TF>(0.5)*(ustar[ij-ii]*most::fm(zsl, z0m, obuk[ij-ii]) + ustar[ij]*most::fm(zsl, z0m, obuk[ij]));
                    vfluxbot[ij] = -(v[ijk]-vbot[ij])*static_cast<TF>(0.5)*(ustar[ij-jj]*most::fm(zsl, z0m, obuk[ij-jj]) + ustar[ij]*most::fm(zsl, z0m, obuk[ij]));
                }

            boundary_cyclic.exec_2d(ufluxbot);
            boundary_cyclic.exec_2d(vfluxbot);
        }
        // the flux is known, calculate the surface value and gradient
        else if (bcbot == Boundary_type::Ustar_type)
        {
            // first redistribute ustar over the two flux components
            const TF minval = 1.e-2;

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;

                    // CvH: now, this is ugly...
                    const TF one  = static_cast<TF>(1);
                    const TF two  = static_cast<TF>(2);
                    const TF four = static_cast<TF>(4);

                    const TF vonu2 = std::max(minval, static_cast<TF>(0.25)*( std::pow(v[ijk-ii]-vbot[ij-ii], two) + std::pow(v[ijk-ii+jj]-vbot[ij-ii+jj], two)
                                                                            + std::pow(v[ijk   ]-vbot[ij   ], two) + std::pow(v[ijk   +jj]-vbot[ij   +jj], two)) );
                    const TF uonv2 = std::max(minval, static_cast<TF>(0.25)*( std::pow(u[ijk-jj]-ubot[ij-jj], two) + std::pow(u[ijk+ii-jj]-ubot[ij+ii-jj], two)
                                                                            + std::pow(u[ijk   ]-ubot[ij   ], two) + std::pow(u[ijk+ii   ]-ubot[ij+ii   ], two)) );
                    const TF u2 = std::max(minval, std::pow(u[ijk]-ubot[ij], two) );
                    const TF v2 = std::max(minval, std::pow(v[ijk]-vbot[ij], two) );
                    const TF ustaronu4 = static_cast<TF>(0.5)*(std::pow(ustar[ij-ii], four) + std::pow(ustar[ij], four));
                    const TF ustaronv4 = static_cast<TF>(0.5)*(std::pow(ustar[ij-jj], four) + std::pow(ustar[ij], four));
                    ufluxbot[ij] = -copysign(one, u[ijk]-ubot[ij]) * std::pow(ustaronu4 / (one + vonu2 / u2), static_cast<TF>(0.5));
                    vfluxbot[ij] = -copysign(one, v[ijk]-vbot[ij]) * std::pow(ustaronv4 / (one + uonv2 / v2), static_cast<TF>(0.5));
                }

            boundary_cyclic.exec_2d(ufluxbot);
            boundary_cyclic.exec_2d(vfluxbot);

            // CvH: I think that the problem is not closed, since both the fluxes and the surface values
            // of u and v are unknown. You have to assume a no slip in order to get the fluxes and therefore
            // should not update the surface values with those that belong to the flux. This procedure needs
            // to be checked more carefully.
            /*
            // calculate the surface values
            for (int j=grid->jstart; j<grid->jend; ++j)
                #pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                ij  = i + j*jj;
                ijk = i + j*jj + kstart*kk;
                // interpolate the whole stability function rather than ustar or obuk
                ubot[ij] = 0.;// ufluxbot[ij] / (0.5*(ustar[ij-ii]*fm(zsl, z0m, obuk[ij-ii]) + ustar[ij]*fm(zsl, z0m, obuk[ij]))) + u[ijk];
                vbot[ij] = 0.;// vfluxbot[ij] / (0.5*(ustar[ij-jj]*fm(zsl, z0m, obuk[ij-jj]) + ustar[ij]*fm(zsl, z0m, obuk[ij]))) + v[ijk];
            }

            grid->boundary_cyclic_2d(ubot);
            grid->boundary_cyclic_2d(vbot);
            */
        }

        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                // use the linearly interpolated grad, rather than the MO grad,
                // to prevent giving unresolvable gradients to advection schemes
                // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0m*ustar[ij]) * phih(zsl/obuk[ij]);
                ugradbot[ij] = (u[ijk]-ubot[ij])/zsl;
                vgradbot[ij] = (v[ijk]-vbot[ij])/zsl;
            }
    }

    template<typename TF>
    void surfs(TF* restrict ustar, TF* restrict obuk, TF* restrict var,
               TF* restrict varbot, TF* restrict vargradbot, TF* restrict varfluxbot,
               const TF zsl, const TF z0m, const TF z0h, const Boundary_type bcbot,
               const int istart, const int iend, const int jstart, const int jend, const int kstart,
               const int icells, const int jcells, const int kk,
               Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int jj = icells;

        // the surface value is known, calculate the flux and gradient
        if (bcbot == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    varfluxbot[ij] = -(var[ijk]-varbot[ij])*ustar[ij]*most::fh(zsl, z0h, obuk[ij]);
                    // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0h*ustar[ij]) * phih(zsl/obuk[ij]);
                    // use the linearly interpolated grad, rather than the MO grad,
                    // to prevent giving unresolvable gradients to advection schemes
                    vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
                }
        }
        else if (bcbot == Boundary_type::Flux_type)
        {
            // the flux is known, calculate the surface value and gradient
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    varbot[ij] = varfluxbot[ij] / (ustar[ij]*most::fh(zsl, z0h, obuk[ij])) + var[ijk];
                    // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0h*ustar[ij]) * phih(zsl/obuk[ij]);
                    // use the linearly interpolated grad, rather than the MO grad,
                    // to prevent giving unresolvable gradients to advection schemes
                    vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
                }
        }
    }
}

template<typename TF>
Boundary_surface<TF>::Boundary_surface(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Boundary<TF>(masterin, gridin, fieldsin, inputin)
{
    swboundary = "surface";

    #ifdef USECUDA
    ustar_g = 0;
    obuk_g  = 0;
    nobuk_g = 0;
    zL_sl_g = 0;
    f_sl_g  = 0;
    #endif
}

template<typename TF>
Boundary_surface<TF>::~Boundary_surface()
{
    #ifdef USECUDA
    clear_device();
    #endif
}

template<typename TF>
void Boundary_surface<TF>::create(Input& input, Netcdf_handle& input_nc, Stats<TF>& stats)
{
    Boundary<TF>::process_time_dependent(input, input_nc);

    // add variables to the statistics
    if (stats.get_switch())
    {
        stats.add_time_series("ustar", "Surface friction velocity", "m s-1");
        stats.add_time_series("obuk", "Obukhov length", "m");
    }
}

template<typename TF>
void Boundary_surface<TF>::init(Input& inputin, Thermo<TF>& thermo)
{
    // 1. Process the boundary conditions now all fields are registered.
    process_bcs(inputin);

    // 2. Read and check the boundary_surface specific settings.
    process_input(inputin, thermo);

    // 3. Allocate and initialize the 2D surface fields.
    init_surface();

    // 4. Initialize the boundary cyclic.
    boundary_cyclic.init();
}

template<typename TF>
void Boundary_surface<TF>::process_input(Input& inputin, Thermo<TF>& thermo)
{
    int nerror = 0;

    z0m = inputin.get_item<TF>("boundary", "z0m", "");
    z0h = inputin.get_item<TF>("boundary", "z0h", "");

    // crash in case fixed gradient is prescribed
    if (mbcbot == Boundary_type::Neumann_type)
    {
        master.print_error("Neumann bc is not supported in surface model\n");
        ++nerror;
    }
    // read the ustar value only if fixed fluxes are prescribed
    else if (mbcbot == Boundary_type::Ustar_type)
        ustarin = inputin.get_item<TF>("boundary", "ustar", "");

    // process the scalars
    for (auto& it : sbc)
    {
        // crash in case fixed gradient is prescribed
        if (it.second.bcbot == Boundary_type::Neumann_type)
        {
            master.print_error("fixed gradient bc is not supported in surface model\n");
            ++nerror;
        }

        // crash in case of fixed momentum flux and dirichlet bc for scalar
        if (it.second.bcbot == Boundary_type::Dirichlet_type && mbcbot == Boundary_type::Ustar_type)
        {
            master.print_error("fixed Ustar bc in combination with Dirichlet bc for scalars is not supported\n");
            ++nerror;
        }
    }

    // check whether the prognostic thermo vars are of the same type
    std::vector<std::string> thermolist;
    thermo.get_prog_vars(thermolist);

    auto it = thermolist.begin();

    // save the bc of the first thermo field in case thermo is enabled
    if (it != thermolist.end())
        thermobc = sbc[*it].bcbot;
    else
        // Set the thermobc to Flux_type to avoid ininitialized errors.
        thermobc = Boundary_type::Flux_type;

    while (it != thermolist.end())
    {
        if (sbc[*it].bcbot != thermobc)
        {
            ++nerror;
            master.print_error("all thermo variables need to have the same bc type\n");
        }
        ++it;
    }

    if (nerror)
        throw 1;

    /*
    // Cross sections
    allowedcrossvars.push_back("ustar");
    allowedcrossvars.push_back("obuk");

    // Read list of cross sections
    nerror += inputin.get_list(&crosslist, "boundary", "crosslist" , "");

    // Get global cross-list from cross.cxx
    std::vector<std::string> *crosslist_global = model->cross->get_crosslist();

    // Check input list of cross variables (crosslist)
    std::vector<std::string>::iterator it2=crosslist_global->begin();
    while (it2 != crosslist_global->end())
    {
        if (std::count(allowedcrossvars.begin(),allowedcrossvars.end(),*it2))
        {
            // Remove variable from global list, put in local list
            crosslist.push_back(*it2);
            crosslist_global->erase(it2); // erase() returns iterator of next element..
        }
        else
            ++it2;
    }
    */
}

template<typename TF>
void Boundary_surface<TF>::init_surface()
{
    auto& gd = grid.get_grid_data();

    obuk.resize(gd.ijcells);
    nobuk.resize(gd.ijcells);
    ustar.resize(gd.ijcells);

    const int jj = gd.icells;

    // Initialize the obukhov length on a small number.
    for (int j=0; j<gd.jcells; ++j)
        #pragma ivdep
        for (int i=0; i<gd.icells; ++i)
        {
            const int ij = i + j*jj;
            obuk[ij]  = Constants::dsmall;
            nobuk[ij] = 0;
        }
}

/*
template<typename TF>
void Boundary_surface<TF>::exec_cross(int iotime)
{
    int nerror = 0;

    for (std::vector<std::string>::const_iterator it=crosslist.begin(); it<crosslist.end(); ++it)
    {
        if (*it == "ustar")
            nerror += model->cross->cross_plane(ustar, fields->atmp["tmp1"]->data, "ustar",iotime);
        else if (*it == "obuk")
            nerror += model->cross->cross_plane(obuk,  fields->atmp["tmp1"]->data, "obuk",iotime);
    }

    if (nerror)
        throw 1;
}
*/

template<typename TF>
void Boundary_surface<TF>::exec_stats(Stats<TF>& stats)
{
    const TF no_offset = 0.;
    stats.calc_stats_2d("obuk", obuk, no_offset, {"mean"});
    stats.calc_stats_2d("ustar", ustar, no_offset, {"mean"});
}

template<typename TF>
void Boundary_surface<TF>::set_values()
{
    auto& gd = grid.get_grid_data();

    // Call the base class function.
    Boundary<TF>::set_values();

    // Override the boundary settings in order to enforce dirichlet BC for surface model.
    set_bc<TF>(fields.mp.at("u")->fld_bot.data(), fields.mp.at("u")->grad_bot.data(), fields.mp.at("u")->flux_bot.data(),
           Boundary_type::Dirichlet_type, ubot, fields.visc, grid.utrans,
           gd.icells, gd.jcells);
    set_bc<TF>(fields.mp.at("v")->fld_bot.data(), fields.mp.at("v")->grad_bot.data(), fields.mp.at("v")->flux_bot.data(),
           Boundary_type::Dirichlet_type, vbot, fields.visc, grid.vtrans,
           gd.icells, gd.jcells);

    // in case the momentum has a fixed ustar, set the value to that of the input
    if (mbcbot == Boundary_type::Ustar_type)
        set_ustar();

    // Prepare the lookup table for the surface solver
    init_solver();
}

template<typename TF>
void Boundary_surface<TF>::set_ustar()
{
    auto& gd = grid.get_grid_data();
    const int jj = gd.icells;

    set_bc<TF>(fields.mp.at("u")->fld_bot.data(), fields.mp.at("u")->grad_bot.data(), fields.mp.at("u")->flux_bot.data(),
           mbcbot, ubot, fields.visc, grid.utrans,
           gd.icells, gd.jcells);
    set_bc<TF>(fields.mp.at("v")->fld_bot.data(), fields.mp.at("v")->grad_bot.data(), fields.mp.at("v")->flux_bot.data(),
           mbcbot, vbot, fields.visc, grid.vtrans,
           gd.icells, gd.jcells);

    for (int j=0; j<gd.jcells; ++j)
        #pragma ivdep
        for (int i=0; i<gd.icells; ++i)
        {
            const int ij = i + j*jj;
            // Limit ustar at 1e-4 to avoid zero divisions.
            ustar[ij] = std::max(static_cast<TF>(0.0001), ustarin);
        }
}

// Prepare the surface layer solver.
template<typename TF>
void Boundary_surface<TF>::init_solver()
{
    auto& gd = grid.get_grid_data();

    zL_sl.resize(nzL);
    f_sl.resize(nzL);

    std::vector<TF> zL_tmp(nzL);

    // Calculate the non-streched part between -5 to 10 z/L with 9/10 of the points,
    // and stretch up to -1e4 in the negative limit.
    // Alter next three values in case the range need to be changed.
    const TF zL_min = -1.e4;
    const TF zLrange_min = -5.;
    const TF zLrange_max = 10.;

    TF dzL = (zLrange_max - zLrange_min) / (9.*nzL/10.-1.);
    zL_tmp[0] = -zLrange_max;
    for (int n=1; n<9*nzL/10; ++n)
        zL_tmp[n] = zL_tmp[n-1] + dzL;

    // Stretch the remainder of the z/L values far down for free convection.
    const TF zLend = -(zL_min - zLrange_min);

    // Find stretching that ends up at the correct value using geometric progression.
    TF r  = 1.01;
    TF r0 = Constants::dhuge;
    while (std::abs( (r-r0)/r0 ) > 1.e-10)
    {
        r0 = r;
        r  = std::pow( 1. - (zLend/dzL)*(1.-r), (1./ (nzL/10.) ) );
    }

    for (int n=9*nzL/10; n<nzL; ++n)
    {
        zL_tmp[n] = zL_tmp[n-1] + dzL;
        dzL *= r;
    }

    // Calculate the final array and delete the temporary array.
    for (int n=0; n<nzL; ++n)
        zL_sl[n] = -zL_tmp[nzL-n-1];

    // Calculate the evaluation function.
    if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Flux_type)
    {
        const TF zsl = gd.z[gd.kstart];
        for (int n=0; n<nzL; ++n)
            f_sl[n] = zL_sl[n] * std::pow(most::fm(zsl, z0m, zsl/zL_sl[n]), 3);
    }
    else if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Dirichlet_type)
    {
        const TF zsl = gd.z[gd.kstart];
        for (int n=0; n<nzL; ++n)
            f_sl[n] = zL_sl[n] * std::pow(most::fm(zsl, z0m, zsl/zL_sl[n]), 2) / most::fh(zsl, z0h, zsl/zL_sl[n]);
    }
}

#ifndef USECUDA
template<typename TF>
void Boundary_surface<TF>::update_bcs(Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();

    // Start with retrieving the stability information.
    if (thermo.get_switch() == "0")
    {
        auto dutot = fields.get_tmp();
        stability_neutral(ustar.data(), obuk.data(),
                          fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(),
                          fields.mp.at("u")->fld_bot.data(), fields.mp.at("v")->fld_bot.data(),
                          dutot->fld.data(), gd.z.data(),
                          z0m,
                          gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart,
                          gd.icells, gd.jcells, gd.ijcells,
                          mbcbot, boundary_cyclic);
        fields.release_tmp(dutot);
    }
    else
    {
        auto buoy = fields.get_tmp();
        auto tmp = fields.get_tmp();

        thermo.get_buoyancy_surf(*buoy, false);
        stability(ustar.data(), obuk.data(), buoy->flux_bot.data(),
                  fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), buoy->fld.data(),
                  fields.mp.at("u")->fld_bot.data(), fields.mp.at("v")->fld_bot.data(), buoy->fld_bot.data(),
                  tmp->fld.data(), gd.z.data(),
                  zL_sl.data(), f_sl.data(), nobuk.data(),
                  z0m, z0h,
                  gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart,
                  gd.icells, gd.jcells, gd.ijcells,
                  mbcbot, thermobc, boundary_cyclic);

        fields.release_tmp(buoy);
        fields.release_tmp(tmp);
    }

    // Calculate the surface value, gradient and flux depending on the chosen boundary condition.
    surfm(ustar.data(), obuk.data(),
          fields.mp.at("u")->fld.data(), fields.mp.at("u")->fld_bot.data(), fields.mp.at("u")->grad_bot.data(), fields.mp.at("u")->flux_bot.data(),
          fields.mp.at("v")->fld.data(), fields.mp.at("v")->fld_bot.data(), fields.mp.at("v")->grad_bot.data(), fields.mp.at("v")->flux_bot.data(),
          gd.z[gd.kstart], z0m, mbcbot,
          gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart,
          gd.icells, gd.jcells, gd.ijcells,
          boundary_cyclic);

    for (auto& it : fields.sp)
    {
        surfs(ustar.data(), obuk.data(), it.second->fld.data(),
              it.second->fld_bot.data(), it.second->grad_bot.data(), it.second->flux_bot.data(),
              gd.z[gd.kstart], z0m, z0h, sbc[it.first].bcbot,
              gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart,
              gd.icells, gd.jcells, gd.ijcells,
              boundary_cyclic);
    }
}
#endif

template<typename TF>
void Boundary_surface<TF>::update_slave_bcs()
{
    // This function does nothing when the surface model is enabled, because
    // the fields are computed by the surface model in update_bcs.
}

template class Boundary_surface<double>;
template class Boundary_surface<float>;
