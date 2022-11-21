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
#include <algorithm>
// #include <math.h>
#include <iostream>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "field3d_operators.h"
#include "timedep.h"
#include "force.h"
#include "defines.h"
#include "finite_difference.h"
#include "timeloop.h"
#include "boundary.h"
#include "netcdf_interface.h"
#include "stats.h"
#include "thermo.h"

using namespace Finite_difference::O2;

namespace
{
    template<typename TF>
    void add_pressure_force(
            TF* restrict ut, const TF fbody,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    ut[ijk] += fbody;
                }
    }


    template<typename TF>
    void enforce_fixed_flux(
            TF* restrict ut,
            const TF u_flux, const TF u_mean, const TF ut_mean, const TF u_grid,
            const TF dt,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const TF fbody = (u_flux - u_mean - u_grid) / dt - ut_mean;

        add_pressure_force(ut, fbody, istart, iend, jstart, jend, kstart, kend, jj, kk);
    }


    template<typename TF>
    void calc_coriolis_2nd(
            TF* const restrict ut, TF* const restrict vt,
            const TF* const restrict u , const TF* const restrict v ,
            const TF* const restrict ug, const TF* const restrict vg, TF const fc,
            const TF ugrid, const TF vgrid,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    ut[ijk] += fc * (TF(0.25)*(v[ijk-ii] + v[ijk] + v[ijk-ii+jj] + v[ijk+jj]) + vgrid - vg[k]);
                }

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    vt[ijk] -= fc * (TF(0.25)*(u[ijk-jj] + u[ijk] + u[ijk+ii-jj] + u[ijk+ii]) + ugrid - ug[k]);
                }
    }

    template<typename TF>
    void calc_coriolis_4th(
            TF* const restrict ut, TF* const restrict vt,
            const TF* const restrict u , const TF* const restrict v ,
            const TF* const restrict ug, const TF* const restrict vg, TF const fc,
            const TF ugrid, const TF vgrid,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using namespace Finite_difference::O4;

        const int ii1 = 1;
        const int ii2 = 2;
        const int jj1 = 1*icells;
        const int jj2 = 2*icells;
        const int kk1 = 1*ijcells;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    ut[ijk] += fc * ( ( ci0<TF>*(ci0<TF>*v[ijk-ii2-jj1] + ci1<TF>*v[ijk-ii1-jj1] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk+ii1-jj1])
                                      + ci1<TF>*(ci0<TF>*v[ijk-ii2    ] + ci1<TF>*v[ijk-ii1    ] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+ii1    ])
                                      + ci2<TF>*(ci0<TF>*v[ijk-ii2+jj1] + ci1<TF>*v[ijk-ii1+jj1] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+ii1+jj1])
                                      + ci3<TF>*(ci0<TF>*v[ijk-ii2+jj2] + ci1<TF>*v[ijk-ii1+jj2] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+ii1+jj2]) )
                                    + vgrid - vg[k] );
                }

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    vt[ijk] -= fc * ( ( ci0<TF>*(ci0<TF>*u[ijk-ii1-jj2] + ci1<TF>*u[ijk-jj2] + ci2<TF>*u[ijk+ii1-jj2] + ci3<TF>*u[ijk+ii2-jj2])
                                      + ci1<TF>*(ci0<TF>*u[ijk-ii1-jj1] + ci1<TF>*u[ijk-jj1] + ci2<TF>*u[ijk+ii1-jj1] + ci3<TF>*u[ijk+ii2-jj1])
                                      + ci2<TF>*(ci0<TF>*u[ijk-ii1    ] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+ii1    ] + ci3<TF>*u[ijk+ii2    ])
                                      + ci3<TF>*(ci0<TF>*u[ijk-ii1+jj1] + ci1<TF>*u[ijk+jj1] + ci2<TF>*u[ijk+ii1+jj1] + ci3<TF>*u[ijk+ii2+jj1]) )
                                    + ugrid - ug[k]);
                }
    }

    template<typename TF>
    void calc_large_scale_source(
            TF* const restrict st, const TF* const restrict sls,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] += sls[k];
                }
    }

    template<typename TF>
    void calc_nudging_tendency(
            TF* const restrict fldtend, const TF* const restrict fldmean,
            const TF* const restrict ref, const TF* const restrict factor,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        for (int k=kstart; k<kend; ++k)
        {
            const TF tend = -factor[k] * (fldmean[k] - ref[k]);
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    fldtend[ijk] += tend;
                }
        }
    }

    template<typename TF>
    void rescale_nudgeprof(TF* const restrict fldmean, const int kinv, const int kstart, const int kend)
    {
        for (int k=kstart+1; k<kinv; ++k)
            fldmean[k] = fldmean[kstart];

        for (int k=kinv+1; k<kend-2; ++k)
            fldmean[k] = fldmean[kend-1];
    }

    template<typename TF>
    void advec_wls_2nd_mean(
            TF* const restrict st, const TF* const restrict s,
            const TF* const restrict wls, const TF* const dzhi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        // use an upwind differentiation
        for (int k=kstart; k<kend; ++k)
        {
            if (wls[k] > 0.)
            {
                for (int j=jstart; j<jend; ++j)
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        st[ijk] -=  wls[k] * (s[k]-s[k-1])*dzhi[k];
                    }
            }
            else
            {
                for (int j=jstart; j<jend; ++j)
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        st[ijk] -=  wls[k] * (s[k+1]-s[k])*dzhi[k+1];
                    }
            }
        }
    }

    template<typename TF>
    void advec_wls_2nd_local(
            TF* const restrict st, const TF* const restrict s,
            const TF* const restrict wls, const TF* const dzhi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        // use an upwind differentiation
        for (int k=kstart; k<kend; ++k)
        {
            if (wls[k] > 0.)
            {
                for (int j=jstart; j<jend; ++j)
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        st[ijk] -=  wls[k] * (s[ijk]-s[ijk-kk])*dzhi[k];
                    }
            }
            else
            {
                for (int j=jstart; j<jend; ++j)
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        st[ijk] -=  wls[k] * (s[ijk+kk]-s[ijk])*dzhi[k+1];
                    }
            }
        }
    }

    // Especially, check this one thoroughly (indices, grid positions, etc.)! - SvdLinden, 28.04.21
    template<typename TF>
    void advec_wls_2nd_local_w(
            TF* const restrict st, const TF* const restrict s,
            const TF* const restrict wls, const TF* const dzi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        // should not be needed to do kend-1 separately ? CHECK?
        // use an upwind differentiation
        for (int k=kstart+1; k<kend; ++k) // for (int k=kstart+2; k<kend-1; ++k)
        {
            if ( interp2( wls[k-1], wls[k] ) > 0.) // formeel ook in conditie interp2
            {
                for (int j=jstart; j<jend; ++j)
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        st[ijk] -=  interp2( wls[k-1], wls[k] ) * (s[ijk]-s[ijk-kk])*dzi[k-1]; // HIER DUS dz !! maar waar begint dz[kstart] +1 of niet??
                    }
            }
            else
            {
                for (int j=jstart; j<jend; ++j)
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        st[ijk] -=  interp2( wls[k-1], wls[k] ) * (s[ijk+kk]-s[ijk])*dzi[k];
                    }
            }
        }
    }

    template<typename TF>
    void add_offset(TF* const restrict prof, const TF offset, const int kstart, const int kend)
    {
        for (int k=kstart; k<kend; ++k)
            prof[k] += offset;
    }
}

template<typename TF>
Force<TF>::Force(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), field3d_operators(masterin, gridin, fieldsin)
{
    std::string swlspres_in = inputin.get_item<std::string>("force", "swlspres", "", "0");
    std::string swwls_in    = inputin.get_item<std::string>("force", "swwls"   , "", "0");
    std::string swnudge_in  = inputin.get_item<std::string>("force", "swnudge" , "", "0");
    std::string swls_in     = inputin.get_item<std::string>("force", "swls"    , "", "0");

    if (swwls_in == "1" || swwls_in == "mean" || swwls_in == "local")
        bool swwls_mom = inputin.get_item<bool>("force", "swwls_mom", "", false);

    // Checks on input:
    if (swwls_in == "1")
    {
        master.print_warning("\"swwls=1\" has been replaced by \"swwls=mean\" or \"swwls=local\". Defaulting to \"swwls=mean\"\n");
        swwls_in = "mean";
    }

    // Set the internal switches and read other required input
    // Large-scale pressure forcing.
    if (swlspres_in == "0")
    {
        swlspres = Large_scale_pressure_type::Disabled;
    }
    else if (swlspres_in == "uflux")
    {
        swlspres = Large_scale_pressure_type::Fixed_flux;
        uflux = inputin.get_item<TF>("force", "uflux", "");
    }
    else if (swlspres_in == "dpdx")
    {
        swlspres = Large_scale_pressure_type::Pressure_gradient;
        dpdx = inputin.get_item<TF>("force", "dpdx", "");
    }
    else if (swlspres_in == "geo")
    {
        swlspres = Large_scale_pressure_type::Geo_wind;
        fc = inputin.get_item<TF>("force", "fc", "");
        tdep_geo.emplace("u_geo", new Timedep<TF>(master, grid, "u_geo", inputin.get_item<bool>("force", "swtimedep_geo", "", false)));
        tdep_geo.emplace("v_geo", new Timedep<TF>(master, grid, "v_geo", inputin.get_item<bool>("force", "swtimedep_geo", "", false)));
    }
    else
    {
        throw std::runtime_error("Invalid option for \"swlspres\"");
    }

    // Large-scale tendencies due to advection and other processes.
    if (swls_in == "0")
        swls = Large_scale_tendency_type::Disabled;
    else if (swls_in == "1")
    {
        swls = Large_scale_tendency_type::Enabled;
        lslist = inputin.get_list<std::string>("force", "lslist", "", std::vector<std::string>());

        if (inputin.get_item<bool>("force", "swtimedep_ls", "", false))
        {
            std::vector<std::string> tdepvars = inputin.get_list<std::string>("force", "timedeplist_ls", "", std::vector<std::string>());
            for (auto& it : tdepvars)
                tdep_ls.emplace(it, new Timedep<TF>(master, grid, it+"_ls", true));
        }
    }
    else
    {
        throw std::runtime_error("Invalid option for \"swls\"");
    }

    // Large-scale subsidence.
    if (swwls_in == "0")
        swwls = Large_scale_subsidence_type::Disabled;
    else if (swwls_in == "mean")
    {

        swwls = Large_scale_subsidence_type::Mean_field;
        fields.set_calc_mean_profs(true);
    }
    else if (swwls_in == "local")
    {
        swwls = Large_scale_subsidence_type::Local_field;
        fields.set_calc_mean_profs(true);
    }
    else
    {
        throw std::runtime_error("Invalid option for \"swwls\"");
    }
    tdep_wls = std::make_unique<Timedep<TF>>(master, grid, "w_ls", inputin.get_item<bool>("force", "swtimedep_wls", "", false));

    // Nudging.
    if (swnudge_in == "0")
        swnudge = Nudging_type::Disabled;
    else if (swnudge_in == "1")
    {
        swnudge = Nudging_type::Enabled;
        nudgelist       = inputin.get_list<std::string>("force", "nudgelist", "", std::vector<std::string>());
        scalednudgelist = inputin.get_list<std::string>("force", "scalednudgelist", "", std::vector<std::string>());

        if (inputin.get_item<bool>("force", "swtimedep_nudge", "", false))
        {
            std::vector<std::string> tdepvars = inputin.get_list<std::string>("force", "timedeplist_nudge", "", std::vector<std::string>());
            for(auto& it : tdepvars)
                tdep_nudge.emplace(it, new Timedep<TF>(master, grid, it+"_nudge", true));
        }
        fields.set_calc_mean_profs(true);
    }
    else
    {
        throw std::runtime_error("Invalid option for \"swnduge\"");
    }
}

template <typename TF>
Force<TF>::~Force()
{
}

template <typename TF>
void Force<TF>::init()
{
    auto& gd = grid.get_grid_data();

    if (swlspres == Large_scale_pressure_type::Geo_wind)
    {
        ug.resize(gd.kcells);
        vg.resize(gd.kcells);
    }

    if (swls == Large_scale_tendency_type::Enabled)
    {
        for (auto& it : lslist)
            lsprofs[it] = std::vector<TF>(gd.kcells);
    }
    if (swwls == Large_scale_subsidence_type::Mean_field || swwls == Large_scale_subsidence_type::Local_field)
        wls.resize(gd.kcells);

    if (swnudge == Nudging_type::Enabled)
    {
        nudge_factor.resize(gd.kcells);
        for (auto& it : nudgelist)
            nudgeprofs[it] = std::vector<TF>(gd.kcells);

    }
}

template <typename TF>
void Force<TF>::create(Input& inputin, Netcdf_handle& input_nc, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();
    Netcdf_group& group_nc = input_nc.get_group("init");

    std::string timedep_dim = "time_ls";

    if (swlspres == Large_scale_pressure_type::Fixed_flux)
    {
        stats.add_tendency(*fields.mt.at("u"), "z", tend_name_pres, tend_longname_pres);
    }
    else if (swlspres == Large_scale_pressure_type::Geo_wind)
    {
        group_nc.get_variable(ug, "u_geo", {0}, {gd.ktot});
        group_nc.get_variable(vg, "v_geo", {0}, {gd.ktot});
        std::rotate(ug.rbegin(), ug.rbegin() + gd.kstart, ug.rend());
        std::rotate(vg.rbegin(), vg.rbegin() + gd.kstart, vg.rend());

        const TF offset = 0;
        for (auto& it : tdep_geo)
            it.second->create_timedep_prof(input_nc, offset, timedep_dim);
        stats.add_tendency(*fields.mt.at("u"), "z", tend_name_cor, tend_longname_cor);
        stats.add_tendency(*fields.mt.at("v"), "z", tend_name_cor, tend_longname_cor);

    }

    if (swls == Large_scale_tendency_type::Enabled)
    {
        // Check whether the fields in the list exist in the prognostic fields.
        for (std::string& it : lslist)
            if (!fields.ap.count(it))
            {
                std::string msg = "field " + it + " in [force][lslist] is illegal";
                throw std::runtime_error(msg);
            }

        // Read the large scale sources, which are the variable names with a "_ls" suffix.
        for (std::string& it : lslist)
        {
            group_nc.get_variable(lsprofs[it], it+"_ls", {0}, {gd.ktot});
            std::rotate(lsprofs[it].rbegin(), lsprofs[it].rbegin() + gd.kstart, lsprofs[it].rend());
        }

        // Process the time dependent data
        const TF offset = 0;
        for (auto& it : tdep_ls)
            it.second->create_timedep_prof(input_nc, offset, timedep_dim);

        for (std::string& it : lslist)
            stats.add_tendency(*fields.at.at(it), "z", tend_name_ls, tend_longname_ls);
    }

    if (swnudge == Nudging_type::Enabled)
    {
        // Get profile with nudging factor as function of height
        group_nc.get_variable(nudge_factor, "nudgefac", {0}, {gd.ktot});
        std::rotate(nudge_factor.rbegin(), nudge_factor.rbegin() + gd.kstart, nudge_factor.rend());

        // check whether the fields in the list exist in the prognostic fields
        for (auto& it : nudgelist)
            if (!fields.ap.count(it))
            {
                std::string msg = "field " + it + " in [force][nudgelist] is illegal";
                throw std::runtime_error(msg);
            }

        // Read the nudging profiles, which are the variable names with a "nudge" suffix
        for (auto& it : nudgelist)
        {
            group_nc.get_variable(nudgeprofs[it], it+"_nudge", {0}, {gd.ktot});
            std::rotate(nudgeprofs[it].rbegin(), nudgeprofs[it].rbegin() + gd.kstart, nudgeprofs[it].rend());

            // Account for the Galilean transformation
            if (it == "u")
                add_offset(nudgeprofs[it].data(), -grid.utrans, gd.kstart, gd.kend);
            else if (it == "v")
                add_offset(nudgeprofs[it].data(), -grid.vtrans, gd.kstart, gd.kend);
        }

        // Process the time dependent data
        for (auto& it : tdep_nudge)
        {
            // Account for the Galilean transformation
            TF offset;
            if (it.first == "u")
                offset = -grid.utrans;
            else if (it.first == "v")
                offset = -grid.vtrans;
            else
                offset = 0;

            it.second->create_timedep_prof(input_nc, offset, timedep_dim);
            stats.add_tendency(*fields.at.at(it.first), "z", tend_name_nudge, tend_longname_nudge);
        }
    }

    // Get the large scale vertical velocity from the input
    if (swwls == Large_scale_subsidence_type::Mean_field || swwls == Large_scale_subsidence_type::Local_field)
    {
        group_nc.get_variable(wls, "w_ls", {0}, {gd.ktot});
        std::rotate(wls.rbegin(), wls.rbegin() + gd.kstart, wls.rend());

        const TF offset = 0;
        tdep_wls->create_timedep_prof(input_nc, offset, timedep_dim);

        for (auto& it : fields.st)
            stats.add_tendency(*it.second, "z", tend_name_subs, tend_longname_subs);

        if (swwls_mom)
        {
            // Initialize statistics output also for u,v,w.
            stats.add_tendency(*fields.mt.at("u"), "z", tend_name_subs, tend_longname_subs);
            stats.add_tendency(*fields.mt.at("v"), "z", tend_name_subs, tend_longname_subs);

            if (swwls == Large_scale_subsidence_type::Local_field)
                stats.add_tendency(*fields.mt.at("w"), "zh", tend_name_subs, tend_longname_subs);
        }
    }
}

#ifndef USECUDA
template <typename TF>
void Force<TF>::exec(double dt, Thermo<TF>& thermo, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    if (swlspres == Large_scale_pressure_type::Fixed_flux)
    {
        const TF u_mean  = field3d_operators.calc_mean(fields.ap.at("u")->fld.data());
        const TF ut_mean = field3d_operators.calc_mean(fields.at.at("u")->fld.data());

        enforce_fixed_flux<TF>(
                fields.at.at("u")->fld.data(), uflux, u_mean, ut_mean, grid.utrans, dt,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);

        stats.calc_tend(*fields.mt.at("u"), tend_name_pres);
    }

    else if (swlspres == Large_scale_pressure_type::Pressure_gradient)
    {
        const TF fbody = TF(-1.)*dpdx;

        add_pressure_force<TF>(
                fields.at.at("u")->fld.data(), fbody,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);

        stats.calc_tend(*fields.mt.at("u"), tend_name_pres);
    }

    else if (swlspres == Large_scale_pressure_type::Geo_wind)
    {
        if (grid.get_spatial_order() == Grid_order::Second)
            calc_coriolis_2nd<TF>(fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), ug.data(), vg.data(), fc,
            grid.utrans, grid.vtrans,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

        else if (grid.get_spatial_order() == Grid_order::Fourth)
            calc_coriolis_4th<TF>(fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), ug.data(), vg.data(), fc,
            grid.utrans, grid.vtrans,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

        stats.calc_tend(*fields.mt.at("u"), tend_name_cor);
        stats.calc_tend(*fields.mt.at("v"), tend_name_cor);

    }

    if (swls == Large_scale_tendency_type::Enabled)
    {
        for (auto& it : lslist)
        {
            calc_large_scale_source<TF>(
                    fields.at.at(it)->fld.data(), lsprofs.at(it).data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

            stats.calc_tend(*fields.at.at(it), tend_name_ls);
        }

    }

    if (swwls == Large_scale_subsidence_type::Mean_field )
    {
        if (swwls_mom)
        {
            // Also apply to the velocity components u,v - SvdLinden, 28.04.21
            advec_wls_2nd_mean<TF>(
                    fields.mt.at("u")->fld.data(), fields.mp.at("u")->fld_mean.data(), wls.data(), gd.dzhi.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
            stats.calc_tend(*fields.mt.at("u"), tend_name_subs);

            advec_wls_2nd_mean<TF>(
                    fields.mt.at("v")->fld.data(), fields.mp.at("v")->fld_mean.data(), wls.data(), gd.dzhi.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
            stats.calc_tend(*fields.mt.at("v"), tend_name_subs);
        }

        for (auto& it : fields.st)
        {
            advec_wls_2nd_mean<TF>(
                    fields.st.at(it.first)->fld.data(), fields.sp.at(it.first)->fld_mean.data(), wls.data(), gd.dzhi.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
            stats.calc_tend(*it.second, tend_name_subs);
        }
    }
    else if ( swwls == Large_scale_subsidence_type::Local_field )
    {
        // New functions for the local subsidence term - SvdLinden, 28.04.21
        if (swwls_mom)
        {
            // Apply to all prognostic scalars, also velocity. Treat w-velocity separately
            advec_wls_2nd_local<TF>(
                    fields.mt.at("u")->fld.data(), fields.mp.at("u")->fld.data(), wls.data(), gd.dzhi.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
            stats.calc_tend(*fields.mt.at("u"), tend_name_subs);

            advec_wls_2nd_local<TF>(
                    fields.mt.at("v")->fld.data(), fields.mp.at("v")->fld.data(), wls.data(), gd.dzhi.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
            stats.calc_tend(*fields.mt.at("v"), tend_name_subs);

            advec_wls_2nd_local_w<TF>(
                    fields.mt.at("w")->fld.data(), fields.mp.at("w")->fld.data(), wls.data(), gd.dzi.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
            stats.calc_tend(*fields.mt.at("w"), tend_name_subs);
        }

        for (auto& it : fields.st)
        {
            advec_wls_2nd_local<TF>(
                    fields.st.at(it.first)->fld.data(), fields.sp.at(it.first)->fld.data(), wls.data(), gd.dzhi.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
            stats.calc_tend(*it.second, tend_name_subs);
        }
    }

    if (swnudge == Nudging_type::Enabled)
    {
        for (auto& it : nudgelist)
        {
            auto it1 = std::find(scalednudgelist.begin(), scalednudgelist.end(), it);
            if (it1 != scalednudgelist.end())
            {
                const int kinv = thermo.get_bl_depth();
                rescale_nudgeprof(nudgeprofs.at(it).data(), kinv, gd.kstart, gd.kend);
            }

            calc_nudging_tendency<TF>(
                    fields.at.at(it)->fld.data(), fields.ap.at(it)->fld_mean.data(),
                    nudgeprofs.at(it).data(), nudge_factor.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

            stats.calc_tend(*fields.at.at(it), tend_name_nudge);
        }
    }
}
#endif

#ifndef USECUDA
template <typename TF>
void Force<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (swls == Large_scale_tendency_type::Enabled)
    {
        for (auto& it : tdep_ls)
            it.second->update_time_dependent_prof(lsprofs.at(it.first), timeloop);
    }

    if (swnudge == Nudging_type::Enabled)
    {
        for (auto& it : tdep_nudge)
            it.second->update_time_dependent_prof(nudgeprofs.at(it.first), timeloop);
    }

    if (swlspres == Large_scale_pressure_type::Geo_wind)
    {
        tdep_geo.at("u_geo")->update_time_dependent_prof(ug, timeloop);
        tdep_geo.at("v_geo")->update_time_dependent_prof(vg, timeloop);
    }

    if (swwls == Large_scale_subsidence_type::Mean_field || swwls == Large_scale_subsidence_type::Local_field )
        tdep_wls->update_time_dependent_prof(wls, timeloop);

    // Idea: could decide to update interpolated wls to full levels here ? - SvdLinden, 28.04.21
}
#endif

template class Force<double>;
template class Force<float>;
