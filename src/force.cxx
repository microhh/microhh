/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
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
#include <math.h>
#include <algorithm>
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
#include "data_block.h"

using namespace Finite_difference::O2;

namespace
{
    template<typename TF>
    void enforce_fixed_flux(
            TF* restrict ut,
            const TF u_flux, const TF u_mean, const TF ut_mean, const TF u_grid,
            const TF dt,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const TF fbody = (u_flux - u_mean - u_grid) / dt - ut_mean;

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
    void rescale_nudgeprof(TF* const restrict fldmean, const int kinv, const int kstart, const int kend)
    {
        for (int k=kstart+1; k<kinv; ++k)
            fldmean[k] = fldmean[kstart];

        for (int k=kinv+1; k<kend-2; ++k)
            fldmean[k] = fldmean[kend-1];
    }

    template<typename TF>
    void advec_wls_2nd(
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
}

template<typename TF>
Force<TF>::Force(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), field3d_operators(masterin, gridin, fieldsin)
{
    std::string swlspres_in = inputin.get_item<std::string>("force", "swlspres", "", "0");
    std::string swls_in     = inputin.get_item<std::string>("force", "swls", "", "0");
    std::string swwls_in    = inputin.get_item<std::string>("force", "swwls", "", "0");
    std::string swnudge_in  = inputin.get_item<std::string>("force", "swnudge", "", "0");

    // Set the internal switches and read other required input

    // Large-scale pressure forcing.
    if (swlspres_in == "0")
    {
        swlspres = Large_scale_pressure_type::disabled;
    }
    else if (swlspres_in == "uflux")
    {
        swlspres = Large_scale_pressure_type::fixed_flux;
        uflux = inputin.get_item<TF>("force", "uflux", "");
    }
    else if (swlspres_in == "geo")
    {
        swlspres = Large_scale_pressure_type::geo_wind;
        fc = inputin.get_item<TF>("force", "fc", "");
        tdep_geo.emplace("ug", new Timedep<TF>(master, grid, "u_geo", inputin.get_item<bool>("force", "swtimedep_geo", "", false)));
        tdep_geo.emplace("vg", new Timedep<TF>(master, grid, "v_geo", inputin.get_item<bool>("force", "swtimedep_geo", "", false)));
    }
    else
    {
        throw std::runtime_error("Invalid option for \"swlspres\"");
    }

    // Large-scale tendencies due to advection and other processes.
    if (swls_in == "0")
        swls = Large_scale_tendency_type::disabled;
    else if (swls_in == "1")
    {
        swls = Large_scale_tendency_type::enabled;
        lslist = inputin.get_list<std::string>("force", "lslist", "", std::vector<std::string>());

        if (inputin.get_item<bool>("force", "swtimedep_ls", "", false))
        {
            std::vector<std::string> tdepvars = inputin.get_list<std::string>("force", "timedeplist_ls", "", std::vector<std::string>());
            for(auto& it : tdepvars)
                tdep_ls.emplace(it, new Timedep<TF>(master, grid, it+"_ls", true));
        }
    }
    else
    {
        throw std::runtime_error("Invalid option for \"swls\"");
    }

    // Large-scale subsidence.
    if (swwls_in == "0")
        swwls = Large_scale_subsidence_type::disabled;
    else if (swwls_in == "1")
    {
        swwls = Large_scale_subsidence_type::enabled;
        fields.set_calc_mean_profs(true);
    }
    else
    {
        throw std::runtime_error("Invalid option for \"swwls\"");
    }
    tdep_wls = std::make_unique<Timedep<TF>>(master, grid, "w_ls", inputin.get_item<bool>("force", "swtimedep_wls", "", false));

    // Nudging.
    if (swnudge_in == "0")
        swnudge = Nudging_type::disabled;
    else if (swnudge_in == "1")
    {
        swnudge = Nudging_type::enabled;
        nudgelist       = inputin.get_list<std::string>("force", "nudgelist", "", std::vector<std::string>());
        scalednudgelist = inputin.get_list<std::string>("force", "scalednudgelist", "", std::vector<std::string>());

        if (inputin.get_item<bool>("force", "swtimedep_nudge", "", false))
        {
            std::vector<std::string> tdepvars = inputin.get_list<std::string>("force", "timedeplist_nudge", "", std::vector<std::string>());
            for(auto& it : tdepvars)
                tdep_ls.emplace(it, new Timedep<TF>(master, grid, it+"_nudge", true));
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

    if (swlspres == Large_scale_pressure_type::geo_wind)
    {
        ug.resize(gd.kcells);
        vg.resize(gd.kcells);
    }

    if (swls == Large_scale_tendency_type::enabled)
    {
        for (auto& it : lslist)
            lsprofs[it] = std::vector<TF>(gd.kcells);
    }
    if (swwls == Large_scale_subsidence_type::enabled)
        wls.resize(gd.kcells);

    if (swnudge == Nudging_type::enabled)
    {
        nudge_factor.resize(gd.kcells);
        for (auto& it : nudgelist)
            nudgeprofs[it] = std::vector<TF>(gd.kcells);

    }
}

template <typename TF>
void Force<TF>::create(Input& inputin, Data_block& profs)
{
    auto& gd = grid.get_grid_data();
    if (swlspres == Large_scale_pressure_type::geo_wind)
    {
        profs.get_vector(ug, "ug", gd.kmax, 0, gd.kstart);
        profs.get_vector(vg, "vg", gd.kmax, 0, gd.kstart);

        for (auto& it : tdep_geo)
            it.second->create_timedep_prof();
    }

    if (swls == Large_scale_tendency_type::enabled)
    {
        // check whether the fields in the list exist in the prognostic fields
        for (auto & it : lslist)
            if (!fields.ap.count(it))
            {
                master.print_error("field %s in [force][lslist] is illegal\n", it.c_str());
            }

        // read the large scale sources, which are the variable names with a "ls" suffix
        for (auto & it : lslist)
            profs.get_vector(lsprofs[it],it+"ls", gd.kmax, 0, gd.kstart);

        // Process the time dependent data
        for (auto& it : tdep_ls)
            it.second->create_timedep_prof();
    }

    if (swnudge == Nudging_type::enabled)
    {
        // Get profile with nudging factor as function of height
        profs.get_vector(nudge_factor,"nudgefac", gd.kmax, 0, gd.kstart);

        // check whether the fields in the list exist in the prognostic fields
        for (auto & it : nudgelist)
            if (!fields.ap.count(it))
            {
                master.print_error("field %s in [force][nudgelist] is illegal\n", it.c_str());
            }

        // read the large scale sources, which are the variable names with a "nudge" suffix
        for (auto & it : nudgelist)
            profs.get_vector(nudgeprofs[it],it+"nudge", gd.kmax, 0, gd.kstart);

        // Process the time dependent data
        for (auto& it : tdep_nudge)
            it.second->create_timedep_prof();
    }

    // Get the large scale vertical velocity from the input
    if (swwls == Large_scale_subsidence_type::enabled)
    {
        profs.get_vector(wls,"wls", gd.kmax, 0, gd.kstart);
        tdep_wls->create_timedep_prof();
    }
}

#ifndef USECUDA
template <typename TF>
void Force<TF>::exec(double dt)
{
    auto& gd = grid.get_grid_data();

    if (swlspres == Large_scale_pressure_type::fixed_flux)
    {
        const TF u_mean  = field3d_operators.calc_mean(fields.ap.at("u")->fld.data());
        const TF ut_mean = field3d_operators.calc_mean(fields.at.at("u")->fld.data());

        enforce_fixed_flux<TF>(fields.at.at("u")->fld.data(), uflux, u_mean, ut_mean, grid.utrans, dt, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
    }

    else if (swlspres == Large_scale_pressure_type::geo_wind)
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
    }

    if (swls == Large_scale_tendency_type::enabled)
    {
        for (auto& it : lslist)
            calc_large_scale_source<TF>(
                    fields.st.at(it)->fld.data(), lsprofs.at(it).data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
    }

    if (swwls == Large_scale_subsidence_type::enabled)
    {
        for (auto& it : fields.st)
            advec_wls_2nd<TF>(
                    fields.st.at(it.first)->fld.data(), fields.sp.at(it.first)->fld_mean.data(), wls.data(), gd.dzhi.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
    }

    if (swnudge == Nudging_type::enabled)
    {
        for (auto& it : nudgelist)
        {
            auto it1 = std::find(scalednudgelist.begin(), scalednudgelist.end(), it);
            if (it1 != scalednudgelist.end())
            {
                const int kinv = calc_zi(fields.sp.at("thl")->fld_mean.data(), gd.kstart, gd.kend, 1);
                rescale_nudgeprof(nudgeprofs.at(it).data(), kinv, gd.kstart, gd.kend);
            }
            calc_nudging_tendency<TF>(
                    fields.st.at(it)->fld.data(), fields.sp.at(it)->fld_mean.data(),
                    nudgeprofs.at(it).data(), nudge_factor.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
        }
    }
}
#endif

#ifndef USECUDA
template <typename TF>
void Force<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (swls == Large_scale_tendency_type::enabled)
    {
        for (auto& it : tdep_ls)
            it.second->update_time_dependent_prof(lsprofs[it.first],timeloop);
    }

    if (swnudge == Nudging_type::enabled)
    {
        for (auto& it : tdep_nudge)
            it.second->update_time_dependent_prof(nudgeprofs[it.first],timeloop);
    }

    if (swlspres == Large_scale_pressure_type::geo_wind)
    {
        tdep_geo.at("ug")->update_time_dependent_prof(ug, timeloop);
        tdep_geo.at("vg")->update_time_dependent_prof(vg, timeloop);
    }

    if (swwls == Large_scale_subsidence_type::enabled)
        tdep_wls->update_time_dependent_prof(wls, timeloop);
}
#endif

template class Force<double>;
template class Force<float>;
