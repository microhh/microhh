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
#include <algorithm>
#include <iostream>
#include <math.h>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "field3d_operators.h"
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
    void enforce_fixed_flux(TF* restrict ut,
                            const TF u_flux, const TF u_mean, const TF ut_mean, const TF u_grid,
                            const TF dt, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int jj, const int kk)
    {
        const TF fbody = (u_flux - u_mean - u_grid) / dt - ut_mean;

        for (int k=kstart; k<kend; ++k)
           for (int j=jstart; j<jend; ++j)
               #pragma ivdep
               for (int i=istart; i<iend; ++i)
               {
                   const int ijk = i + j*jj + k*kk;
                   ut[ijk] +=fbody;
               }
    }
    template<typename TF>
    void calc_coriolis_2nd(
            TF* const restrict ut, TF* const restrict vt,
            const TF* const restrict u , const TF* const restrict v ,
            const TF* const restrict ug, const TF* const restrict vg, TF const fc,
            const TF ugrid, const TF vgrid, const int istart, const int iend, const int icells,
            const int jstart, const int jend, const int ijcells, const int kstart, const int kend)
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
    void calc_coriolis_4th(TF* const restrict ut, TF* const restrict vt,
                                  const TF* const restrict u , const TF* const restrict v ,
                                  const TF* const restrict ug, const TF* const restrict vg, TF const fc,
                                  const TF ugrid, const TF vgrid, const int istart, const int iend, const int icells,
                                  const int jstart, const int jend, const int ijcells, const int kstart, const int kend)
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
                    ut[ijk] += fc * ( ( ci0*(ci0*v[ijk-ii2-jj1] + ci1*v[ijk-ii1-jj1] + ci2*v[ijk-jj1] + ci3*v[ijk+ii1-jj1])
                                      + ci1*(ci0*v[ijk-ii2    ] + ci1*v[ijk-ii1    ] + ci2*v[ijk    ] + ci3*v[ijk+ii1    ])
                                      + ci2*(ci0*v[ijk-ii2+jj1] + ci1*v[ijk-ii1+jj1] + ci2*v[ijk+jj1] + ci3*v[ijk+ii1+jj1])
                                      + ci3*(ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii1+jj2] + ci2*v[ijk+jj2] + ci3*v[ijk+ii1+jj2]) )
                                    + vgrid - vg[k] );
                }

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    vt[ijk] -= fc * ( ( ci0*(ci0*u[ijk-ii1-jj2] + ci1*u[ijk-jj2] + ci2*u[ijk+ii1-jj2] + ci3*u[ijk+ii2-jj2])
                                      + ci1*(ci0*u[ijk-ii1-jj1] + ci1*u[ijk-jj1] + ci2*u[ijk+ii1-jj1] + ci3*u[ijk+ii2-jj1])
                                      + ci2*(ci0*u[ijk-ii1    ] + ci1*u[ijk    ] + ci2*u[ijk+ii1    ] + ci3*u[ijk+ii2    ])
                                      + ci3*(ci0*u[ijk-ii1+jj1] + ci1*u[ijk+jj1] + ci2*u[ijk+ii1+jj1] + ci3*u[ijk+ii2+jj1]) )
                                    + ugrid - ug[k]);
                }
    }

    template<typename TF>
    void calc_large_scale_source(
            TF* const restrict st, const TF* const restrict sls,
            const int istart, const int iend, const int icells, const int jstart, const int jend,
            const int ijcells, const int kstart, const int kend)
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
            const int istart, const int iend, const int icells, const int jstart, const int jend,
            const int ijcells, const int kstart, const int kend)
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
    void advec_wls_2nd(TF* const restrict st, const TF* const restrict s,
                              const TF* const restrict wls, const TF* const dzhi,
                              const int istart, const int iend, const int icells, const int jstart, const int jend,
                              const int ijcells, const int kstart, const int kend)
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
    // Read the switches from the input
    std::string swlspres_in = inputin.get_item<std::string>("force", "swlspres", "", "0");
    std::string swls_in     = inputin.get_item<std::string>("force", "swls", "", "0");
    std::string swwls_in    = inputin.get_item<std::string>("force", "swwls", "", "0");
    std::string swnudge_in  = inputin.get_item<std::string>("force", "swnudge", "", "0");

    // Set the internal switches and read other required input

    // Large-scale pressure forcing.
    if (swlspres_in == "0")
        swlspres = Large_scale_pressure_type::disabled;
    else if (swlspres_in == "uflux")
    {
        swlspres = Large_scale_pressure_type::fixed_flux;
        uflux = inputin.get_item<TF>("force", "uflux", "");
    }
    else if (swlspres_in == "geo")
    {
        swlspres = Large_scale_pressure_type::geo_wind;
        tdep_geo.sw = inputin.get_item<bool>("force", "swtimedep_geo",   "", "0");
        tdep_geo.vars = {"ug", "vg"};
    }
    else
        throw std::runtime_error("Invalid option for \"swlspres\"");

    // Large-scale tendencies due to advection and other processes.
    if (swls_in == "0")
        swls = Large_scale_tendency_type::disabled;
    else if (swls_in == "1")
    {
        swls = Large_scale_tendency_type::enabled;
        tdep_ls.sw = inputin.get_item<bool>("force", "swtimedep_ls",   "", "0");
        lslist = inputin.get_list<std::string>("force", "lslist", "", std::vector<std::string>());
        tdep_ls.vars = inputin.get_list<std::string>("force", "timedeptime_ls", "", std::vector<std::string>());
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
        tdep_wls.sw = inputin.get_item<bool>("force", "swtimedep_wls",   "", "0");
        fields.set_calc_mean_profs(true);
    }
    else
    {
        throw std::runtime_error("Invalid option for \"swwls\"");
    }

    // Nudging.
    if (swnudge_in == "0")
        swnudge = Nudging_type::disabled;
    else if (swnudge_in == "1")
    {
        swnudge = Nudging_type::enabled;
        tdep_nudge.sw   = inputin.get_item<bool>("force", "swtimedep_nudge",   "", "0");
        tdep_nudge.vars = inputin.get_list<std::string>("force", "timedeptime_nudge", "", std::vector<std::string>());
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
    #ifdef USECUDA
    clear_device();
    #endif
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
        for (auto& it : lsprofs)
            it.second.resize(gd.kcells);
    }
    if (swwls == Large_scale_subsidence_type::enabled)
        wls.resize(gd.kcells);

    if (swnudge == Nudging_type::enabled)
    {
        for (auto& it : nudgeprofs)
            it.second.resize(gd.kcells);
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

        if (tdep_geo.sw)
            create_timedep(tdep_geo,"g");
    }

    if (swls == Large_scale_tendency_type::enabled)
    {
        // check whether the fields in the list exist in the prognostic fields
        for (auto & it : lslist)
            if (!fields.ap.count(it))
            {
                throw std::runtime_error("field  %s in [force][lslist] is illegal\n");
            }

        // read the large scale sources, which are the variable names with a "ls" suffix
        for (auto & it : lslist)
            profs.get_vector(lsprofs[it],it+"ls", gd.kmax, 0, gd.kstart);

        // Process the time dependent data
        if (tdep_ls.sw)
            create_timedep(tdep_geo,"ls");
    }

    if (swnudge == Nudging_type::enabled)
    {
        // Get profile with nudging factor as function of height
        profs.get_vector(nudge_factor,"nudgefac", gd.kmax, 0, gd.kstart);

        // check whether the fields in the list exist in the prognostic fields
        for (auto & it : nudgelist)
            if (!fields.ap.count(it))
            {
                throw std::runtime_error("field %s in [force][nudgelist] is illegal\n");
            }

        // read the large scale sources, which are the variable names with a "nudge" suffix
        for (auto & it : nudgelist)
            profs.get_vector(nudgeprofs[it],it+"nudge", gd.kmax, 0, gd.kstart);

        // Process the time dependent data
        if (tdep_nudge.sw)
            create_timedep(tdep_nudge,"nudge");
    }

    // Get the large scale vertical velocity from the input
    if (swwls == Large_scale_subsidence_type::enabled)
    {
        profs.get_vector(wls,"wls", gd.kmax, 0, gd.kstart);

        if (tdep_wls.sw)
            create_timedep(tdep_wls,"wls");
    }
}

template <typename TF>
void Force<TF>::create_timedep(Force<TF>::tdep timedep, std::string suffix)
{
    auto& gd = grid.get_grid_data();
    std::vector<TF> tmp;
    for (auto& var : timedep.vars)
    {
        Data_block data_block(master, var+suffix);
        std::vector<std::string>  headers = data_block.get_headers();
        // Sort the times
        std::sort (headers.begin()+1, headers.end());


        for (auto& it : headers)
        {
            timedep.time[var].push_back(std::stod(it));
            data_block.get_vector(tmp, it, gd.kmax, 0, gd.kstart);
            timedep.data[var].insert(timedep.data[var].end(),tmp.begin(),tmp.end());
        }
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

    else if (swlspres== Large_scale_pressure_type::geo_wind)
    {
        if (grid.swspatialorder == "2")
            calc_coriolis_2nd<TF>(fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), ug.data(), vg.data(), fc,
            grid.utrans, grid.vtrans, gd.istart, gd.iend, gd.icells, gd.jstart, gd.jend,
            gd.ijcells, gd.kstart, gd.kend);
        else if (grid.swspatialorder == "4")
            calc_coriolis_4th<TF>(fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), ug.data(), vg.data(), fc,
            grid.utrans, grid.vtrans, gd.istart, gd.iend, gd.icells, gd.jstart, gd.jend,
            gd.ijcells, gd.kstart, gd.kend);
    }

    if (swls == Large_scale_tendency_type::enabled)
    {
        for (auto& it : lslist)
            calc_large_scale_source<TF>(fields.st.at(it)->fld.data(), lsprofs.at(it).data(),
            gd.istart, gd.iend, gd.icells, gd.jstart, gd.jend,
            gd.ijcells, gd.kstart, gd.kend);
    }

    if (swwls == Large_scale_subsidence_type::enabled)
    {
        for (auto& it : fields.st)
            advec_wls_2nd<TF>(fields.st.at(it.first)->fld.data(), fields.sp.at(it.first)->fld_mean.data(), wls.data(), gd.dzhi.data(),
            gd.istart, gd.iend, gd.icells, gd.jstart, gd.jend,
            gd.ijcells, gd.kstart, gd.kend);
    }

    if (swnudge == Nudging_type::enabled)
    {
        for (auto& it : nudgelist)
            calc_nudging_tendency<TF>(fields.st.at(it)->fld.data(), fields.sp.at(it)->fld_mean.data(),
            nudgeprofs.at(it).data(), nudge_factor.data(),
            gd.istart, gd.iend, gd.icells, gd.jstart, gd.jend,
            gd.ijcells, gd.kstart, gd.kend);
    }
}
#endif

#ifndef USECUDA
template <typename TF>
void Force<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (tdep_ls.sw)
        update_time_dependent_profs(timeloop, lsprofs, tdep_ls);
    if (tdep_nudge.sw)
        update_time_dependent_profs(timeloop, nudgeprofs, tdep_nudge);
    if (tdep_geo.sw)
    {
        update_time_dependent_prof(timeloop, ug, tdep_geo,"u");
        update_time_dependent_prof(timeloop, vg, tdep_geo,"v");
    }
    if (tdep_wls.sw)
        update_time_dependent_prof(timeloop, wls, tdep_wls,"w");
}
#endif

template <typename TF>
void Force<TF>::update_time_dependent_profs(Timeloop<TF>& timeloop, std::map<std::string, std::vector<TF>> profiles, tdep timedep)
{
    auto& gd = grid.get_grid_data();
    const int kk = gd.kmax;
    const int kgc = gd.kgc;

    // Loop over all profiles which might be time dependent
    for (auto& it : timedep.data)
    {
        // Get/calculate the interpolation indexes/factors. Assing to zero to avoid compiler warnings.
        int index0 = 0, index1 = 0;
        TF fac0 = 0., fac1 = 0.;

        timeloop.get_interpolation_factors(index0, index1, fac0, fac1, timedep.time[it.first]);

        // Calculate the new vertical profile
        for (int k=0; k<gd.kmax; ++k    )
            profiles[it.first][k+kgc] = fac0 * it.second[index0*kk+k] + fac1 * it.second[index1*kk+k];
    }
}

template <typename TF>
void Force<TF>::update_time_dependent_prof(Timeloop<TF>& timeloop, std::vector<TF> prof, tdep timedep, const std::string& name)
{
    auto& gd = grid.get_grid_data();
    const int kk = gd.kmax;
    const int kgc = gd.kgc;

    // Get/calculate the interpolation indexes/factors
    int index0 = 0, index1 = 0;
    TF fac0 = 0., fac1 = 0.;
    timeloop.get_interpolation_factors(index0, index1, fac0, fac1, timedep.time[name]);

    // Calculate the new vertical profile
    for (int k=0; k<gd.kmax; ++k)
        prof[k+kgc] = fac0 * timedep.data[name][index0*kk+k] + fac1 * timedep.data[name][index1*kk+k];
}

template class Force<double>;
template class Force<float>;
