/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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
#include <cmath>
// #include <math.h>
#include <iostream>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "field3d_operators.h"
#include "timedep.h"
#include "force.h"
#include "constants.h"
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
    void interpolate_geo_3d(
            TF* const restrict ug,
            TF* const restrict vg,
            const TF* const restrict ug_prev,
            const TF* const restrict vg_prev,
            const TF* const restrict ug_next,
            const TF* const restrict vg_next,
            const TF tfac,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        const TF tfac1 = TF(1) - tfac;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jstride + k*kstride;

                    ug[ijk] = tfac * ug_prev[ijk] + tfac1 * ug_next[ijk];
                    vg[ijk] = tfac * vg_prev[ijk] + tfac1 * vg_next[ijk];
                }
    }


    template<typename TF>
    void calc_coriolis_ls(
            TF* const restrict ut, TF* const restrict vt,
            const TF* const restrict u, const TF* const restrict v,
            const TF* const restrict ug, const TF* const restrict vg, const TF* const restrict fc_2d,
            const TF* const rhoref,
            const TF ugrid, const TF vgrid,
            const TF dxi, const TF dyi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        for (int k=kstart; k<kend; ++k)
        {
            const TF rhorefi = TF(1.) / rhoref[k];

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    const int ijk = i + j*jj + k*kk;
                    const TF fc = TF(0.5)*(fc_2d[ij-ii] + fc_2d[ij]);
                    ut[ijk] += fc * (TF(0.25)*(v[ijk-ii] + v[ijk] + v[ijk-ii+jj] + v[ijk+jj]) + vgrid - vg[ijk]);
                }
        }

        for (int k=kstart; k<kend; ++k)
        {
            const TF rhorefi = TF(1.) / rhoref[k];

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    const int ijk = i + j*jj + k*kk;
                    const TF fc = TF(0.5)*(fc_2d[ij-jj] + fc_2d[ij]);
                    vt[ijk] += - fc * (TF(0.25)*(u[ijk-jj] + u[ijk] + u[ijk+ii-jj] + u[ijk+ii]) + ugrid - ug[ijk]);
                }
        }
    }


    template<typename TF>
    void calc_rotation(
            TF* const restrict ut,
            TF* const restrict vt,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict fc_2d,
            const TF ugrid, const TF vgrid,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
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
                    const int ij = i + j*jj;
                    const int ijk = i + j*jj + k*kk;

                    const TF fc_u = TF(0.5)*(fc_2d[ij-ii] + fc_2d[ij]);
                    const TF fc_v = TF(0.5)*(fc_2d[ij-jj] + fc_2d[ij]);

                    ut[ijk] +=  fc_u * (TF(0.25)*(v[ijk-ii] + v[ijk] + v[ijk-ii+jj] + v[ijk+jj]) + vgrid);
                    vt[ijk] += -fc_v * (TF(0.25)*(u[ijk-jj] + u[ijk] + u[ijk+ii-jj] + u[ijk+ii]) + ugrid);
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
            TF* const restrict st,
            const TF* const restrict s,
            const TF* const restrict wls,
            const TF* const dzi,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        for (int k=kstart+1; k<kend; ++k)
        {
            if ( interp2( wls[k-1], wls[k] ) > 0.)
            {
                for (int j=jstart; j<jend; ++j)
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        st[ijk] -=  interp2( wls[k-1], wls[k] ) * (s[ijk]-s[ijk-kk])*dzi[k-1];
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
    master(masterin), grid(gridin), fields(fieldsin), field3d_operators(masterin, gridin, fieldsin), boundary_cyclic(masterin, gridin)
{
    std::string swlspres_in = inputin.get_item<std::string>("force", "swlspres", "", "0");
    std::string swwls_in    = inputin.get_item<std::string>("force", "swwls"   , "", "0");
    std::string swnudge_in  = inputin.get_item<std::string>("force", "swnudge" , "", "0");
    std::string swls_in     = inputin.get_item<std::string>("force", "swls"    , "", "0");

    swrotation_2d = inputin.get_item<bool>("force", "swrotation_2d", "", false);

    if (swwls_in == "mean" || swwls_in == "local")
        swwls_mom = inputin.get_item<bool>("force", "swwls_mom", "", false);

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
        fc = inputin.get_item<TF>("force", "fc", "", -1.);
        // Test whether latitude is available in the input file.
        if (fc < 0.)
            inputin.get_item<TF>("grid", "lat", "");

        swtimedep_geo = inputin.get_item<bool>("force", "swtimedep_geo", "", false);

        tdep_geo.emplace("u_geo", new Timedep<TF>(master, grid, "u_geo", swtimedep_geo));
        tdep_geo.emplace("v_geo", new Timedep<TF>(master, grid, "v_geo", swtimedep_geo));
    }
    else if (swlspres_in == "geo3d")
    {
        swlspres = Large_scale_pressure_type::Geo_wind_3d;

        fc = inputin.get_item<TF>("force", "fc", ""); // CvH I still load the mean fc, maybe this is asking for trouble?
        swtimedep_geo = inputin.get_item<TF>("force", "swtimedep_geo", "", false);

        if (swtimedep_geo)
            ugeo_loadtime = inputin.get_item<int>("force", "ugeo_loadtime", "");
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
        swtimedep_ls = inputin.get_item<bool>("force", "swtimedep_ls", "", false);

        if (swtimedep_ls)
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

    swtimedep_wls = inputin.get_item<bool>("force", "swtimedep_wls", "", false);
    tdep_wls = std::make_unique<Timedep<TF>>(master, grid, "w_ls", swtimedep_wls);

    // Nudging.
    if (swnudge_in == "0")
        swnudge = Nudging_type::Disabled;
    else if (swnudge_in == "1")
    {
        swnudge = Nudging_type::Enabled;
        nudgelist       = inputin.get_list<std::string>("force", "nudgelist", "", std::vector<std::string>());
        scalednudgelist = inputin.get_list<std::string>("force", "scalednudgelist", "", std::vector<std::string>());
        swtimedep_nudge = inputin.get_item<bool>("force", "swtimedep_nudge", "", false);

        if (swtimedep_nudge)
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
    else if (swlspres == Large_scale_pressure_type::Geo_wind_3d)
    {
        ug.resize(gd.ncells);
        vg.resize(gd.ncells);
        fc_2d.resize(gd.ijcells);
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

    if (swrotation_2d)
        fc_2d.resize(gd.ijcells);

    boundary_cyclic.init();
}

template <typename TF>
void Force<TF>::create(Input& inputin, Netcdf_handle& input_nc, Stats<TF>& stats, Timeloop<TF>& timeloop)
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
        const TF offset = 0;
        if (!swtimedep_geo)
        {
            group_nc.get_variable(ug, "u_geo", {0}, {gd.ktot});
            std::rotate(ug.rbegin(), ug.rbegin() + gd.kstart, ug.rend());
        }
        else
            tdep_geo["u_geo"]->create_timedep_prof(input_nc, offset, timedep_dim);

        if (!swtimedep_geo)
        {
            group_nc.get_variable(vg, "v_geo", {0}, {gd.ktot});
            std::rotate(vg.rbegin(), vg.rbegin() + gd.kstart, vg.rend());
        }
        else
            tdep_geo["v_geo"]->create_timedep_prof(input_nc, offset, timedep_dim);

        stats.add_tendency(*fields.mt.at("u"), "z", tend_name_cor, tend_longname_cor);
        stats.add_tendency(*fields.mt.at("v"), "z", tend_name_cor, tend_longname_cor);

    }
    else if (swlspres == Large_scale_pressure_type::Geo_wind_3d)
    {
        constexpr TF no_offset = TF(0);

        constexpr int n = 0;
        auto tmp1 = fields.get_tmp();
        auto tmp2 = fields.get_tmp();

        // CvH: for now assume full 3D size.
        Field3d_io field3d_io(master, grid);
        int nerror = 0;

        // Load the coriolis parameter.
        char filename[256];
        std::sprintf(filename, "fc.%07d", n);
        master.print_message("Loading \"%s\" ... ", filename);
        if (field3d_io.load_xy_slice(fc_2d.data(), tmp1->fld.data(), filename))
            master.print_message("FAILED\n");
        else
            master.print_message("OK\n");

        // I do here a cyclic BC, but that is (slightly) incorrect at the edges, but those tendencies we do not need with open BCs.
        boundary_cyclic.exec_2d(fc_2d.data());

        // Read the geowind 3d files.
        auto read_3d_binary = [&](const std::string& name, TF* fld, int iotime)
        {
            char filename[256];
            std::sprintf(filename, "%s.%07d", name.c_str(), iotime);
            master.print_message("Loading \"%s\" ... ", filename);

            // Load the ug geowind.
            if (field3d_io.load_field3d(
                        fld,
                        tmp1->fld.data(), tmp2->fld.data(),
                        filename, no_offset,
                        gd.kstart, gd.kend))
            {
                master.print_message("FAILED\n");
                ++nerror;
            }
            else
                master.print_message("OK\n");
        };

        if (swtimedep_geo)
        {
            // Find previous and next times.
            unsigned long itime = timeloop.get_itime();
            unsigned long iiotimeprec = timeloop.get_iiotimeprec();

            // Read first two input times
            unsigned long ugeo_iloadtime = convert_to_itime(ugeo_loadtime);
            itime_ugeo_prev = itime / ugeo_iloadtime * ugeo_iloadtime;
            itime_ugeo_next = itime_ugeo_prev + ugeo_iloadtime;

            // IO time accounting for iotimeprec
            const unsigned long iotime_prev = int(itime_ugeo_prev / iiotimeprec);
            const unsigned long iotime_next = int(itime_ugeo_next / iiotimeprec);

            ug_prev.resize(gd.ncells);
            vg_prev.resize(gd.ncells);

            ug_next.resize(gd.ncells);
            vg_next.resize(gd.ncells);

            // Read the first two w_top fields.
            read_3d_binary("ug", ug_prev.data(), iotime_prev);
            read_3d_binary("vg", vg_prev.data(), iotime_prev);

            read_3d_binary("ug", ug_next.data(), iotime_next);
            read_3d_binary("vg", vg_next.data(), iotime_next);
        }
        else
        {
            read_3d_binary("ug", ug.data(), 0);
            read_3d_binary("vg", vg.data(), 0);
        }

        fields.release_tmp(tmp1);
        fields.release_tmp(tmp2);
    }

    if (swrotation_2d)
    {
        auto tmp1 = fields.get_tmp();

        Field3d_io field3d_io(master, grid);

        // Load the coriolis parameter.
        char filename[256];
        int n = 0;
        std::sprintf(filename, "fc.%07d", n);
        master.print_message("Loading \"%s\" ... ", filename);
        if (field3d_io.load_xy_slice(fc_2d.data(), tmp1->fld.data(), filename))
            master.print_message("FAILED\n");
        else
            master.print_message("OK\n");

        boundary_cyclic.exec_2d(fc_2d.data());

        fields.release_tmp(tmp1);
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
            if (tdep_ls.find(it) == tdep_ls.end())
            {
                group_nc.get_variable(lsprofs[it], it+"_ls", {0}, {gd.ktot});
                std::rotate(lsprofs[it].rbegin(), lsprofs[it].rbegin() + gd.kstart, lsprofs[it].rend());
            }
            else
            {
                const TF offset = 0;
                tdep_ls[it]->create_timedep_prof(input_nc, offset, timedep_dim);
            }
        }

        for (std::string& it : lslist)
            stats.add_tendency(*fields.at.at(it), "z", tend_name_ls, tend_longname_ls);
    }

    if (swnudge == Nudging_type::Enabled)
    {
        // Get profile with nudging factor as function of height
        group_nc.get_variable(nudge_factor, "nudgefac", {0}, {gd.ktot});
        TF minnudge = 1e-6;
        TF max = *std::max_element(nudge_factor.begin(), nudge_factor.end());
        if (max < minnudge)
        {
            std::string msg = "The maximum value of the nudging factor is smaller than the minimum allowed value of " + std::to_string(minnudge);
            throw std::runtime_error(msg);
        }
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
            if (tdep_nudge.find(it) == tdep_nudge.end())
            {
                group_nc.get_variable(nudgeprofs[it], it+"_nudge", {0}, {gd.ktot});
                std::rotate(nudgeprofs[it].rbegin(), nudgeprofs[it].rbegin() + gd.kstart, nudgeprofs[it].rend());

                // Account for the Galilean transformation
                if (it == "u")
                    add_offset(nudgeprofs[it].data(), -gd.utrans, gd.kstart, gd.kend);
                else if (it == "v")
                    add_offset(nudgeprofs[it].data(), -gd.vtrans, gd.kstart, gd.kend);
            }
            else
            {
                // Account for the Galilean transformation
                TF offset;
                if (it == "u")
                    offset = -gd.utrans;
                else if (it == "v")
                    offset = -gd.vtrans;
                else
                    offset = 0;

                tdep_nudge.at(it)->create_timedep_prof(input_nc, offset, timedep_dim);
                stats.add_tendency(*fields.at.at(it), "z", tend_name_nudge, tend_longname_nudge);
            }
        }
    }

    // Get the large scale vertical velocity from the input
    if (swwls == Large_scale_subsidence_type::Mean_field || swwls == Large_scale_subsidence_type::Local_field)
    {
         const TF offset = 0;
        if (!swtimedep_wls)
        {
            group_nc.get_variable(wls, "w_ls", {0}, {gd.ktot});
            std::rotate(wls.rbegin(), wls.rbegin() + gd.kstart, wls.rend());
        }
        else
        {
            const TF offset = 0;
            tdep_wls->create_timedep_prof(input_nc, offset, timedep_dim);
        }

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
                fields.at.at("u")->fld.data(), uflux, u_mean, ut_mean, gd.utrans, dt,
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
        TF fc_loc = fc;
        if (fc_loc < 0)
            fc_loc = 2. * Constants::e_rot<TF> * std::sin(gd.lat * TF(M_PI) / 180.);

        if (grid.get_spatial_order() == Grid_order::Second)
            calc_coriolis_2nd<TF>(fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), ug.data(), vg.data(), fc_loc,
            gd.utrans, gd.vtrans,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

        else if (grid.get_spatial_order() == Grid_order::Fourth)
            calc_coriolis_4th<TF>(fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), ug.data(), vg.data(), fc_loc,
            gd.utrans, gd.vtrans,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

        stats.calc_tend(*fields.mt.at("u"), tend_name_cor);
        stats.calc_tend(*fields.mt.at("v"), tend_name_cor);
    }
    else if (swlspres == Large_scale_pressure_type::Geo_wind_3d)
    {
        if (grid.get_spatial_order() == Grid_order::Second)
        {
            calc_coriolis_ls<TF>(
                    fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(),
                    fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(),
                    ug.data(), vg.data(), fc_2d.data(),
                    fields.rhoref.data(),
                    gd.utrans, gd.vtrans,
                    gd.dxi, gd.dyi,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
        }
    }


    if (swrotation_2d)
    {
        calc_rotation<TF>(
                fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(),
                fc_2d.data(),
                gd.utrans,
                gd.vtrans,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
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
            // Also apply to the velocity components u,v.
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
    else if (swlspres == Large_scale_pressure_type::Geo_wind_3d && swtimedep_geo)
    {
        auto& gd = grid.get_grid_data();

        constexpr int n = 0;
        constexpr TF no_offset = TF(0);

        unsigned long itime = timeloop.get_itime();

        if (itime > itime_ugeo_next)
        {
            // Read new w_top field
            unsigned long iiotimeprec = timeloop.get_iiotimeprec();
            unsigned long ugeo_iloadtime = convert_to_itime(ugeo_loadtime);

            itime_ugeo_prev = itime_ugeo_next;
            itime_ugeo_next = itime_ugeo_prev + ugeo_iloadtime;

            const int iotime1 = int(itime_ugeo_next / iiotimeprec);

            // Copy of data from next to prev time
            ug_prev = ug_next;
            vg_prev = vg_next;

            Field3d_io field3d_io(master, grid);
            int nerror = 0;

            auto tmp1 = fields.get_tmp();
            auto tmp2 = fields.get_tmp();

            auto read_3d_binary = [&](const std::string& name, TF* fld)
            {
                char filename[256];
                std::sprintf(filename, "%s.%07d", name.c_str(), iotime1);
                master.print_message("Loading \"%s\" ... ", filename);

                // Load the ug geowind.
                if (field3d_io.load_field3d(
                            fld,
                            tmp1->fld.data(), tmp2->fld.data(),
                            filename, no_offset,
                            gd.kstart, gd.kend))
                {
                    master.print_message("FAILED\n");
                    ++nerror;
                }
                else
                    master.print_message("OK\n");
            };

            read_3d_binary("ug", ug_next.data());
            read_3d_binary("vg", vg_next.data());

            master.sum(&nerror, 1);

            if (nerror)
                throw std::runtime_error("Error in reading time dependent geowind.");

            fields.release_tmp(tmp1);
            fields.release_tmp(tmp2);
        }

        // Update 3D geo-wind field.
        const TF tfac = TF(1) - (TF(itime - itime_ugeo_prev) / TF(itime_ugeo_next - itime_ugeo_prev));

        interpolate_geo_3d(
                ug.data(),
                vg.data(),
                ug_prev.data(),
                vg_prev.data(),
                ug_next.data(),
                vg_next.data(),
                tfac,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }

    if (swwls == Large_scale_subsidence_type::Mean_field || swwls == Large_scale_subsidence_type::Local_field )
        tdep_wls->update_time_dependent_prof(wls, timeloop);
}
#endif


#ifdef FLOAT_SINGLE
template class Force<float>;
#else
template class Force<double>;
#endif
