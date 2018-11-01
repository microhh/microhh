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

    template<typename TF> //EW: enforce B.C. to be from a mesoscale model, input format is pending
    void enforce_BC_from_mesoscale_model(
    		TF* const restrict data, const int fixed_value,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int jcells, const int kcells, const int igc, const int jgc, const int jtot,
			Edge edge)
    {
        const int jj = icells;
        const int kk = icells*jcells;
        if (edge == Edge::East_west_edge || edge == Edge::Both_edges)
        {
            // first, east west boundaries
            for (int k=0; k<kcells; ++k)
                for (int j=0; j<jcells; ++j)
                    #pragma ivdep
                    for (int i=0; i<igc+1; ++i)
                    {
                        const int ijk0 = i          + j*jj + k*kk;
                        data[ijk0] = fixed_value;
                    }

            for (int k=0; k<kcells; ++k)
                for (int j=0; j<jcells; ++j)
                    #pragma ivdep
                    for (int i=0; i<igc+1; ++i)
                    {
                        const int ijk0 = i+iend   + j*jj + k*kk;
                        data[ijk0] = fixed_value;
                    }
        }
        if (edge == Edge::North_south_edge || edge == Edge::Both_edges)
        {
            // if the run is 3D, apply the BCs
            if (jtot > 1)
            {
                // second, send and receive the ghost cells in the north-south direction
                for (int k=0; k<kcells; ++k)
                    for (int j=0; j<jgc+1; ++j)
                        #pragma ivdep
                        for (int i=0; i<icells; ++i)
                        {
                            const int ijk0 = i + j                 *jj + k*kk;
                            data[ijk0] = fixed_value;
                        }

                for (int k=0; k<kcells; ++k)
                    for (int j=0; j<jgc+1; ++j)
                        #pragma ivdep
                        for (int i=0; i<icells; ++i)
                        {
                            const int ijk0 = i + (j+jend  )*jj + k*kk;
                            data[ijk0] = fixed_value;
                        }
            }
            // in case of 2D, fill all the ghost cells with the current value
            else
            {
                for (int k=kstart; k<kend; ++k)
                    for (int j=0; j<jgc+1; ++j)
                        #pragma ivdep
                        for (int i=0; i<icells; ++i)
                        {
                            const int ijknorth = i + j*jj           + k*kk;
                            const int ijksouth = i + (jend+j)*jj + k*kk;
                            data[ijknorth] = fixed_value;
                            data[ijksouth] = fixed_value;
                        }
            }
        }
    }
    //ToDo: initialize 1d, 2d, 3d arrays properly
    //sunray for tau and swn
    //zenith for mu
    template<typename TF> //EW: simplified radiative parameterization for LW and SW fluxes for DYCOMS
    void gcss_rad(TF* const restrict st, const TF* const restrict s,
            const TF* const restrict wls, const TF* const dzi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells){
    	
        const int jj = icells;
        const int kk = ijcells;
        const TF xka = 85.;
        const TF fr0 = 70.;
        const TF fr1 = 22.;
        const TF rho_l = 1000.;
        const TF reff = 1.E-5;
        const TF cp = 1005; //can read this from constant.h
        const TF div = 3.75E-6; //fix divergence for now
        TF tauc; //single or double precision??
        TF fact;
        int ki; //PBLH index
        auto lwp = fields.get_tmp(); //how do I set lwp = 0?
        auto flx = fields.get_tmp();
        auto tau = fields.get_tmp();
        const TF mu = 0.05;//zenith(32.5,time); //zenith
        for (int j=jstart; j<jend; ++j){
            for (int i=istart; i<iend; ++i){
                ki = kend; //set to top of domain
                for (int k=kstart; k<kend; ++k)
                {
                    const int ij   = i + j*jj;
                    const int ijk  = i + j*jj + k*kk;
                    const int km1 = std::max(1,k-1);
                    lwp[ij] = lwp[ij]+std::max(0.,ql[ijk]*fields.rhoref[k]*(dzi[k]-dzi[km1]));
                    flx[ijk] = fr1*std::exp(-1.*xka*lwp[ij]);
                    if ((ql[ijk]>0.01E-3)&&(qt[ijk]>=0.008)) ki = k; //this is the PBLH index
                }
				
                if (mu>0.035){
                    tauc = 0.0;
                    for (k=kstart;k<kend;++k){
                        const int ij   = i + j*jj;
                        const int ijk  = i + j*jj + k*kk;
                        const int km1 = std::max(1,k-1);
                        tau[k] = 0.0;
                        if (ql[ijk]>1.E-5){
                            tau[k] = std::max(0.,1.5*ql[ijk]*fields.rhoref[k]*(dzi[k]-dzi[km1])/reff/rho_l);
                            tauc = tauc + tau[k];
                        }
                    }
                    //sunray
                    swn = 1.0; //SW
                }
                fact = div*cp*fields.rhoref[ki];
                flx[i+j*jj+kstart*kk] = flx[i+j*jj+kstart*kk]+fr0*std::exp(-1.*xka*lwp[ij]);
                for (int k=kstart+1;k<kend;++k){
                    const int ij   = i + j*jj;
                    const int ijk  = i + j*jj + k*kk;
                    const int km1 = std::max(kstart+1,k-1);
                    const int ijkm = i + j*jj + km1*kk;
                    lwp[ij] = lwp[ij]-std::max(0.,ql[ijk]*fields.rhoref[k]*(dzi[k]-dzi[km1]));
                    flx[ijk] = flx[ijk]+fr0*std::exp(-1.*xka*lwp[ij]);
                    if ((k>ki)&&(ki>1)&&(fact>0.)){ //above PBLH
                        flx[ijk] = flx[ijk] + fact*(0.25*std::pow(z[k]-z[km],1.333)+z[km]*std::pow(z[k]-z[ki],0.33333))
                    }
                    tt[ijk] = tt[ijk]-(flx[ijk]-flx[ijkm])*dzh[k]/(fields.rhoref[k]*cp);
                    tt[ijk] = tt[ijk]+(swn[ijk]-swn[ijkm])*dzh[k]/(fields.rhoref[k]*cp);
                }
                //subsidence part
                if (div!=0.){
                    for (int k=kstart+1;k<kend-2,++k){
                        const int ijk  = i + j*jj + k*kk;
                        const int kp1 = k+1;
                        const int ijkp1 = i + j*jj + kp1*kk;
                        tt[ijk] = tt[ijk] + div*z[k]*(tl[ijkp1]-tl[ijk])*dzi[k];
                        qtt[ijk] = qtt[ijk] + div*z[k]*(qt[ijkp1]-qt[ijk])*dzi[k];
                    }
                }
				
            } // end of i
        } // end of j
    } // end of gcss_rad
}

template<typename TF>
Force<TF>::Force(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), field3d_operators(masterin, gridin, fieldsin)
{
    std::string swlspres_in = inputin.get_item<std::string>("force", "swlspres", "", "0");
    std::string swls_in     = inputin.get_item<std::string>("force", "swls", "", "0");
    std::string swwls_in    = inputin.get_item<std::string>("force", "swwls", "", "0");
    std::string swnudge_in  = inputin.get_item<std::string>("force", "swnudge", "", "0");
    std::string mesoscalebc_in = inputin.get_item<std::string>("force", "mesoscalebc", "", "0");

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
        nudgelist = inputin.get_list<std::string>("force", "nudgelist", "", std::vector<std::string>());

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
    // Enforce BC from some mesoscale model
    if (mesoscalebc_in == "0")
    	mesoscalebc = Mesoscale_BC_type::disabled;
    else if (mesoscalebc_in == "1")
    {
    	mesoscalebc = Mesoscale_BC_type::enabled;
    }
    else
    {
    	throw std::runtime_error("Invalid option for \"mesoscalebc\"");
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

        for (auto& it : tdep_geo)
            it.second->create_timedep_prof();
    }

    if (swls == Large_scale_tendency_type::enabled)
    {
        // check whether the fields in the list exist in the prognostic fields
        for (auto & it : lslist)
            if (!fields.ap.count(it))
            {
                throw std::runtime_error("field %s in [force][lslist] is illegal\n");
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
                throw std::runtime_error("field %s in [force][nudgelist] is illegal\n");
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
            calc_nudging_tendency<TF>(
                    fields.st.at(it)->fld.data(), fields.sp.at(it)->fld_mean.data(),
                    nudgeprofs.at(it).data(), nudge_factor.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
    }

    if (mesoscalebc == Mesoscale_BC_type::enabled)
    {   //hard-coding a fixed u value at the BC for testing purposes
    	enforce_BC_from_mesoscale_model(fields.sp.at("th")->fld.data(), 285,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
				gd.icells, gd.jcells, gd.kcells, gd.igc, gd.jgc, gd.jtot,
				Edge::Both_edges);
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