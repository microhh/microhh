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
}

template<typename TF>
Force<TF>::Force(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), field3d_operators(masterin, gridin, fieldsin)
{
    // Read the switches from the input
    std::string swlspres_in = inputin.get_item<std::string>("force", "swlspres", "", "0");

    // Set the internal switches and read other required input
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
    }
    else
        throw std::runtime_error("Invalid option for \"swlspres\"");
}

template <typename TF>
Force<TF>::~Force()
{
}

template <typename TF>
void Force<TF>::init()
{
//    if (swlspres == "geo")
//    {
//        ug = new double[grid->kcells];
//        vg = new double[grid->kcells];
//    }
//
//    if (swls == "1")
//    {
//        for (std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
//            lsprofs[*it] = new double[grid->kcells];
//    }
//
//    if (swwls == "1")
//        wls = new double[grid->kcells];
}

template <typename TF>
void Force<TF>::create(Input& inputin)
{
//    int nerror = 0;
//
//    if (swlspres == "geo")
//    {
//        nerror += inputin->get_prof(&ug[grid->kstart], "ug", grid->kmax);
//        nerror += inputin->get_prof(&vg[grid->kstart], "vg", grid->kmax);
//    }
//
//    if (swls == "1")
//    {
//        // check whether the fields in the list exist in the prognostic fields
//        for (std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
//            if (!fields->ap.count(*it))
//            {
//                master->print_error("field %s in [force][lslist] is illegal\n", it->c_str());
//                ++nerror;
//            }
//
//        // read the large scale sources, which are the variable names with a "ls" suffix
//        for (std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
//            nerror += inputin->get_prof(&lsprofs[*it][grid->kstart], *it+"ls", grid->kmax);
//    }
//
//    if (swwls == "1")
//        nerror += inputin->get_prof(&wls[grid->kstart], "wls", grid->kmax);
//
//    // process the profiles for the time dependent data
//    if (swtimedep == "1")
//    {
//        // create temporary list to check which entries are used
//        std::vector<std::string> tmplist = timedeplist;
//
//        // process time dependent bcs for the large scale forcings
//        for (std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
//        {
//            // \TODO make sure to give each element its own time series and remove the clear()
//            timedeptime.clear();
//            std::string name = *it + "ls";
//            if (std::find(timedeplist.begin(), timedeplist.end(), *it) != timedeplist.end())
//            {
//                nerror += inputin->get_time_prof(&timedepdata[name], &timedeptime, name, grid->kmax);
//
//                // remove the item from the tmplist
//                std::vector<std::string>::iterator ittmp = std::find(tmplist.begin(), tmplist.end(), *it);
//                if (ittmp != tmplist.end())
//                    tmplist.erase(ittmp);
//            }
//        }
//
//        // display a warning for the non-supported
//        for (std::vector<std::string>::const_iterator ittmp=tmplist.begin(); ittmp!=tmplist.end(); ++ittmp)
//            master->print_warning("%s is not supported (yet) as a time dependent parameter\n", ittmp->c_str());
//    }
//
//    if (nerror)
//        throw 1;
}

//#ifndef USECUDA
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

    //else if (swlspres == "geo")
    //{
    //    if (grid->swspatialorder == "2")
    //        calc_coriolis_2nd(fields->ut->data, fields->vt->data, fields->u->data, fields->v->data, ug, vg);
    //    else if (grid->swspatialorder == "4")
    //        calc_coriolis_4th(fields->ut->data, fields->vt->data, fields->u->data, fields->v->data, ug, vg);
    //}

    //if (swls == "1")
    //{
    //    for (std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
    //        calc_large_scale_source(fields->st[*it]->data, lsprofs[*it]);
    //}

    //if (swwls == "1")
    //{
    //    for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); ++it)
    //        advec_wls_2nd(it->second->data, fields->sp[it->first]->fld_mean, wls, grid->dzhi);
    //}
}
//#endif

//Force::Force(Model* modelin, Input* inputin)
//{


    //model  = modelin;
    //grid   = model->grid;
    //fields = model->fields;
    //master = model->master;

    //ug  = 0;
    //vg  = 0;
    //wls = 0;

    //ug_g  = 0;
    //vg_g  = 0;
    //wls_g = 0;

    //int nerror = 0;
    //nerror += inputin->get_item(&swlspres, "force", "swlspres", "", "0");
    //nerror += inputin->get_item(&swls    , "force", "swls"    , "", "0");
    //nerror += inputin->get_item(&swwls   , "force", "swwls"   , "", "0");

    //if (swlspres != "0")
    //{
    //    if (swlspres == "uflux")
    //        nerror += inputin->get_item(&uflux, "force", "uflux", "");
    //    else if (swlspres == "geo")
    //        nerror += inputin->get_item(&fc, "force", "fc", "");
    //    else
    //    {
    //        ++nerror;
    //        master->print_error("\"%s\" is an illegal option for swlspres\n", swlspres.c_str());
    //    }
    //}

    //if (swls == "1")
    //    nerror += inputin->get_list(&lslist, "force", "lslist", "");
    //else if (swls != "0")
    //{
    //    ++nerror;
    //    master->print_error("\"%s\" is an illegal option for swls\n", swls.c_str());
    //}

    //if (swwls == "1")
    //    fields->set_calc_mean_profs(true);
    //else if (swwls != "0")
    //{
    //    ++nerror;
    //    master->print_error("\"%s\" is an illegal option for swwls\n", swwls.c_str());
    //}

    //// get the list of time varying variables
    //nerror += inputin->get_item(&swtimedep  , "force", "swtimedep"  , "", "0");
    //nerror += inputin->get_list(&timedeplist, "force", "timedeplist", "");

    //if (nerror)
    //    throw 1;
//}

//Force::~Force()
//{
//    delete[] ug;
//    delete[] vg;
//    delete[] wls;
//
//    if (swls == "1")
//    {
//        for (std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
//            delete[] lsprofs[*it];
//    }
//
//    // clean up time dependent data
//    for (std::map<std::string, double*>::const_iterator it=timedepdata.begin(); it!=timedepdata.end(); ++it)
//        delete[] it->second;
//
//#ifdef USECUDA
//    clear_device();
//#endif
//}
//
//
//
//void Force::update_time_dependent()
//{
//    if (swtimedep == "0")
//        return;
//
//    // first find the index for the time entries
//    unsigned int index0 = 0;
//    unsigned int index1 = 0;
//    for (std::vector<double>::const_iterator it=timedeptime.begin(); it!=timedeptime.end(); ++it)
//    {
//        if (model->timeloop->get_time() < *it)
//            break;
//        else
//            ++index1;
//    }
//
//    // second, calculate the weighting factor
//    double fac0, fac1;
//
//    // correct for out of range situations where the simulation is longer than the time range in input
//    if (index1 == 0)
//    {
//        fac0 = 0.;
//        fac1 = 1.;
//        index0 = 0;
//    }
//    else if (index1 == timedeptime.size())
//    {
//        fac0 = 1.;
//        fac1 = 0.;
//        index0 = index1-1;
//        index1 = index0;
//    }
//    else
//    {
//        index0 = index1-1;
//        double timestep;
//        timestep = timedeptime[index1] - timedeptime[index0];
//        fac0 = (timedeptime[index1] - model->timeloop->get_time()) / timestep;
//        fac1 = (model->timeloop->get_time() - timedeptime[index0]) / timestep;
//    }
//
//    update_time_dependent_profs(fac0, fac1, index0, index1);
//}
//
//#ifndef USECUDA
//void Force::update_time_dependent_profs(const double fac0, const double fac1, const int index0, const int index1)
//{
//    // process time dependent bcs for the large scale forcings
//    const int kk = grid->kmax;
//    const int kgc = grid->kgc;
//
//    for (std::vector<std::string>::const_iterator it1=lslist.begin(); it1!=lslist.end(); ++it1)
//    {
//        std::string name = *it1 + "ls";
//        std::map<std::string, double*>::const_iterator it2 = timedepdata.find(name);
//
//        // update the profile
//        if (it2 != timedepdata.end())
//            for (int k=0; k<grid->kmax; ++k)
//                lsprofs[*it1][k+kgc] = fac0*it2->second[index0*kk+k] + fac1*it2->second[index1*kk+k];
//    }
//}
//#endif
//
//void Force::calc_flux(double* const restrict ut, const double* const restrict u,
//                      const double* const restrict dz, const double dt)
//{
//    const int jj = grid->icells;
//    const int kk = grid->ijcells;
//
//    const double ugrid = grid->utrans;
//
//    double uavg  = 0.;
//    double utavg = 0.;
//
//    for (int k=grid->kstart; k<grid->kend; ++k)
//        for (int j=grid->jstart; j<grid->jend; ++j)
//#pragma ivdep
//            for (int i=grid->istart; i<grid->iend; ++i)
//            {
//                const int ijk = i + j*jj + k*kk;
//                uavg  += u [ijk]*dz[k];
//                utavg += ut[ijk]*dz[k];
//            }
//
//    grid->get_sum(&uavg);
//    grid->get_sum(&utavg);
//
//    uavg  = uavg  / (grid->itot*grid->jtot*grid->zsize);
//    utavg = utavg / (grid->itot*grid->jtot*grid->zsize);
//
//    double fbody;
//    fbody = (uflux - uavg - ugrid) / dt - utavg;
//
//    for (int n=0; n<grid->ncells; n++)
//        ut[n] += fbody;
//}
//
//void Force::calc_coriolis_2nd(double* const restrict ut, double* const restrict vt,
//                              const double* const restrict u , const double* const restrict v ,
//                              const double* const restrict ug, const double* const restrict vg)
//{
//    const int ii = 1;
//    const int jj = grid->icells;
//    const int kk = grid->ijcells;
//
//    const double ugrid = grid->utrans;
//    const double vgrid = grid->vtrans;
//
//    for (int k=grid->kstart; k<grid->kend; ++k)
//        for (int j=grid->jstart; j<grid->jend; ++j)
//#pragma ivdep
//            for (int i=grid->istart; i<grid->iend; ++i)
//            {
//                const int ijk = i + j*jj + k*kk;
//                ut[ijk] += fc * (0.25*(v[ijk-ii] + v[ijk] + v[ijk-ii+jj] + v[ijk+jj]) + vgrid - vg[k]);
//            }
//
//    for (int k=grid->kstart; k<grid->kend; ++k)
//        for (int j=grid->jstart; j<grid->jend; ++j)
//#pragma ivdep
//            for (int i=grid->istart; i<grid->iend; ++i)
//            {
//                const int ijk = i + j*jj + k*kk;
//                vt[ijk] -= fc * (0.25*(u[ijk-jj] + u[ijk] + u[ijk+ii-jj] + u[ijk+ii]) + ugrid - ug[k]);
//            }
//}
//
//void Force::calc_coriolis_4th(double* const restrict ut, double* const restrict vt,
//                              const double* const restrict u , const double* const restrict v ,
//                              const double* const restrict ug, const double* const restrict vg)
//{
//    using namespace Finite_difference::O4;
//
//    const int ii1 = 1;
//    const int ii2 = 2;
//    const int jj1 = 1*grid->icells;
//    const int jj2 = 2*grid->icells;
//    const int kk1 = 1*grid->ijcells;
//
//    const double ugrid = grid->utrans;
//    const double vgrid = grid->vtrans;
//
//    for (int k=grid->kstart; k<grid->kend; ++k)
//        for (int j=grid->jstart; j<grid->jend; ++j)
//#pragma ivdep
//            for (int i=grid->istart; i<grid->iend; ++i)
//            {
//                const int ijk = i + j*jj1 + k*kk1;
//                ut[ijk] += fc * ( ( ci0*(ci0*v[ijk-ii2-jj1] + ci1*v[ijk-ii1-jj1] + ci2*v[ijk-jj1] + ci3*v[ijk+ii1-jj1])
//                                  + ci1*(ci0*v[ijk-ii2    ] + ci1*v[ijk-ii1    ] + ci2*v[ijk    ] + ci3*v[ijk+ii1    ])
//                                  + ci2*(ci0*v[ijk-ii2+jj1] + ci1*v[ijk-ii1+jj1] + ci2*v[ijk+jj1] + ci3*v[ijk+ii1+jj1])
//                                  + ci3*(ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii1+jj2] + ci2*v[ijk+jj2] + ci3*v[ijk+ii1+jj2]) )
//                                + vgrid - vg[k] );
//            }
//
//    for (int k=grid->kstart; k<grid->kend; ++k)
//        for (int j=grid->jstart; j<grid->jend; ++j)
//#pragma ivdep
//            for (int i=grid->istart; i<grid->iend; ++i)
//            {
//                const int ijk = i + j*jj1 + k*kk1;
//                vt[ijk] -= fc * ( ( ci0*(ci0*u[ijk-ii1-jj2] + ci1*u[ijk-jj2] + ci2*u[ijk+ii1-jj2] + ci3*u[ijk+ii2-jj2])
//                                  + ci1*(ci0*u[ijk-ii1-jj1] + ci1*u[ijk-jj1] + ci2*u[ijk+ii1-jj1] + ci3*u[ijk+ii2-jj1])
//                                  + ci2*(ci0*u[ijk-ii1    ] + ci1*u[ijk    ] + ci2*u[ijk+ii1    ] + ci3*u[ijk+ii2    ])
//                                  + ci3*(ci0*u[ijk-ii1+jj1] + ci1*u[ijk+jj1] + ci2*u[ijk+ii1+jj1] + ci3*u[ijk+ii2+jj1]) )
//                                + ugrid - ug[k]);
//            }
//}
//
//void Force::calc_large_scale_source(double* const restrict st, const double* const restrict sls)
//{
//    const int jj = grid->icells;
//    const int kk = grid->ijcells;
//
//    for (int k=grid->kstart; k<grid->kend; ++k)
//        for (int j=grid->jstart; j<grid->jend; ++j)
//            for (int i=grid->istart; i<grid->iend; ++i)
//            {
//                const int ijk = i + j*jj + k*kk;
//                st[ijk] += sls[k];
//            }
//}
//
//void Force::advec_wls_2nd(double* const restrict st, const double* const restrict s,
//                          const double* const restrict wls, const double* const dzhi)
//{
//    const int jj = grid->icells;
//    const int kk = grid->ijcells;
//
//    // use an upwind differentiation
//    for (int k=grid->kstart; k<grid->kend; ++k)
//    {
//        if (wls[k] > 0.)
//        {
//            for (int j=grid->jstart; j<grid->jend; ++j)
//                for (int i=grid->istart; i<grid->iend; ++i)
//                {
//                    const int ijk = i + j*jj + k*kk;
//                    st[ijk] -=  wls[k] * (s[k]-s[k-1])*dzhi[k];
//                }
//        }
//        else
//        {
//            for (int j=grid->jstart; j<grid->jend; ++j)
//                for (int i=grid->istart; i<grid->iend; ++i)
//                {
//                    const int ijk = i + j*jj + k*kk;
//                    st[ijk] -=  wls[k] * (s[k+1]-s[k])*dzhi[k+1];
//                }
//        }
//    }
//}


template class Force<double>;
template class Force<float>;
