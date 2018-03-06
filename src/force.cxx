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
#include "force.h"
#include "defines.h"
#include "finite_difference.h"
#include "model.h"
#include "timeloop.h"
#include "boundary.h"

using namespace Finite_difference::O2;

Force::Force(Model* modelin, Input* inputin)
{
    model  = modelin;
    grid   = model->grid;
    fields = model->fields;
    master = model->master;

    ug  = 0;
    vg  = 0;
    wls = 0;
    nudge_factor = 0;

    ug_g  = 0;
    vg_g  = 0;
    wls_g = 0;
    nudge_factor_g = 0;

    int nerror = 0;
    nerror += inputin->get_item(&swlspres, "force", "swlspres", "", "0");
    nerror += inputin->get_item(&swls    , "force", "swls"    , "", "0");
    nerror += inputin->get_item(&swwls   , "force", "swwls"   , "", "0");
    nerror += inputin->get_item(&swnudge , "force", "swnudge" , "", "0");

    if (swlspres != "0")
    {
        if (swlspres == "uflux")
            nerror += inputin->get_item(&uflux, "force", "uflux", "");
        else if (swlspres == "geo")
        {
            nerror += inputin->get_item(&fc, "force", "fc", "");
            nerror += inputin->get_item(&swtimedep_geo, "force", "swtimedep_geo", "", "0");
        }
        else
        {
            ++nerror;
            master->print_error("\"%s\" is an illegal option for swlspres\n", swlspres.c_str());
        }
    }

    if (swls == "1")
    {
        nerror += inputin->get_item(&swtimedep_ls,   "force", "swtimedep_ls",   "", "0");
        nerror += inputin->get_list(&lslist,         "force", "lslist",         "");
        nerror += inputin->get_list(&timedeplist_ls, "force", "timedeplist_ls", "");
    }
    else if (swls != "0")
    {
        ++nerror;
        master->print_error("\"%s\" is an illegal option for swls\n", swls.c_str());
    }

    if (swwls == "1")
    {
        nerror += inputin->get_item(&swtimedep_wls, "force", "swtimedep_wls", "", "0");
        fields->set_calc_mean_profs(true);
    }
    else if (swwls != "0")
    {
        ++nerror;
        master->print_error("\"%s\" is an illegal option for swwls\n", swwls.c_str());
    }

    if (swnudge == "1")
    {
        nerror += inputin->get_list(&nudgelist,         "force", "nudgelist", "");
        nerror += inputin->get_item(&swtimedep_nudge,   "force", "swtimedep_nudge",   "", "0");
        nerror += inputin->get_list(&timedeplist_nudge, "force", "timedeplist_nudge", "");
        fields->set_calc_mean_profs(true);
    }
    else if (swnudge != "0")
    {
        ++nerror;
        master->print_error("\"%s\" is an illegal option for swnudge\n", swnudge.c_str());
    }

    if (nerror)
        throw 1;
}

Force::~Force()
{
    delete[] ug;
    delete[] vg;
    delete[] wls;
    delete[] timedepdata_wls;
    delete[] nudge_factor;

    if (swls == "1")
    {
        for (std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
            delete[] lsprofs[*it];

        // Clean up time dependent data
        for (std::map<std::string, double*>::const_iterator it=timedepdata_ls.begin(); it!=timedepdata_ls.end(); ++it)
            delete[] it->second;
    }

    if (swnudge == "1")
    {
        for (std::vector<std::string>::const_iterator it=nudgelist.begin(); it!=nudgelist.end(); ++it)
            delete[] nudgeprofs[*it];

        // Clean up time dependent data
        for (std::map<std::string, double*>::const_iterator it=timedepdata_nudge.begin(); it!=timedepdata_nudge.end(); ++it)
            delete[] it->second;
    }
        // Clean up time dependent data
    for (std::map<std::string, double*>::const_iterator it=timedepdata_geo.begin(); it!=timedepdata_geo.end(); ++it)
        delete[] it->second;


    #ifdef USECUDA
    clear_device();
    #endif
}

void Force::init()
{
    if (swlspres == "geo")
    {
        ug = new double[grid->kcells];
        vg = new double[grid->kcells];
    }

    if (swls == "1")
    {
        for (std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
            lsprofs[*it] = new double[grid->kcells];
    }

    if (swwls == "1")
        wls = new double[grid->kcells];

    if (swnudge == "1")
    {
        nudge_factor = new double[grid->kcells];
        for (std::vector<std::string>::const_iterator it=nudgelist.begin(); it!=nudgelist.end(); ++it)
            nudgeprofs[*it] = new double[grid->kcells];
    }
}

void Force::create(Input *inputin)
{
    int nerror = 0;

    if (swlspres == "geo")
    {
        nerror += inputin->get_prof(&ug[grid->kstart], "ug", grid->kmax);
        nerror += inputin->get_prof(&vg[grid->kstart], "vg", grid->kmax);

        if (swtimedep_geo == "1")
        {
            // Read the time dependent geostrophic wind components
            nerror += model->input->get_time_prof(&timedepdata_geo["ug"], &timedeptime_geo["ug"], "ug", grid->kmax);
            nerror += model->input->get_time_prof(&timedepdata_geo["vg"], &timedeptime_geo["vg"], "vg", grid->kmax);
        }
    }

    if (swls == "1")
    {
        // check whether the fields in the list exist in the prognostic fields
        for (std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
            if (!fields->ap.count(*it))
            {
                master->print_error("field %s in [force][lslist] is illegal\n", it->c_str());
                ++nerror;
            }

        // read the large scale sources, which are the variable names with a "ls" suffix
        for (std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
            nerror += inputin->get_prof(&lsprofs[*it][grid->kstart], *it+"ls", grid->kmax);

        // Process the time dependent data
        if (swtimedep_ls == "1")
            nerror += create_timedep(timedepdata_ls, timedeptime_ls, timedeplist_ls, lslist, "ls");
    }

    if (swnudge == "1")
    {
        // Get profile with nudging factor as function of height
        nerror += inputin->get_prof(&nudge_factor[grid->kstart], "nudgefac", grid->kmax);

        // check whether the fields in the list exist in the prognostic fields
        for (std::vector<std::string>::const_iterator it=nudgelist.begin(); it!=nudgelist.end(); ++it)
            if (!fields->ap.count(*it))
            {
                master->print_error("field %s in [force][nudgelist] is illegal\n", it->c_str());
                ++nerror;
            }

        // read the large scale sources, which are the variable names with a "nudge" suffix
        for (std::vector<std::string>::const_iterator it=nudgelist.begin(); it!=nudgelist.end(); ++it)
            nerror += inputin->get_prof(&nudgeprofs[*it][grid->kstart], *it+"nudge", grid->kmax);

        // Process the time dependent data
        if (swtimedep_nudge == "1")
        {
            nerror += create_timedep(timedepdata_nudge, timedeptime_nudge, timedeplist_nudge, nudgelist, "nudge");
        }
    }

    // Get the large scale vertical velocity from the input
    if (swwls == "1")
    {
        nerror += inputin->get_prof(&wls[grid->kstart], "wls", grid->kmax);

        if (swtimedep_wls == "1")
            nerror += model->input->get_time_prof(&timedepdata_wls, &timedeptime_wls, "wls", grid->kmax);
    }

    if (nerror)
        throw 1;
}

int Force::create_timedep(std::map<std::string, double*>& data, std::map<std::string, std::vector<double>>& time,
                           std::vector<std::string>& timedep_variables, std::vector<std::string> all_variables, std::string suffix)
{
    int nerror = 0;

    // Create temporary copy of the time dependent variables, to check which ones are used
    std::vector<std::string> tmplist = timedep_variables;

    // Loop over all variables which might be time dependent
    for (auto& it : all_variables)
    {
        // Check if the current variable is time dependent
        if (std::find(timedep_variables.begin(), timedep_variables.end(), it) != timedep_variables.end())
        {
            std::string name = it + suffix;

            // Get the time dependent data and time
            nerror += model->input->get_time_prof(&data[name], &time[name], name, grid->kmax);
            // Remove from tmplist
            std::vector<std::string>::iterator ittmp = std::find(tmplist.begin(), tmplist.end(), it);
            if (ittmp != tmplist.end())
                tmplist.erase(ittmp);
        }
    }

    // Print a warning for the non-supported time dependent variables, and remove them from list of time dependent variables
    for (auto& it : tmplist)
    {
        master->print_warning("%s is not supported (yet) as a time dependent parameter\n", it.c_str());
        timedep_variables.erase(std::remove(timedep_variables.begin(), timedep_variables.end(), it), timedep_variables.end());
    }

    return nerror;
}

#ifndef USECUDA
void Force::exec(double dt)
{
    if (swlspres == "uflux")
        calc_flux(fields->ut->data, fields->u->data, grid->dz, dt);

    else if (swlspres == "geo")
    {
        if (grid->swspatialorder == "2")
            calc_coriolis_2nd(fields->ut->data, fields->vt->data, fields->u->data, fields->v->data, ug, vg);
        else if (grid->swspatialorder == "4")
            calc_coriolis_4th(fields->ut->data, fields->vt->data, fields->u->data, fields->v->data, ug, vg);
    }

    if (swls == "1")
    {
        for (std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
            calc_large_scale_source(fields->at[*it]->data, lsprofs[*it]);
    }

    if (swwls == "1")
    {
        for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); ++it)
            advec_wls_2nd(it->second->data, fields->sp[it->first]->datamean, wls, grid->dzhi);
    }

    if (swnudge == "1")
    {
        for (std::vector<std::string>::const_iterator it=nudgelist.begin(); it!=nudgelist.end(); ++it)
            calc_nudging_tendency(fields->at[*it]->data, fields->ap[*it]->datamean, nudgeprofs[*it], nudge_factor);
    }
}
#endif

void Force::update_time_dependent()
{
    if (swtimedep_ls == "1")
    {
        #ifndef USECUDA
        update_time_dependent_profs(lsprofs, timedepdata_ls, timedeptime_ls, "ls");
        #elif USECUDA
        update_time_dependent_profs(lsprofs_g, timedepdata_ls_g, timedeptime_ls, "ls");
        #endif
    }

    if (swtimedep_nudge == "1")
    {
        #ifndef USECUDA
        update_time_dependent_profs(nudgeprofs, timedepdata_nudge, timedeptime_nudge, "nudge");
        #elif USECUDA
        update_time_dependent_profs(nudgeprofs_g, timedepdata_nudge_g, timedeptime_nudge, "nudge");
        #endif
    }

    if (swtimedep_geo == "1")
    {
        #ifndef USECUDA
        update_time_dependent_prof(ug, timedepdata_geo["ug"], timedeptime_geo["ug"]);
        update_time_dependent_prof(vg, timedepdata_geo["vg"], timedeptime_geo["vg"]);
        #elif USECUDA
        update_time_dependent_prof(ug_g, timedepdata_geo_g["ug"], timedeptime_geo["ug"]);
        update_time_dependent_prof(vg_g, timedepdata_geo_g["vg"], timedeptime_geo["vg"]);
        #endif
    }

    if (swtimedep_wls == "1")
    {
        #ifndef USECUDA
        update_time_dependent_prof(wls, timedepdata_wls, timedeptime_wls);
        #elif USECUDA
        update_time_dependent_prof(wls_g, timedepdata_wls_g, timedeptime_wls);
        #endif
    }
}

#ifndef USECUDA
void Force::update_time_dependent_profs(std::map<std::string, double*>& profiles, std::map<std::string, double*> time_profiles,
                                        std::map<std::string, std::vector<double>> times, std::string suffix)
{
    const int kk = grid->kmax;
    const int kgc = grid->kgc;

    // Loop over all profiles which might be time dependent
    for (auto& it : profiles)
    {
        std::string name = it.first + suffix;

        // Check if they have time dependent data
        if (time_profiles.find(name) != time_profiles.end())
        {
            // Get/calculate the interpolation indexes/factors
            int index0, index1;
            double fac0, fac1;

            model->timeloop->get_interpolation_factors(index0, index1, fac0, fac1, times[name]);

            // Calculate the new vertical profile
            for (int k=0; k<grid->kmax; ++k)
                it.second[k+kgc] = fac0 * time_profiles[name][index0*kk+k] + fac1 * time_profiles[name][index1*kk+k];
        }
    }
}
#endif
#ifndef USECUDA
void Force::update_time_dependent_prof(double* const restrict prof, const double* const restrict data,
                                       const std::vector<double> times)
{
    const int kk = grid->kmax;
    const int kgc = grid->kgc;

    // Get/calculate the interpolation indexes/factors
    int index0, index1;
    double fac0, fac1;
    model->timeloop->get_interpolation_factors(index0, index1, fac0, fac1, times);

    // Calculate the new vertical profile
    for (int k=0; k<grid->kmax; ++k)
        prof[k+kgc] = fac0 * data[index0*kk+k] + fac1 * data[index1*kk+k];
}
#endif
void Force::calc_flux(double* const restrict ut, const double* const restrict u,
                      const double* const restrict dz, const double dt)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double ugrid = grid->utrans;

    double uavg  = 0.;
    double utavg = 0.;

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                uavg  += u [ijk]*dz[k];
                utavg += ut[ijk]*dz[k];
            }

    grid->get_sum(&uavg);
    grid->get_sum(&utavg);

    uavg  = uavg  / (grid->itot*grid->jtot*grid->zsize);
    utavg = utavg / (grid->itot*grid->jtot*grid->zsize);

    double fbody;
    fbody = (uflux - uavg - ugrid) / dt - utavg;

    for (int n=0; n<grid->ncells; n++)
        ut[n] += fbody;
}

void Force::calc_coriolis_2nd(double* const restrict ut, double* const restrict vt,
                              const double* const restrict u , const double* const restrict v ,
                              const double* const restrict ug, const double* const restrict vg)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double ugrid = grid->utrans;
    const double vgrid = grid->vtrans;

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                ut[ijk] += fc * (0.25*(v[ijk-ii] + v[ijk] + v[ijk-ii+jj] + v[ijk+jj]) + vgrid - vg[k]);
            }

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                vt[ijk] -= fc * (0.25*(u[ijk-jj] + u[ijk] + u[ijk+ii-jj] + u[ijk+ii]) + ugrid - ug[k]);
            }
}

void Force::calc_coriolis_4th(double* const restrict ut, double* const restrict vt,
                              const double* const restrict u , const double* const restrict v ,
                              const double* const restrict ug, const double* const restrict vg)
{
    using namespace Finite_difference::O4;

    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*grid->icells;
    const int jj2 = 2*grid->icells;
    const int kk1 = 1*grid->ijcells;

    const double ugrid = grid->utrans;
    const double vgrid = grid->vtrans;

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                ut[ijk] += fc * ( ( ci0*(ci0*v[ijk-ii2-jj1] + ci1*v[ijk-ii1-jj1] + ci2*v[ijk-jj1] + ci3*v[ijk+ii1-jj1])
                                  + ci1*(ci0*v[ijk-ii2    ] + ci1*v[ijk-ii1    ] + ci2*v[ijk    ] + ci3*v[ijk+ii1    ])
                                  + ci2*(ci0*v[ijk-ii2+jj1] + ci1*v[ijk-ii1+jj1] + ci2*v[ijk+jj1] + ci3*v[ijk+ii1+jj1])
                                  + ci3*(ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii1+jj2] + ci2*v[ijk+jj2] + ci3*v[ijk+ii1+jj2]) )
                                + vgrid - vg[k] );
            }

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                vt[ijk] -= fc * ( ( ci0*(ci0*u[ijk-ii1-jj2] + ci1*u[ijk-jj2] + ci2*u[ijk+ii1-jj2] + ci3*u[ijk+ii2-jj2])
                                  + ci1*(ci0*u[ijk-ii1-jj1] + ci1*u[ijk-jj1] + ci2*u[ijk+ii1-jj1] + ci3*u[ijk+ii2-jj1])
                                  + ci2*(ci0*u[ijk-ii1    ] + ci1*u[ijk    ] + ci2*u[ijk+ii1    ] + ci3*u[ijk+ii2    ])
                                  + ci3*(ci0*u[ijk-ii1+jj1] + ci1*u[ijk+jj1] + ci2*u[ijk+ii1+jj1] + ci3*u[ijk+ii2+jj1]) )
                                + ugrid - ug[k]);
            }
}

void Force::calc_large_scale_source(double* const restrict st, const double* const restrict sls)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                st[ijk] += sls[k];
            }
}

void Force::calc_nudging_tendency(double* const restrict fldtend, const double* const restrict fldmean,
                                  const double* const restrict ref, const double* const restrict factor)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=grid->kstart; k<grid->kend; ++k)
    {
        const double tend = -factor[k] * (fldmean[k] - ref[k]);
        for (int j=grid->jstart; j<grid->jend; ++j)
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                fldtend[ijk] += tend;
            }
    }
}

void Force::advec_wls_2nd(double* const restrict st, const double* const restrict s,
                          const double* const restrict wls, const double* const dzhi)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    // use an upwind differentiation
    for (int k=grid->kstart; k<grid->kend; ++k)
    {
        if (wls[k] > 0.)
        {
            for (int j=grid->jstart; j<grid->jend; ++j)
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] -=  wls[k] * (s[k]-s[k-1])*dzhi[k];
                }
        }
        else
        {
            for (int j=grid->jstart; j<grid->jend; ++j)
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] -=  wls[k] * (s[k+1]-s[k])*dzhi[k+1];
                }
        }
    }
}
