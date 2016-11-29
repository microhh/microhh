/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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
#include <sstream>
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "model.h"
#include "stats.h"
#include "diff_smag2.h"
#include "master.h"
#include "cross.h"
#include "dump.h"

#include "thermo_dry.h"
#include "thermo_moist_functions.h"  // For Exner function

using Finite_difference::O2::interp2;
using Finite_difference::O4::interp4;
using namespace Constants;
using namespace Thermo_moist_functions;

Thermo_dry::Thermo_dry(Model *modelin, Input *inputin) : Thermo(modelin, inputin)
{
    swthermo = "dry";

    thref     = 0;
    threfh    = 0;
    pref      = 0;
    prefh     = 0;
    exnref    = 0;
    exnrefh   = 0;

    thref_g   = 0;
    threfh_g  = 0;
    pref_g    = 0;
    prefh_g   = 0;
    exnref_g  = 0;
    exnrefh_g = 0;

    int nerror = 0;

    fields->init_prognostic_field("th", "Potential Temperature", "K");

    nerror += inputin->get_item(&fields->sp["th"]->visc, "fields", "svisc", "th");

    // Get base state option (boussinesq or anelastic)
    nerror += inputin->get_item(&swbasestate, "thermo", "swbasestate", "", "");

    if (!(swbasestate == "boussinesq" || swbasestate == "anelastic"))
    {
        master->print_error("\"%s\" is an illegal value for swbasestate\n", swbasestate.c_str());
        throw 1;
    }

    if (grid->swspatialorder == "4" && swbasestate == "anelastic")
    {
        master->print_error("Anelastic mode is not supported for swspatialorder=4\n");
        throw 1;
    }

    // Remove the data from the input that is not used, to avoid warnings.
    if (master->mode == "init")
    {
        inputin->flag_as_used("thermo", "thref0");
        inputin->flag_as_used("thermo", "pbot");
    }

    if (nerror)
        throw 1;
}

Thermo_dry::~Thermo_dry()
{
    delete[] thref;
    delete[] threfh;
    delete[] pref;
    delete[] prefh;
    delete[] exnref;
    delete[] exnrefh;

#ifdef USECUDA
    clear_device();
#endif
}

void Thermo_dry::init()
{
    // copy pointers
    stats = model->stats;

    // fields for Boussinesq and anelastic solver
    thref   = new double[grid->kcells];
    threfh  = new double[grid->kcells];
    pref    = new double[grid->kcells];
    prefh   = new double[grid->kcells];
    exnref  = new double[grid->kcells];
    exnrefh = new double[grid->kcells];

    init_cross();
    init_dump(); 
}

void Thermo_dry::create(Input *inputin)
{
    /* Setup base state: 
       For anelastic setup, calculate reference density and temperature from input sounding
       For boussinesq, reference density and temperature are fixed */
    if (swbasestate == "anelastic")
    {
        if (inputin->get_item(&pbot, "thermo", "pbot", ""))
            throw 1;
        if (inputin->get_prof(&thref[grid->kstart], "th", grid->kmax))
            throw 1;
        calc_base_state(fields->rhoref, fields->rhorefh, pref, prefh, exnref, exnrefh, thref, threfh, pbot);
    }
    else
    {
        if (inputin->get_item(&thref0, "thermo", "thref0", ""))
            throw 1;

        // Set entire column to reference value. Density is already initialized at 1.0 in fields.cxx
        for (int k=0; k<grid->kcells; ++k)
        {
            thref[k]  = thref0;
            threfh[k] = thref0;
        }
    }

    init_stat();
}

#ifndef USECUDA
void Thermo_dry::exec()
{
    if (grid->swspatialorder== "2")
        calc_buoyancy_tend_2nd(fields->wt->data, fields->sp["th"]->data, threfh);
    else if (grid->swspatialorder == "4")
        calc_buoyancy_tend_4th(fields->wt->data, fields->sp["th"]->data, threfh);
}
#endif

unsigned long Thermo_dry::get_time_limit(unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}

void Thermo_dry::exec_stats(Mask *m)
{
    const double NoOffset = 0.;

    // calculate the buoyancy and its surface flux for the profiles
    calc_buoyancy(fields->atmp["tmp1"]->data, fields->sp["th"]->data, thref);
    calc_buoyancy_fluxbot(fields->atmp["tmp1"]->datafluxbot, fields->sp["th"]->datafluxbot, threfh);

    // define the location
    const int sloc[] = {0,0,0};

    // calculate the mean
    stats->calc_mean(m->profs["b"].data, fields->atmp["tmp1"]->data, NoOffset, sloc,
                     fields->atmp["tmp3"]->data, stats->nmask);

    // calculate the moments
    for (int n=2; n<5; ++n)
    {
        std::stringstream ss;
        ss << n;
        std::string sn = ss.str();
        stats->calc_moment(fields->atmp["tmp1"]->data, m->profs["b"].data, m->profs["b"+sn].data, n, sloc,
                           fields->atmp["tmp3"]->data, stats->nmask);
    }

    // calculate the gradients
    if (grid->swspatialorder == "2")
        stats->calc_grad_2nd(fields->atmp["tmp1"]->data, m->profs["bgrad"].data, grid->dzhi, sloc,
                             fields->atmp["tmp4"]->data, stats->nmaskh);
    else if (grid->swspatialorder == "4")
        stats->calc_grad_4th(fields->atmp["tmp1"]->data, m->profs["bgrad"].data, grid->dzhi4, sloc,
                             fields->atmp["tmp4"]->data, stats->nmaskh);

    // calculate turbulent fluxes
    if (grid->swspatialorder == "2")
        stats->calc_flux_2nd(fields->atmp["tmp1"]->data, m->profs["b"].data, fields->w->data, m->profs["w"].data,
                             m->profs["bw"].data, fields->atmp["tmp2"]->data, sloc,
                             fields->atmp["tmp4"]->data, stats->nmaskh);
    else if (grid->swspatialorder == "4")
        stats->calc_flux_4th(fields->atmp["tmp1"]->data, fields->w->data, m->profs["bw"].data, fields->atmp["tmp2"]->data, sloc,
                             fields->atmp["tmp4"]->data, stats->nmaskh);

    // calculate diffusive fluxes
    if (grid->swspatialorder == "2")
    {
        if (model->diff->get_switch() == "smag2")
        {
            Diff_smag_2* diffptr = static_cast<Diff_smag_2*>(model->diff);
            stats->calc_diff_2nd(fields->atmp["tmp1"]->data, fields->w->data, fields->sd["evisc"]->data,
                                 m->profs["bdiff"].data, grid->dzhi,
                                 fields->atmp["tmp1"]->datafluxbot, fields->atmp["tmp1"]->datafluxtop, diffptr->tPr, sloc,
                                 fields->atmp["tmp4"]->data, stats->nmaskh);
        }
        else
            stats->calc_diff_2nd(fields->atmp["tmp1"]->data, m->profs["bdiff"].data, grid->dzhi, fields->sp["th"]->visc, sloc,
                                 fields->atmp["tmp4"]->data, stats->nmaskh);
    }
    else if (grid->swspatialorder == "4")
    {
        stats->calc_diff_4th(fields->atmp["tmp1"]->data, m->profs["bdiff"].data, grid->dzhi4, fields->sp["th"]->visc, sloc,
                             fields->atmp["tmp4"]->data, stats->nmaskh);
    }

    // calculate the total fluxes
    stats->add_fluxes(m->profs["bflux"].data, m->profs["bw"].data, m->profs["bdiff"].data);

    // calculate the sorted buoyancy profile
    //stats->calc_sorted_prof(fields->sd["tmp1"]->data, fields->sd["tmp2"]->data, m->profs["bsort"].data);
}

void Thermo_dry::exec_cross()
{
    int nerror = 0;

    Cross* cross = model->cross;

    // With one additional temp field, we wouldn't have to re-calculate the ql or b field for simple,lngrad,path, etc.
    for (std::vector<std::string>::iterator it=crosslist.begin(); it<crosslist.end(); ++it)
    {
        /* BvS: for now, don't call getThermoField() or getBuoyancySurf(), but directly the function itself. With CUDA enabled, 
           statistics etc. is done on the host, while getThermoField() is executed on the GPU */ 

        if (*it == "b")
        {
            //getThermoField(fields->s["tmp1"], fields->s["tmp2"], *it);
            calc_buoyancy(fields->atmp["tmp1"]->data, fields->sp["th"]->data, thref);
            nerror += cross->cross_simple(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, *it);
        }
        else if (*it == "blngrad")
        {
            //getThermoField(fields->s["tmp1"], fields->s["tmp2"], "b");
            calc_buoyancy(fields->atmp["tmp1"]->data, fields->sp["th"]->data, thref);
            // Note: tmp1 twice used as argument -> overwritten in crosspath()
            nerror += cross->cross_lngrad(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, grid->dzi4, *it);
        }
        else if (*it == "bbot" or *it == "bfluxbot")
        {
            //getBuoyancySurf(fields->s["tmp1"]);
            calc_buoyancy_bot(fields->atmp["tmp1"]->data, fields->atmp["tmp1"]->databot, fields->sp["th"]->data, fields->sp["th"]->databot, thref, threfh);
            calc_buoyancy_fluxbot(fields->atmp["tmp1"]->datafluxbot, fields->sp["th"]->datafluxbot, threfh);

            if (*it == "bbot")
                nerror += cross->cross_plane(fields->atmp["tmp1"]->databot, fields->atmp["tmp1"]->data, "bbot");
            else if (*it == "bfluxbot")
                nerror += cross->cross_plane(fields->atmp["tmp1"]->datafluxbot, fields->atmp["tmp1"]->data, "bfluxbot");
        }
    }

    if (nerror)
        throw 1;
}

void Thermo_dry::exec_dump()
{
    for (std::vector<std::string>::const_iterator it=dumplist.begin(); it<dumplist.end(); ++it)
    {
        // TODO BvS restore getThermoField(), the combination of checkThermoField with getThermoField is more elegant... 
        if (*it == "b")
            calc_buoyancy(fields->atmp["tmp2"]->data, fields->sp["th"]->data, thref);
        else
            throw 1;

        model->dump->save_dump(fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, *it);
    }
}

bool Thermo_dry::check_field_exists(std::string name)
{
    if (name == "b")
        return true;
    else
        return false;
}

#ifndef USECUDA
void Thermo_dry::get_thermo_field(Field3d *fld, Field3d *tmp, std::string name, bool cyclic)
{
    if (name == "b")
        calc_buoyancy(fld->data, fields->sp["th"]->data, thref);
    else if (name == "N2")
        calc_N2(fld->data, fields->sp["th"]->data, grid->dzi, thref);
    else
        throw 1;

    if (cyclic)
        grid->boundary_cyclic(fld->data);
}
#endif

#ifndef USECUDA
void Thermo_dry::get_buoyancy_fluxbot(Field3d* bfield)
{
    calc_buoyancy_fluxbot(bfield->datafluxbot, fields->sp["th"]->datafluxbot, threfh);
}
#endif

#ifndef USECUDA
void Thermo_dry::get_buoyancy_surf(Field3d *bfield)
{
    calc_buoyancy_bot(bfield->data, bfield->databot,
            fields->sp["th"]->data, fields->sp["th"]->databot, thref, threfh);
    calc_buoyancy_fluxbot(bfield->datafluxbot, fields->sp["th"]->datafluxbot, threfh);
}
#endif

void Thermo_dry::get_prog_vars(std::vector<std::string>* list)
{
    list->push_back("th");
}

double Thermo_dry::get_buoyancy_diffusivity()
{
    // Use the diffusivity from theta
    return fields->sp["th"]->visc; 
}

void Thermo_dry::calc_buoyancy(double* restrict b, double* restrict th, double* restrict thref)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=0; k<grid->kcells; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                b[ijk] = grav/thref[k] * (th[ijk] - thref[k]);
            }
}

void Thermo_dry::calc_N2(double* restrict N2, double* restrict th, double* restrict dzi, double* restrict thref)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                N2[ijk] = grav/thref[k]*0.5*(th[ijk+kk] - th[ijk-kk])*dzi[k];
            }
}

void Thermo_dry::calc_buoyancy_bot(double* restrict b , double* restrict bbot,
                                   double* restrict th, double* restrict thbot,
                                   double* restrict thref, double* restrict threfh)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
        for (int i=0; i<grid->icells; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            bbot[ij] = grav/threfh[kstart] * (thbot[ij] - threfh[kstart]);
            b[ijk]   = grav/thref [kstart] * (th[ijk]   - thref [kstart]);
        }
}

void Thermo_dry::calc_buoyancy_fluxbot(double* restrict bfluxbot, double* restrict thfluxbot, double* restrict threfh)
{
    const int jj = grid->icells;
    const int kstart = grid->kstart;

    for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
        for (int i=0; i<grid->icells; ++i)
        {
            const int ij = i + j*jj;
            bfluxbot[ij] = grav/threfh[kstart]*thfluxbot[ij];
        }
}

void Thermo_dry::calc_buoyancy_tend_2nd(double* restrict wt, double* restrict th, double* restrict threfh)
{
    using namespace Finite_difference::O2;

    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=grid->kstart+1; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                wt[ijk] += grav/threfh[k] * (interp2(th[ijk-kk], th[ijk]) - threfh[k]);
            }
}

void Thermo_dry::calc_buoyancy_tend_4th(double* restrict wt, double* restrict th, double* restrict threfh)
{
    const int jj  = grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    for (int k=grid->kstart+1; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk1;
                wt[ijk] += grav/threfh[k] * (interp4(th[ijk-kk2], th[ijk-kk1], th[ijk], th[ijk+kk1]) - threfh[k]);
            }
}

// Initialize the base state for the anelastic solver
void Thermo_dry::calc_base_state(double* restrict rhoref, double* restrict rhorefh,
                                 double* restrict pref,   double* restrict prefh,
                                 double* restrict exnref, double* restrict exnrefh,
                                 double* restrict thref,  double* restrict threfh,
                                 const double pbot)
{
    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    // extrapolate the input sounding to get the bottom value
    threfh[kstart] = thref[kstart] - grid->z[kstart]*(thref[kstart+1]-thref[kstart])*grid->dzhi[kstart+1];

    // extrapolate the input sounding to get the top value
    threfh[kend] = thref[kend-1] + (grid->zh[kend]-grid->z[kend-1])*(thref[kend-1]-thref[kend-2])*grid->dzhi[kend-1];

    // set the ghost cells for the reference potential temperature
    thref[kstart-1] = 2.*threfh[kstart] - thref[kstart];
    thref[kend]     = 2.*threfh[kend]   - thref[kend-1];

    // interpolate the input sounding to half levels
    for (int k=kstart+1; k<kend; ++k)
        threfh[k] = 0.5*(thref[k-1] + thref[k]);

    // Calculate pressure
    prefh[kstart] = pbot;
    pref [kstart] = pbot * std::exp(-grav * grid->z[kstart] / (Rd * threfh[kstart] * exner(prefh[kstart])));
    for (int k=kstart+1; k<kend+1; ++k)
    {
        prefh[k] = prefh[k-1] * std::exp(-grav * grid->dz[k-1] / (Rd * thref[k-1] * exner(pref[k-1])));
        pref [k] = pref [k-1] * std::exp(-grav * grid->dzh[k ] / (Rd * threfh[k ] * exner(prefh[k ])));
    }
    pref[kstart-1] = 2.*prefh[kstart] - pref[kstart];

    // Calculate density and exner
    for (int k=0; k<grid->kcells; ++k)
    {
        exnref[k]  = exner(pref[k] );
        exnrefh[k] = exner(prefh[k]);
        rhoref[k]  = pref[k]  / (Rd * thref[k]  * exnref[k] );
        rhorefh[k] = prefh[k] / (Rd * threfh[k] * exnrefh[k]);
    }
}

void Thermo_dry::init_stat()
{
    if (stats->get_switch() == "1")
    {
        // Add base state profiles to statistics
        stats->add_fixed_prof("rhoref",  "Full level basic state density",  "kg m-3", "z",  fields->rhoref);
        stats->add_fixed_prof("rhorefh", "Half level basic state density",  "kg m-3", "zh", fields->rhorefh);
        stats->add_fixed_prof("thref",   "Full level basic state potential temperature", "K", "z", thref);
        stats->add_fixed_prof("threfh",  "Half level basic state potential temperature", "K", "zh",thref);
        if (swbasestate == "anelastic")
        {
            stats->add_fixed_prof("ph",    "Full level hydrostatic pressure", "Pa",     "z",  pref);
            stats->add_fixed_prof("phh",   "Half level hydrostatic pressure", "Pa",     "zh", prefh);
        }

        stats->add_prof("b", "Buoyancy", "m s-2", "z");
        for (int n=2; n<5; ++n)
        {
            std::stringstream ss;
            ss << n;
            std::string sn = ss.str();
            stats->add_prof("b"+sn, "Moment " +sn+" of the buoyancy", "(m s-2)"+sn,"z");
        }

        stats->add_prof("bgrad", "Gradient of the buoyancy", "s-2", "zh");
        stats->add_prof("bw"   , "Turbulent flux of the buoyancy", "m2 s-3", "zh");
        stats->add_prof("bdiff", "usive flux of the buoyancy", "m2 s-3", "zh");
        stats->add_prof("bflux", "Total flux of the buoyancy", "m2 s-3", "zh");

        stats->add_prof("bsort", "Sorted buoyancy", "m s-2", "z");
    }
}

void Thermo_dry::init_cross()
{
    if (model->cross->get_switch() == "1")
    {
        // Populate list with allowed cross-section variables
        allowedcrossvars.push_back("b");
        allowedcrossvars.push_back("bbot");
        allowedcrossvars.push_back("bfluxbot");
        if (grid->swspatialorder == "4")
            allowedcrossvars.push_back("blngrad");

        // Get global cross-list from cross.cxx
        std::vector<std::string> *crosslist_global = model->cross->get_crosslist(); 

        // Check input list of cross variables (crosslist)
        std::vector<std::string>::iterator it=crosslist_global->begin();
        while (it != crosslist_global->end())
        {
            if (std::count(allowedcrossvars.begin(),allowedcrossvars.end(),*it))
            {
                // Remove variable from global list, put in local list
                crosslist.push_back(*it);
                crosslist_global->erase(it); // erase() returns iterator of next element..
            }
            else
                ++it;
        }

        // Sort crosslist to group ql and b variables
        std::sort(crosslist.begin(),crosslist.end());
    }
}

void Thermo_dry::init_dump()
{
    if (model->dump->get_switch() == "1")
    {
        // Get global cross-list from cross.cxx
        std::vector<std::string>* dumplist_global = model->dump->get_dumplist(); 

        // Check if fields in dumplist are retrievable thermo fields
        std::vector<std::string>::iterator dumpvar=dumplist_global->begin();
        while (dumpvar != dumplist_global->end())
        {
            if (check_field_exists(*dumpvar))
            {
                // Remove variable from global list, put in local list
                dumplist.push_back(*dumpvar);
                dumplist_global->erase(dumpvar); // erase() returns iterator of next element..
            }
            else
                ++dumpvar;
        }
    }
}
