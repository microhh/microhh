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
#include <sstream>
#include <algorithm>
#include <netcdf>
#include "grid.h"
#include "fields.h"
#include "thermo_vapor.h"
#include "diff_smag2.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "model.h"
#include "stats.h"
#include "master.h"
#include "cross.h"
#include "dump.h"
#include "column.h"
#include "thermo_moist_functions.h"
#include "timeloop.h"

using Finite_difference::O2::interp2;
using Finite_difference::O4::interp4;
using namespace Constants;
using namespace Thermo_moist_functions;


Thermo_vapor::Thermo_vapor(Model* modelin, Input* inputin) : Thermo(modelin, inputin)
{
    swthermo = "vapor";

    thl0 = 0;
    qt0  = 0;

    thvref  = 0;
    thvrefh = 0;
    exnref  = 0;
    exnrefh = 0;
    pref    = 0;
    prefh   = 0;

    thvref_g  = 0;
    thvrefh_g = 0;
    exnref_g  = 0;
    exnrefh_g = 0;
    pref_g    = 0;
    prefh_g   = 0;

    int nerror = 0;

    // Option to overrule the prognostic variable
    nerror += inputin->get_item(&thvar, "thermo", "progvar", "", "thl");  // defaults to thl

    // Initialize the prognostic fields
    fields->init_prognostic_field(thvar, "Liquid water potential temperature", "K");
    fields->init_prognostic_field("qt", "Total water mixing ratio", "kg kg-1");

    // Get the diffusivities of temperature and moisture
    nerror += inputin->get_item(&fields->sp[thvar]->visc, "fields", "svisc", thvar );
    nerror += inputin->get_item(&fields->sp["qt"]->visc, "fields", "svisc", "qt");

    // Test if the diffusivities of theta and qt are equal, else throw error
    if (fields->sp[thvar]->visc != fields->sp["qt"]->visc)
    {
        master->print_error("The diffusivities of temperature and moisture must be equal\n");
        throw 1;
    }

    nerror += inputin->get_item(&pbot, "thermo", "pbot", "");

    // Get base state option (boussinesq or anelastic)
    nerror += inputin->get_item(&swbasestate, "thermo", "swbasestate", "", "");

    if (!(swbasestate == "boussinesq" || swbasestate == "anelastic"))
    {
        master->print_error("\"%s\" is an illegal value for swbasestate\n", swbasestate.c_str());
        throw 1;
    }

    // Note BvS: the base state calculation was never tested/supported for 4th order DNS,
    // and has therefor been removed for now.
    if (grid->swspatialorder == "4") // && swbasestate == "anelastic")
    {
        //master->print_error("Anelastic mode is not supported for swspatialorder=4\n");
        master->print_error("Moist thermodynamics are not supported for swspatialorder=4\n");
        throw 1;
    }

    // BvS test for updating hydrostatic prssure during run
    // swupdate..=0 -> initial base state pressure used in saturation calculation
    // swupdate..=1 -> base state pressure updated before saturation calculation
    nerror += inputin->get_item(&swupdatebasestate, "thermo", "swupdatebasestate", "");
    
    // Time variable surface pressure
    nerror += inputin->get_item(&swtimedep_pbot, "thermo", "swtimedep_pbot", "", 0);

    // Remove the data from the input that is not used, to avoid warnings.
    if (master->mode == "init")
    {
        inputin->flag_as_used("thermo", "thvref0");
        inputin->flag_as_used("thermo", "pbot");
    }

    if (nerror)
        throw 1;
}

Thermo_vapor::~Thermo_vapor()
{
    delete[] thl0;
    delete[] qt0;
    delete[] thvref;
    delete[] thvrefh;
    delete[] exnref;
    delete[] exnrefh;
    delete[] pref;
    delete[] prefh;

    #ifdef USECUDA
    clear_device();
    #endif
}

void Thermo_vapor::init()
{
    stats = model->stats;

    thl0    = new double[grid->kcells];
    qt0     = new double[grid->kcells];
    thvref  = new double[grid->kcells];
    thvrefh = new double[grid->kcells];
    exnref  = new double[grid->kcells];
    exnrefh = new double[grid->kcells];
    pref    = new double[grid->kcells];
    prefh   = new double[grid->kcells];

    for (int k=0; k<grid->kcells; ++k)
    {
        thl0   [k] = 0.;
        qt0    [k] = 0.;
        thvref [k] = 0.;
        thvrefh[k] = 0.;
        exnref [k] = 0.;
        exnrefh[k] = 0.;
        pref   [k] = 0.;
        prefh  [k] = 0.;
    }

    init_cross();
    init_dump();
}

void Thermo_vapor::create(Input* inputin)
{
    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    // Enable automated calculation of horizontally averaged fields
    if (swupdatebasestate)
        fields->set_calc_mean_profs(true);

    // Calculate the base state profiles. With swupdatebasestate=1, these profiles are updated on every iteration.
    // 1. Take the initial profile as the reference
    if (inputin->get_prof(&thl0[grid->kstart], thvar, grid->kmax))
        throw 1;
    if (inputin->get_prof(&qt0[grid->kstart], "qt", grid->kmax))
        throw 1;

    // 2. Calculate surface and model top values thl and qt
    double thl0s, qt0s, thl0t, qt0t;
    thl0s = thl0[kstart] - grid->z[kstart]*(thl0[kstart+1]-thl0[kstart])*grid->dzhi[kstart+1];
    qt0s  = qt0[kstart]  - grid->z[kstart]*(qt0[kstart+1] -qt0[kstart] )*grid->dzhi[kstart+1];
    thl0t = thl0[kend-1] + (grid->zh[kend]-grid->z[kend-1])*(thl0[kend-1]-thl0[kend-2])*grid->dzhi[kend-1];
    qt0t  = qt0[kend-1]  + (grid->zh[kend]-grid->z[kend-1])*(qt0[kend-1]- qt0[kend-2] )*grid->dzhi[kend-1];

    // 3. Set the ghost cells for the reference temperature and moisture
    thl0[kstart-1]  = 2.*thl0s - thl0[kstart];
    thl0[kend]      = 2.*thl0t - thl0[kend-1];
    qt0[kstart-1]   = 2.*qt0s  - qt0[kstart];
    qt0[kend]       = 2.*qt0t  - qt0[kend-1];

    // 4. Calculate the initial/reference base state
    calc_base_state(pref, prefh, fields->rhoref, fields->rhorefh, thvref, thvrefh, exnref, exnrefh, thl0, qt0);

    // 5. In Boussinesq mode, overwrite reference temperature and density
    if (swbasestate == "boussinesq")
    {
        if (inputin->get_item(&thvref0, "thermo", "thvref0", ""))
            throw 1;

        for (int k=0; k<grid->kcells; ++k)
        {
            fields->rhoref[k]  = 1.;
            fields->rhorefh[k] = 1.;
            thvref[k]          = thvref0;
            thvrefh[k]         = thvref0;
        }
    }

    // 6. Process the time dependent surface pressure
    if (swtimedep_pbot == 1)
    {
        const int nerror = inputin->get_time(&timedeppbot, &timedeptime, "pbot");
        if (nerror > 0)
            throw 1;
    }

    init_stat();
    init_column();
}

#ifndef USECUDA
void Thermo_vapor::exec()
{
    const int kk = grid->ijcells;
    const int kcells = grid->kcells;


    // Re-calculate hydrostatic pressure and exner, pass dummy as rhoref,thvref to prevent overwriting base state
    double *tmp2 = fields->atmp["tmp2"]->data;
    if (swupdatebasestate)
        calc_base_state(pref, prefh,
                        &tmp2[0*kcells], &tmp2[1*kcells], &tmp2[2*kcells], &tmp2[3*kcells],
                        exnref, exnrefh, fields->sp[thvar]->datamean, fields->sp["qt"]->datamean);

    // extend later for gravity vector not normal to surface
    if (grid->swspatialorder == "2")
    {
        calc_buoyancy_tend_2nd(fields->wt->data, fields->sp[thvar]->data, fields->sp["qt"]->data, prefh,
                               &fields->atmp["tmp2"]->data[0*kk], &fields->atmp["tmp2"]->data[1*kk],
                               thvrefh);
    }
    //else if (grid->swspatialorder == "4")
    //{
    //    calc_buoyancy_tend_4th(fields->wt->data, fields->sp[thvar]->data, fields->sp["qt"]->data, prefh,
    //                           &fields->atmp["tmp2"]->data[0*kk], &fields->atmp["tmp2"]->data[1*kk],
    //                           &fields->atmp["tmp2"]->data[2*kk],
    //                           thvrefh);
    //}

}
#endif

unsigned long Thermo_vapor::get_time_limit(unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}

void Thermo_vapor::exec_stats(Mask *m)
{
    const double NoOffset = 0.;

    // calc the buoyancy and its surface flux for the profiles
    calc_buoyancy(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, thvref);
    calc_buoyancy_fluxbot(fields->atmp["tmp1"]->datafluxbot, fields->sp[thvar]->databot, fields->sp[thvar]->datafluxbot, fields->sp["qt"]->databot, fields->sp["qt"]->datafluxbot, thvrefh);

    // define location
    const int sloc[] = {0,0,0};

    // mean
    stats->calc_mean(m->profs["b"].data, fields->atmp["tmp1"]->data, NoOffset, sloc,
                     fields->atmp["tmp3"]->data, stats->nmask);

    // moments
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
            Diff_smag_2 *diffptr = static_cast<Diff_smag_2 *>(model->diff);
            stats->calc_diff_2nd(fields->atmp["tmp1"]->data, fields->w->data, fields->sd["evisc"]->data,
                                 m->profs["bdiff"].data, grid->dzhi,
                                 fields->atmp["tmp1"]->datafluxbot, fields->atmp["tmp1"]->datafluxtop, diffptr->tPr, sloc,
                                 fields->atmp["tmp4"]->data, stats->nmaskh);
        }
        else
        {
            stats->calc_diff_2nd(fields->atmp["tmp1"]->data, m->profs["bdiff"].data, grid->dzhi, fields->sp[thvar]->visc, sloc,
                                 fields->atmp["tmp4"]->data, stats->nmaskh);
        }
    }
    else if (grid->swspatialorder == "4")
    {
        // take the diffusivity of temperature for that of buoyancy
        stats->calc_diff_4th(fields->atmp["tmp1"]->data, m->profs["bdiff"].data, grid->dzhi4, fields->sp[thvar]->visc, sloc,
                             fields->atmp["tmp4"]->data, stats->nmaskh);
    }

    // calculate the total fluxes
    stats->add_fluxes(m->profs["bflux"].data, m->profs["bw"].data, m->profs["bdiff"].data);
}


void Thermo_vapor::exec_column()
{
    const double NoOffset = 0.;

    // Buoyancy mean
    model->column->calc_column(model->column->profs["b"].data, fields->atmp["tmp1"]->data, NoOffset);
}

void Thermo_vapor::exec_cross(int iotime)
{
    int nerror = 0;

    Cross* cross = model->cross;

    // With one additional temp field, we wouldn't have to re-calculate the b field for simple,lngrad,path, etc.
    for (std::vector<std::string>::iterator it=crosslist.begin(); it<crosslist.end(); ++it)
    {
        /* BvS: for now, don't call getThermoField() or getBuoyancySurf(), but directly the function itself. With CUDA enabled,
           statistics etc. is done on the host, while getThermoField() is executed on the GPU */

        if (*it == "b")
        {
            calc_buoyancy(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, thvref);
            nerror += cross->cross_simple(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, *it, iotime);
        }
        else if (*it == "blngrad")
        {
            calc_buoyancy(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, thvref);
            // Note: tmp1 twice used as argument -> overwritten in crosspath()
            nerror += cross->cross_lngrad(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, grid->dzi4, *it, iotime);
        }
    }

    if (nerror)
        throw 1;
}

void Thermo_vapor::exec_dump(int iotime)
{
    for (std::vector<std::string>::const_iterator it=dumplist.begin(); it<dumplist.end(); ++it)
    {
        // TODO BvS restore getThermoField(), the combination of checkThermoField with getThermoField is more elegant...
        if (*it == "b")
            calc_buoyancy(fields->atmp["tmp2"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, thvref);
        else
            throw 1;

        model->dump->save_dump(fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, *it, iotime);
    }
}

bool Thermo_vapor::check_field_exists(const std::string name)
{
    if (name == "b")
        return true;
    else
        return false;
}

void Thermo_vapor::update_time_dependent()
{
    if (swtimedep_pbot == 0)
        return;

    int index0, index1;
    double fac0, fac1;

    model->timeloop->get_interpolation_factors(index0, index1, fac0, fac1, timedeptime);

    pbot = fac0 * timedeppbot[index0] + fac1 * timedeppbot[index1];
}

#ifndef USECUDA
void Thermo_vapor::get_thermo_field(Field3d* fld, Field3d* tmp, const std::string name, bool cyclic)
{
    const int kcells = grid->kcells;

    // BvS: getThermoField() is called from subgrid-model, before thermo(), so re-calculate the hydrostatic pressure
    // Pass dummy as rhoref,thvref to prevent overwriting base state
    double* restrict tmp2 = fields->atmp["tmp2"]->data;
    if (swupdatebasestate)
        calc_base_state(pref, prefh, &tmp2[0*kcells], &tmp2[1*kcells], &tmp2[2*kcells], &tmp2[3*kcells], exnref, exnrefh,
                fields->sp[thvar]->datamean, fields->sp["qt"]->datamean);

    if (name == "b")
        calc_buoyancy(fld->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, thvref);
    else if (name == "N2")
        calc_N2(fld->data, fields->sp[thvar]->data, grid->dzi, thvref);
    else
        throw 1;

    if (cyclic)
        grid->boundary_cyclic(fld->data);
}
#endif

#ifndef USECUDA
void Thermo_vapor::get_buoyancy_surf(Field3d* bfield)
{
    calc_buoyancy_bot(bfield->data, bfield->databot,
                      fields->sp[thvar]->data, fields->sp[thvar]->databot,
                      fields->sp["qt"]->data, fields->sp["qt"]->databot,
                      thvref, thvrefh);
    calc_buoyancy_fluxbot(bfield->datafluxbot, fields->sp[thvar]->databot, fields->sp[thvar]->datafluxbot, fields->sp["qt"]->databot, fields->sp["qt"]->datafluxbot, thvrefh);
}
#endif

#ifndef USECUDA
void Thermo_vapor::get_buoyancy_fluxbot(Field3d *bfield)
{
    calc_buoyancy_fluxbot(bfield->datafluxbot, fields->sp[thvar]->databot, fields->sp[thvar]->datafluxbot, fields->sp["qt"]->databot, fields->sp["qt"]->datafluxbot, thvrefh);
}
#endif

void Thermo_vapor::get_prog_vars(std::vector<std::string> *list)
{
    list->push_back(thvar);
    list->push_back("qt");
}

double Thermo_vapor::get_buoyancy_diffusivity()
{
    // Use the diffusivity from the liquid water potential temperature
    return fields->sp[thvar]->visc;
}


/**
 * This function calculates the hydrostatic pressure at full and half levels,
 * with option to return base state profiles like reference density and temperature
 * Solves: dpi/dz=-g/thv with pi=cp*(p/p0)**(rd/cp)
 * @param pref Pointer to output hydrostatic pressure array (full level)
 * @param prefh Pointer to output hydrostatic pressure array (half level)
 * @param rho Pointer to output density array (full level)
 * @param rhoh Pointer to output density array (half level)
 * @param thv Pointer to output virtual potential temperature array (full level)
 * @param thvh Pointer to output virtual potential temperature array (half level)
 * @param ex Pointer to output exner array (full level)
 * @param exh Pointer to output exner array (half level)
 * @param thlmean Pointer to input liq. water potential temperature array (horizontal mean, full level)
 * @param qtmean Pointer to input tot. moisture mix. ratio  array (horizontal mean, full level)
 */
void  Thermo_vapor::calc_base_state(double* restrict pref,    double* restrict prefh,
                                    double* restrict rho,     double* restrict rhoh,
                                    double* restrict thv,     double* restrict thvh,
                                    double* restrict ex,      double* restrict exh,
                                    double* restrict thlmean, double* restrict qtmean)
{
    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    const double thlsurf = interp2(thlmean[kstart-1], thlmean[kstart]);
    const double qtsurf  = interp2(qtmean[kstart-1],  qtmean[kstart]);

    // Calculate the values at the surface (half level == kstart)
    prefh[kstart] = pbot;
    exh[kstart]   = exner(prefh[kstart]);

    thvh[kstart]  = virtual_temperature_no_ql(exh[kstart], thlsurf, qtsurf);
    rhoh[kstart]  = pbot / (Rd * exh[kstart] * thvh[kstart]);

    // Calculate the first full level pressure
    pref[kstart]  = prefh[kstart] * std::exp(-grav * grid->z[kstart] / (Rd * exh[kstart] * thvh[kstart]));

    for (int k=kstart+1; k<kend+1; ++k)
    {
        // 1. Calculate remaining values (thv and rho) at full-level[k-1]
        ex[k-1]  = exner(pref[k-1]);
        thv[k-1] = virtual_temperature_no_ql(ex[k-1], thlmean[k-1], qtmean[k-1]);
        rho[k-1] = pref[k-1] / (Rd * ex[k-1] * thv[k-1]);

        // 2. Calculate pressure at half-level[k]
        prefh[k] = prefh[k-1] * std::exp(-grav * grid->dz[k-1] / (Rd * ex[k-1] * thv[k-1]));
        exh[k]   = exner(prefh[k]);

        // 3. Use interpolated conserved quantities to calculate half-level[k] values
        const double thli = interp2(thlmean[k-1], thlmean[k]);
        const double qti  = interp2(qtmean [k-1], qtmean [k]);


        thvh[k]  = virtual_temperature_no_ql(exh[k], thli, qti);
        rhoh[k]  = prefh[k] / (Rd * exh[k] * thvh[k]);

        // 4. Calculate pressure at full-level[k]
        pref[k] = pref[k-1] * std::exp(-grav * grid->dzh[k] / (Rd * exh[k] * thvh[k]));
    }

    pref[kstart-1] = 2.*prefh[kstart] - pref[kstart];
}

void Thermo_vapor::calc_buoyancy_tend_2nd(double* restrict wt, double* restrict thl, double* restrict qt,
                                          double* restrict ph, double* restrict thlh, double* restrict qth,
                                          double* restrict thvrefh)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    
    for (int k=grid->kstart+1; k<grid->kend; k++)
    {
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ij  = i + j*jj;
                thlh[ij] = interp2(thl[ijk-kk], thl[ijk]);
                qth[ij]  = interp2(qt[ijk-kk], qt[ijk]);
                wt[ijk] += buoyancy_no_ql(thlh[ij], qth[ij], thvrefh[k]);
            }
    }
}


void Thermo_vapor::calc_buoyancy(double* restrict b, double* restrict thl, double* restrict qt,
                                 double* restrict p, double* restrict thvref)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=0; k<grid->kcells; k++)
    {
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ij  = i + j*jj;
                b[ijk] = buoyancy_no_ql(thl[ijk], qt[ijk], thvref[k]);
            }
    }

    grid->boundary_cyclic(b);
}

void Thermo_vapor::calc_N2(double* restrict N2, double* restrict thl, double* restrict dzi,
                           double* restrict thvref)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                N2[ijk] = grav/thvref[k]*0.5*(thl[ijk+kk] - thl[ijk-kk])*dzi[k];
            }
}

void Thermo_vapor::calc_buoyancy_bot(double* restrict b,      double* restrict bbot,
                                     double* restrict thl,    double* restrict thlbot,
                                     double* restrict qt,     double* restrict qtbot,
                                     double* restrict thvref, double* restrict thvrefh)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    // assume no liquid water at the lowest model level
    for (int j=0; j<grid->jcells; j++)
        #pragma ivdep
        for (int i=0; i<grid->icells; i++)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            bbot[ij ] = buoyancy_no_ql(thlbot[ij], qtbot[ij], thvrefh[kstart]);
            b   [ijk] = buoyancy_no_ql(thl[ijk], qt[ijk], thvref[kstart]);
        }
}

void Thermo_vapor::calc_buoyancy_fluxbot(double* restrict bfluxbot, double* restrict thlbot, double* restrict thlfluxbot, double* restrict qtbot, double* restrict qtfluxbot,
        double* restrict thvrefh)
{
    const int jj = grid->icells;
    const int kstart = grid->kstart;

    // assume no liquid water at the lowest model level
    for (int j=0; j<grid->jcells; j++)
        #pragma ivdep
        for (int i=0; i<grid->icells; i++)
        {
            const int ij = i + j*jj;
            bfluxbot[ij] = buoyancy_flux_no_ql(thlbot[ij], thlfluxbot[ij], qtbot[ij], qtfluxbot[ij], thvrefh[kstart]);
        }
}

void Thermo_vapor::init_stat()
{
    // Add variables to the statistics
    if (stats->get_switch() == "1")
    {
        /* Add fixed base-state density and temperature profiles. Density should probably be in fields (?), but
           there the statistics are initialized before thermo->create() is called */
        stats->add_fixed_prof("rhoref",  "Full level basic state density", "kg m-3", "z",  fields->rhoref );
        stats->add_fixed_prof("rhorefh", "Half level basic state density", "kg m-3", "zh", fields->rhorefh);
        stats->add_fixed_prof("thvref",  "Full level basic state virtual potential temperature", "K", "z", thvref );
        stats->add_fixed_prof("thvrefh", "Half level basic state virtual potential temperature", "K", "zh",thvrefh);

        if (swupdatebasestate == 1)
        {
            stats->add_prof("ph",   "Full level hydrostatic pressure", "Pa",     "z" );
            stats->add_prof("phh",  "Half level hydrostatic pressure", "Pa",     "zh");
            stats->add_prof("rho",  "Full level density",  "kg m-3", "z" );
            stats->add_prof("rhoh", "Half level density",  "kg m-3", "zh");
        }
        else
        {
            stats->add_fixed_prof("ph",  "Full level hydrostatic pressure", "Pa", "z",  pref );
            stats->add_fixed_prof("phh", "Half level hydrostatic pressure", "Pa", "zh", prefh);
        }

        stats->add_prof("b", "Buoyancy", "m s-2", "z");
        for (int n=2; n<5; ++n)
        {
            std::stringstream ss;
            ss << n;
            std::string sn = ss.str();
            stats->add_prof("b"+sn, "Moment " +sn+" of the buoyancy", "(m s-2)"+sn,"z");
        }

        stats->add_prof("bgrad", "Gradient of the buoyancy", "m s-3", "zh");
        stats->add_prof("bw"   , "Turbulent flux of the buoyancy", "m2 s-3", "zh");
        stats->add_prof("bdiff", "Diffusive flux of the buoyancy", "m2 s-3", "zh");
        stats->add_prof("bflux", "Total flux of the buoyancy", "m2 s-3", "zh");
    }
}


void Thermo_vapor::init_column()
{
    // Add variables to the statistics
    if (model->column->get_switch() == "1")
    {

        model->column->add_prof("b", "Buoyancy", "m s-2", "z");

    }
}

void Thermo_vapor::init_cross()
{
    if (model->cross->get_switch() == "1")
    {
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

    }
}

void Thermo_vapor::init_dump()
{
    if (model->dump->get_switch() == "1")
    {
        // Get global cross-list from cross.cxx
        std::vector<std::string> *dumplist_global = model->dump->get_dumplist();

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

void Thermo_vapor::calc_buoyancy_tend_4th(double* restrict wt, double* restrict thl,  double* restrict qt,
                                          double* restrict ph, double* restrict thlh, double* restrict qth,
                                          double* restrict thvrefh)
{
    const int jj  = grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    for (int k=grid->kstart+1; k<grid->kend; k++)
    {
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk1;
                const int ij  = i + j*jj;
                thlh[ij] = interp4(thl[ijk-kk2], thl[ijk-kk1], thl[ijk], thl[ijk+kk1]);
                qth[ij]  = interp4(qt[ijk-kk2],  qt[ijk-kk1],  qt[ijk],  qt[ijk+kk1]);
                wt[ijk] += buoyancy_no_ql(thlh[ij], qth[ij], thvrefh[k]);
            }
    }
}

