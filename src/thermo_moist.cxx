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
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "thermo_moist.h"
#include "diff_smag2.h"
#include "defines.h"
#include "constants.h"
#include "fd.h"
#include "model.h"
#include "stats.h"
#include "master.h"
#include "cross.h"
#include "dump.h"

using fd::o2::interp2;
using fd::o4::interp4;
using namespace constants;

Thermo_moist::Thermo_moist(Model* modelin, Input* inputin) : Thermo(modelin, inputin)
{
    swthermo = "moist";

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

    nerror += inputin->get_item(&fields->sp[thvar]->visc, "fields", "svisc", thvar );
    nerror += inputin->get_item(&fields->sp["qt"]->visc, "fields", "svisc", "qt");
    nerror += inputin->get_item(&pbot, "thermo", "pbot", "");

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

    // BvS test for updating hydrostatic prssure during run
    // swupdate..=0 -> initial base state pressure used in saturation calculation
    // swupdate..=1 -> base state pressure updated before saturation calculation
    nerror += inputin->get_item(&swupdatebasestate, "thermo", "swupdatebasestate", ""); 

    // Remove the data from the input that is not used, to avoid warnings.
    if (master->mode == "init")
    {
        inputin->flag_as_used("thermo", "thvref0");
        inputin->flag_as_used("thermo", "pbot");
    }

    if (nerror)
        throw 1;
}

Thermo_moist::~Thermo_moist()
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

void Thermo_moist::init()
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

void Thermo_moist::create(Input* inputin)
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

    init_stat();
}

#ifndef USECUDA
void Thermo_moist::exec()
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
                               &fields->atmp["tmp2"]->data[2*kk], thvrefh);
    }
    else if (grid->swspatialorder == "4")
    {
        calc_buoyancy_tend_4th(fields->wt->data, fields->sp[thvar]->data, fields->sp["qt"]->data, prefh,
                               &fields->atmp["tmp2"]->data[0*kk], &fields->atmp["tmp2"]->data[1*kk],
                               &fields->atmp["tmp2"]->data[2*kk],
                               thvrefh);
    }
}
#endif

void Thermo_moist::get_mask(Field3d *mfield, Field3d *mfieldh, Mask *m)
{
    if (m->name == "ql")
    {
        calc_liquid_water(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
        calc_mask_ql(mfield->data, mfieldh->data, mfieldh->databot,
                     stats->nmask, stats->nmaskh, &stats->nmaskbot,
                     fields->atmp["tmp1"]->data);
    }
    else if (m->name == "qlcore")
    {
        calc_buoyancy(fields->atmp["tmp2"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, fields->atmp["tmp1"]->data,thvref);
        // calculate the mean buoyancy to determine positive buoyancy
        grid->calc_mean(fields->atmp["tmp2"]->datamean, fields->atmp["tmp2"]->data, grid->kcells);

        calc_liquid_water(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
        calc_mask_qlcore(mfield->data, mfieldh->data, mfieldh->databot,
                         stats->nmask, stats->nmaskh, &stats->nmaskbot,
                         fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp2"]->datamean);
    }
}

void Thermo_moist::calc_mask_ql(double* restrict mask, double* restrict maskh, double* restrict maskbot,
                                int* restrict nmask, int* restrict nmaskh, int* restrict nmaskbot,
                                double* restrict ql)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    for (int k=grid->kstart; k<grid->kend; k++)
    {
        nmask[k] = 0;
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ntmp = ql[ijk] > 0.;
                nmask[k] += ntmp;
                mask[ijk] = (double)ntmp;
            }
    }

    for (int k=grid->kstart; k<grid->kend+1; k++)
    {
        nmaskh[k] = 0;
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ntmp = (ql[ijk-kk] + ql[ijk]) > 0.;

                nmaskh[k] += ntmp;
                maskh[ijk] = (double)ntmp;
            }
    }

    // Set the mask for surface projected quantities
    // In this case: ql at surface
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            maskbot[ij] = maskh[ijk];
        }

    grid->boundary_cyclic(mask);
    grid->boundary_cyclic(maskh);
    grid->boundary_cyclic_2d(maskbot);

    master->sum(nmask , grid->kcells);
    master->sum(nmaskh, grid->kcells);
    *nmaskbot = nmaskh[grid->kstart];

    // BvS: should no longer be necessary now that the ql ghost cells are set to zero
    //nmaskh[kstart] = 0;
    //nmaskh[kend  ] = 0;
}

void Thermo_moist::calc_mask_qlcore(double* restrict mask, double* restrict maskh, double* restrict maskbot,
                                    int* restrict nmask, int* restrict nmaskh, int* restrict nmaskbot,
                                    double* restrict ql, double* restrict b, double* restrict bmean)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    for (int k=grid->kstart; k<grid->kend; k++)
    {
        nmask[k] = 0;
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ntmp = (ql[ijk] > 0.)*(b[ijk]-bmean[k] > 0.);
                nmask[k] += ntmp;
                mask[ijk] = (double)ntmp;
            }
    }

    for (int k=grid->kstart; k<grid->kend+1; k++)
    {
        nmaskh[k] = 0;
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ntmp = (ql[ijk-kk]+ql[ijk] > 0.)*(b[ijk-kk]+b[ijk]-bmean[k-1]-bmean[k] > 0.);
                nmaskh[k] += ntmp;
                maskh[ijk] = (double)ntmp;
            }
    }

    // Set the mask for surface projected quantities
    // In this case: qlcore at surface
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            maskbot[ij] = maskh[ijk];
        }

    grid->boundary_cyclic(mask);
    grid->boundary_cyclic(maskh);
    grid->boundary_cyclic_2d(maskbot);

    master->sum(nmask , grid->kcells);
    master->sum(nmaskh, grid->kcells);
    *nmaskbot = nmaskh[grid->kstart];
}

void Thermo_moist::exec_stats(Mask *m)
{
    const double NoOffset = 0.;

    // calc the buoyancy and its surface flux for the profiles
    calc_buoyancy(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, fields->atmp["tmp2"]->data, thvref);
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
        if (model->diff->get_name() == "smag2")
        {
            Diff_smag_2 *diffptr = static_cast<Diff_smag_2 *>(model->diff);
            stats->calc_diff_2nd(fields->atmp["tmp1"]->data, fields->w->data, fields->sd["evisc"]->data,
                                 m->profs["bdiff"].data, grid->dzhi,
                                 fields->atmp["tmp1"]->datafluxbot, fields->atmp["tmp1"]->datafluxtop, diffptr->tPr, sloc,
                                 fields->atmp["tmp4"]->data, stats->nmaskh);
        }
        else
        {
            stats->calc_diff_2nd(fields->atmp["tmp1"]->data, m->profs["bdiff"].data, grid->dzhi, fields->sp["th"]->visc, sloc,
                                 fields->atmp["tmp4"]->data, stats->nmaskh);
        }
    }
    else if (grid->swspatialorder == "4")
    {
        // take the diffusivity of temperature for that of buoyancy
        stats->calc_diff_4th(fields->atmp["tmp1"]->data, m->profs["bdiff"].data, grid->dzhi4, fields->sp["th"]->visc, sloc,
                             fields->atmp["tmp4"]->data, stats->nmaskh);
    }

    // calculate the total fluxes
    stats->add_fluxes(m->profs["bflux"].data, m->profs["bw"].data, m->profs["bdiff"].data);

    // calculate the liquid water stats
    calc_liquid_water(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
    stats->calc_mean(m->profs["ql"].data, fields->atmp["tmp1"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);
    stats->calc_count(fields->atmp["tmp1"]->data, m->profs["cfrac"].data, 0.,
                      fields->atmp["tmp3"]->data, stats->nmask);

    stats->calc_cover(fields->atmp["tmp1"]->data, fields->atmp["tmp4"]->databot, &stats->nmaskbot, &m->tseries["ccover"].data, 0.);
    stats->calc_path (fields->atmp["tmp1"]->data, fields->atmp["tmp4"]->databot, &stats->nmaskbot, &m->tseries["lwp"].data);

    // Calculate base state in tmp array
    if (swupdatebasestate == 1)
    {
        const int kcells = grid->kcells;
        double* restrict tmp2 = fields->atmp["tmp2"]->data;
        calc_base_state(&tmp2[0*kcells], &tmp2[1*kcells], &tmp2[2*kcells], &tmp2[3*kcells], 
                        &tmp2[4*kcells], &tmp2[5*kcells], &tmp2[6*kcells], &tmp2[7*kcells], 
                        fields->sp[thvar]->datamean, fields->sp["qt"]->datamean);

        for (int k=0; k<kcells; ++k)
        {
            m->profs["ph"  ].data[k] = tmp2[0*kcells+k];
            m->profs["phh" ].data[k] = tmp2[1*kcells+k];
            m->profs["rho" ].data[k] = tmp2[2*kcells+k];
            m->profs["rhoh"].data[k] = tmp2[3*kcells+k];
        } 
    }
}

void Thermo_moist::exec_cross()
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
            calc_buoyancy(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, fields->atmp["tmp2"]->data, thvref);
            nerror += cross->cross_simple(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, *it);
        }
        else if (*it == "ql")
        {
            calc_liquid_water(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
            nerror += cross->cross_simple(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, *it);
        }
        else if (*it == "blngrad")
        {
            calc_buoyancy(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, fields->atmp["tmp2"]->data, thvref);
            // Note: tmp1 twice used as argument -> overwritten in crosspath()
            nerror += cross->cross_lngrad(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, grid->dzi4, *it);
        }
        else if (*it == "qlpath")
        {
            calc_liquid_water(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
            // Note: tmp1 twice used as argument -> overwritten in crosspath()
            nerror += cross->cross_path(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, "qlpath");
        }
        else if (*it == "bbot" or *it == "bfluxbot")
        {
            calc_buoyancy_bot(fields->atmp["tmp1"]->data, fields->atmp["tmp1"]->databot, fields->sp[thvar ]->data, fields->sp[thvar]->databot, fields->sp["qt"]->data, fields->sp["qt"]->databot, thvref, thvrefh);
            calc_buoyancy_fluxbot(fields->atmp["tmp1"]->datafluxbot, fields->sp[thvar]->databot, fields->sp[thvar]->datafluxbot, fields->sp["qt"]->databot, fields->sp["qt"]->datafluxbot, thvrefh);

            if (*it == "bbot")
                nerror += cross->cross_plane(fields->atmp["tmp1"]->databot, fields->atmp["tmp1"]->data, "bbot");
            else if (*it == "bfluxbot")
                nerror += cross->cross_plane(fields->atmp["tmp1"]->datafluxbot, fields->atmp["tmp1"]->data, "bfluxbot");
        }
    }

    if (nerror)
        throw 1;
}

void Thermo_moist::exec_dump()
{
    for (std::vector<std::string>::const_iterator it=dumplist.begin(); it<dumplist.end(); ++it)
    {
        // TODO BvS restore getThermoField(), the combination of checkThermoField with getThermoField is more elegant... 
        if (*it == "b")
            calc_buoyancy(fields->atmp["tmp2"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, fields->atmp["tmp1"]->data, thvref);
        else if (*it == "ql")
            calc_liquid_water(fields->atmp["tmp2"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
        else
            throw 1;

        model->dump->save_dump(fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, *it);
    }
}

bool Thermo_moist::check_field_exists(const std::string name)
{
    if (name == "b" || name == "ql")
        return true;
    else
        return false;
}

#ifndef USECUDA
void Thermo_moist::get_thermo_field(Field3d* fld, Field3d* tmp, const std::string name)
{
    const int kcells = grid->kcells;

    // BvS: getThermoField() is called from subgrid-model, before thermo(), so re-calculate the hydrostatic pressure
    // Pass dummy as rhoref,thvref to prevent overwriting base state 
    double* restrict tmp2 = fields->atmp["tmp2"]->data;
    if (swupdatebasestate)
        calc_base_state(pref, prefh, &tmp2[0*kcells], &tmp2[1*kcells], &tmp2[2*kcells], &tmp2[3*kcells], exnref, exnrefh, 
                fields->sp[thvar]->datamean, fields->sp["qt"]->datamean);

    if (name == "b")
        calc_buoyancy(fld->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, tmp->data, thvref);
    else if (name == "ql")
        calc_liquid_water(fld->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
    else if (name == "N2")
        calc_N2(fld->data, fields->sp[thvar]->data, grid->dzi, thvref);
    else
        throw 1;
}
#endif

#ifndef USECUDA
void Thermo_moist::get_buoyancy_surf(Field3d* bfield)
{
    calc_buoyancy_bot(bfield->data, bfield->databot,
                      fields->sp[thvar]->data, fields->sp[thvar]->databot,
                      fields->sp["qt"]->data, fields->sp["qt"]->databot,
                      thvref, thvrefh);
    calc_buoyancy_fluxbot(bfield->datafluxbot, fields->sp[thvar]->databot, fields->sp[thvar]->datafluxbot, fields->sp["qt"]->databot, fields->sp["qt"]->datafluxbot, thvrefh);
}
#endif

#ifndef USECUDA
void Thermo_moist::get_buoyancy_fluxbot(Field3d *bfield)
{
    calc_buoyancy_fluxbot(bfield->datafluxbot, fields->sp[thvar]->databot, fields->sp[thvar]->datafluxbot, fields->sp["qt"]->databot, fields->sp["qt"]->datafluxbot, thvrefh);
}
#endif

void Thermo_moist::get_prog_vars(std::vector<std::string> *list)
{
    list->push_back(thvar);
    list->push_back("qt");
}

namespace
{
    // INLINE FUNCTIONS
    inline double buoyancy(const double exn, const double thl, const double qt, const double ql, const double thvref)
    {
        return grav * ((thl + Lv*ql/(cp*exn)) * (1. - (1. - Rv/Rd)*qt - Rv/Rd*ql) - thvref) / thvref;
    }

    inline double buoyancy_no_ql(const double thl, const double qt, const double thvref)
    {
        return grav * (thl * (1. - (1. - Rv/Rd)*qt) - thvref) / thvref;
    }

    inline double buoyancy_flux_no_ql(const double thl, const double thlflux, const double qt, const double qtflux, const double thvref)
    {
        return grav/thvref * (thlflux * (1. - (1.-Rv/Rd)*qt) - (1.-Rv/Rd)*thl*qtflux);
    }

    inline double esat(const double T)
    {
        const double x=std::max(-80.,T-T0);
        return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))));
    }

    inline double qsat(const double p, const double T)
    {
        return ep*esat(T)/(p-(1-ep)*esat(T));
    }

    inline double sat_adjust(const double thl, const double qt, const double p, const double exn)
    {
        int niter = 0, nitermax = 30;
        double ql, tl, tnr_old = 1.e9, tnr, qs=0;
        tl = thl * exn;
        tnr = tl;
        while (std::fabs(tnr-tnr_old)/tnr_old> 1e-5 && niter < nitermax)
        {
            ++niter;
            tnr_old = tnr;
            qs = qsat(p,tnr);
            tnr = tnr - (tnr+(Lv/cp)*qs-tl-(Lv/cp)*qt)/(1+(std::pow(Lv,2)*qs)/ (Rv*cp*std::pow(tnr,2)));
        }

        if (niter == nitermax)
        {  
            printf("Saturation adjustment not converged!! [thl=%f K, qt=%f kg/kg, p=%f p]\n",thl,qt,p);
            throw 1;
        }  

        ql = std::max(0.,qt - qs);
        return ql;
    }

    inline double exner(const double p)
    {
        return pow((p/p0),(Rd/cp));
    }
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
void  Thermo_moist::calc_base_state(double* restrict pref,    double* restrict prefh,
                                    double* restrict rho,     double* restrict rhoh,
                                    double* restrict thv,     double* restrict thvh,
                                    double* restrict ex,      double* restrict exh,
                                    double* restrict thlmean, double* restrict qtmean)
{
    const double rdcp = Rd/cp;

    const int kstart = grid->kstart;
    const int kend = grid->kend;

    double thlsurf = constants::dhuge;
    double qtsurf = constants::dhuge;

    if (grid->swspatialorder == "2")
    {
        thlsurf  = interp2(thlmean[kstart-1], thlmean[kstart]);
        qtsurf   = interp2(qtmean[kstart-1],  qtmean[kstart]);
    }
    else if (grid->swspatialorder == "4")
    {
        thlsurf  = interp4(thlmean[kstart-2], thlmean[kstart-1], thlmean[kstart], thlmean[kstart+1]);
        qtsurf   = interp4(qtmean[kstart-2],  qtmean[kstart-1],  qtmean[kstart],  qtmean[kstart+1]);
    }

    double ql,thli=0,qti=0,qli;

    // Calculate surface (half=kstart) values
    exh[kstart]   = exner(pbot);
    ql            = sat_adjust(thlsurf,qtsurf,pbot,exh[kstart]); 
    thvh[kstart]  = (thlsurf + Lv*ql/(cp*exh[kstart])) * (1. - (1. - Rv/Rd)*qtsurf - Rv/Rd*ql);
    prefh[kstart] = pbot;
    rhoh[kstart]  = pbot / (Rd * exh[kstart] * thvh[kstart]);

    // First full grid level pressure
    pref[kstart] = pow((pow(pbot,rdcp) - grav * pow(p0,rdcp) * grid->z[kstart] / (cp * thvh[kstart])),(1./rdcp)); 

    for (int k=kstart+1; k<kend+1; k++)
    {
        // 1. Calculate values at full level below zh[k] 
        ex[k-1]  = exner(pref[k-1]);
        ql       = sat_adjust(thlmean[k-1],qtmean[k-1],pref[k-1],ex[k-1]); 
        thv[k-1] = (thlmean[k-1] + Lv*ql/(cp*ex[k-1])) * (1. - (1. - Rv/Rd)*qtmean[k-1] - Rv/Rd*ql); 
        rho[k-1] = pref[k-1] / (Rd * ex[k-1] * thv[k-1]);

        // 2. Calculate half level pressure at zh[k] using values at z[k-1]
        prefh[k] = pow((pow(prefh[k-1],rdcp) - grav * pow(p0,rdcp) * grid->dz[k-1] / (cp * thv[k-1])),(1./rdcp));

        // 3. Interpolate conserved variables to zh[k] and calculate virtual temp and ql
        if (grid->swspatialorder == "2")
        {
            thli   = interp2(thlmean[k-1],thlmean[k]);
            qti    = interp2(qtmean[k-1],qtmean[k]);
        }
        else if (grid->swspatialorder == "4")
        {
            thli   = interp4(thlmean[k-2],thlmean[k-1],thlmean[k],thlmean[k+1]);
            qti    = interp4(qtmean[k-2],qtmean[k-1],qtmean[k],qtmean[k+1]);
        }

        exh[k]   = exner(prefh[k]);
        qli      = sat_adjust(thli,qti,prefh[k],exh[k]);
        thvh[k]  = (thli + Lv*qli/(cp*exh[k])) * (1. - (1. - Rv/Rd)*qti - Rv/Rd*qli); 
        rhoh[k]  = prefh[k] / (Rd * exh[k] * thvh[k]); 

        // 4. Calculate full level pressure at z[k]
        pref[k]  = pow((pow(pref[k-1],rdcp) - grav * pow(p0,rdcp) * grid->dzh[k] / (cp * thvh[k])),(1./rdcp)); 
    }

    // Fill bottom and top full level ghost cells 
    if (grid->swspatialorder == "2")
    {
        pref[kstart-1] = 2.*prefh[kstart] - pref[kstart];
        pref[kend]     = 2.*prefh[kend]   - pref[kend-1];
    }
    else if (grid->swspatialorder == "4")
    {
        pref[kstart-1] = (8./3.)*prefh[kstart] - 2.*pref[kstart] + (1./3.)*pref[kstart+1];
        pref[kstart-2] = 8.*prefh[kstart]      - 9.*pref[kstart] + 2.*pref[kstart+1];
        pref[kend]     = (8./3.)*prefh[kend]   - 2.*pref[kend-1] + (1./3.)*pref[kend-2];
        pref[kend+1]   = 8.*prefh[kend]        - 9.*pref[kend-1] + 2.*pref[kend-2];
    }
}

void Thermo_moist::calc_buoyancy_tend_2nd(double* restrict wt, double* restrict thl, double* restrict qt,
                                          double* restrict ph, double* restrict thlh, double* restrict qth,
                                          double* restrict ql, double* restrict thvrefh)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    double tl, exnh;

    for (int k=grid->kstart+1; k<grid->kend; k++)
    {
        exnh = exner(ph[k]);
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ij  = i + j*jj;
                thlh[ij] = interp2(thl[ijk-kk], thl[ijk]);
                qth[ij]  = interp2(qt[ijk-kk], qt[ijk]);
                tl       = thlh[ij] * exnh;
                // Calculate first estimate of ql using Tl
                // if ql(Tl)>0, saturation adjustment routine needed
                ql[ij]  = qth[ij]-qsat(ph[k],tl);
            }
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ij  = i + j*jj;
                if (ql[ij]>0)   // already doesn't vectorize because of iteration in sat_adjust()
                {
                    ql[ij] = sat_adjust(thlh[ij], qth[ij], ph[k], exnh);
                }
                else
                    ql[ij] = 0.;
            }
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ij  = i + j*jj;
                wt[ijk] += buoyancy(exnh, thlh[ij], qth[ij], ql[ij], thvrefh[k]);
            }
    }
}

void Thermo_moist::calc_buoyancy_tend_4th(double* restrict wt, double* restrict thl,  double* restrict qt, 
                                          double* restrict ph, double* restrict thlh, double* restrict qth,
                                          double* restrict ql, double* restrict thvrefh)
{
    const int jj  = grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    double tl, exnh;

    for (int k=grid->kstart+1; k<grid->kend; k++)
    {
        exnh = exner(ph[k]);
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk1;
                const int ij  = i + j*jj;
                thlh[ij] = interp4(thl[ijk-kk2], thl[ijk-kk1], thl[ijk], thl[ijk+kk1]);
                qth[ij]  = interp4(qt[ijk-kk2],  qt[ijk-kk1],  qt[ijk],  qt[ijk+kk1]);
                tl       = thlh[ij] * exnh;
                // Calculate first estimate of ql using Tl
                // if ql(Tl)>0, saturation adjustment routine needed
                ql[ij]  = qth[ij]-qsat(ph[k],tl);   
            }
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ij = i + j*jj;
                if (ql[ij]>0)   // already doesn't vectorize because of iteration in sat_adjust()
                    ql[ij] = sat_adjust(thlh[ij], qth[ij], ph[k], exnh);
                else
                    ql[ij] = 0.;
            }
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk1;
                const int ij  = i + j*jj;
                wt[ijk] += buoyancy(exnh, thlh[ij], qth[ij], ql[ij], thvrefh[k]);
            }
    }
}

void Thermo_moist::calc_buoyancy(double* restrict b, double* restrict thl, double* restrict qt,
                                 double* restrict p, double* restrict ql, double* restrict thvref)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    double tl, ex;

    for (int k=0; k<grid->kcells; k++)
    {
        ex = exner(p[k]);
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ij  = i + j*jj;
                tl  = thl[ijk] * ex;
                ql[ij]  = qt[ijk]-qsat(p[k],tl);   // not real ql, just estimate
            }

        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ij  = i + j*jj;
                if (ql[ij] > 0)
                    ql[ij] = sat_adjust(thl[ijk], qt[ijk], p[k], ex);
                else
                    ql[ij] = 0.;
            }

        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ij  = i + j*jj;
                b[ijk] = buoyancy(ex, thl[ijk], qt[ijk], ql[ij], thvref[k]);
            }
    }
}

void Thermo_moist::calc_liquid_water(double* restrict ql, double* restrict thl, double* restrict qt, double* restrict p)
{
    double ex;

    const int jj = grid->icells;
    const int kk = grid->ijcells;

    // Fill ghost cells with zeros to prevent problems in calculating ql or qlcore masks 
    for (int k=0; k<grid->kstart; k++)
    {
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                ql[ijk] = 0.;
            }
    }

    for (int k=grid->kend; k<grid->kcells; k++)
    {
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                ql[ijk] = 0.;
            }
    }

    // Calculate the ql field
    for (int k=grid->kstart; k<grid->kend; k++)
    {
        ex = exner(p[k]);
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                ql[ijk] = sat_adjust(thl[ijk], qt[ijk], p[k], ex);
            }
    }
}

void Thermo_moist::calc_N2(double* restrict N2, double* restrict thl, double* restrict dzi,
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

void Thermo_moist::calc_buoyancy_bot(double* restrict b,      double* restrict bbot,
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

void Thermo_moist::calc_buoyancy_fluxbot(double* restrict bfluxbot, double* restrict thlbot, double* restrict thlfluxbot, double* restrict qtbot, double* restrict qtfluxbot,
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

void Thermo_moist::init_stat()
{
    // Add variables to the statistics
    if (stats->getSwitch() == "1")
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

        stats->add_prof("ql", "Liquid water mixing ratio", "kg kg-1", "z");
        stats->add_prof("cfrac", "Cloud fraction", "-","z");

        stats->add_time_series("lwp", "Liquid water path", "kg m-2");
        stats->add_time_series("ccover", "Projected cloud cover", "-");
    }
}

void Thermo_moist::init_cross()
{
    if (model->cross->get_switch() == "1")
    {
        allowedcrossvars.push_back("b");
        allowedcrossvars.push_back("bbot");
        allowedcrossvars.push_back("bfluxbot");
        if (grid->swspatialorder == "4")
            allowedcrossvars.push_back("blngrad");
        allowedcrossvars.push_back("ql");
        allowedcrossvars.push_back("qlpath");

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

void Thermo_moist::init_dump()
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
