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

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <sstream>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "finite_difference.h"
#include "model.h"
#include "stats.h"
#include "cross.h"
#include "dump.h"
#include "diff_smag2.h"

Fields::Fields(Model *modelin, Input *inputin)
{
    model  = modelin;
    grid   = model->grid;
    master = model->master;

    calc_mean_profs = false;

    // Initialize the pointers.
    rhoref  = 0;
    rhorefh = 0;
    umodel  = 0;
    vmodel  = 0;

    // Initialize GPU pointers
    rhoref_g  = 0;
    rhorefh_g = 0;

    // input parameters
    int nerror = 0;

    // obligatory parameters
    nerror += inputin->get_item(&visc, "fields", "visc", "");

    // read the name of the passive scalars
    std::vector<std::string> slist;
    nerror += inputin->get_list(&slist, "fields", "slist", "");

    // initialize the scalars
    for (std::vector<std::string>::const_iterator it=slist.begin(); it!=slist.end(); ++it)
    {
        init_prognostic_field(*it, *it, "-");
        nerror += inputin->get_item(&sp[*it]->visc, "fields", "svisc", *it);
    }

    if (nerror)
        throw 1;

    // initialize the basic set of fields
    init_momentum_field(u, ut, "u", "U velocity", "m s-1");
    init_momentum_field(v, vt, "v", "V velocity", "m s-1");
    init_momentum_field(w, wt, "w", "Vertical velocity", "m s-1");
    init_diagnostic_field("p", "Pressure", "Pa");

    // Set a default of 4 temporary fields. Other classes can increase this number
    // before the init phase, where they are initialized in Fields::init()
    n_tmp_fields = 4;

    // Remove the data from the input that is not used in run mode, to avoid warnings.
    if (master->mode == "run")
    {
        inputin->flag_as_used("fields", "rndamp");
        inputin->flag_as_used("fields", "rndexp");
        inputin->flag_as_used("fields", "rndseed");
        inputin->flag_as_used("fields", "rndz");

        inputin->flag_as_used("fields", "vortexnpair");
        inputin->flag_as_used("fields", "vortexamp"  );
        inputin->flag_as_used("fields", "vortexaxis" );
    }
}

Fields::~Fields()
{
    // DEALLOCATE ALL THE FIELDS
    // deallocate the prognostic velocity fields
    for (FieldMap::iterator it=mp.begin(); it!=mp.end(); ++it)
        delete it->second;

    // deallocate the velocity tendency fields
    for (FieldMap::iterator it=mt.begin(); it!=mt.end(); ++it)
        delete it->second;

    // deallocate the prognostic scalar fields
    for (FieldMap::iterator it=sp.begin(); it!=sp.end(); ++it)
        delete it->second;

    // deallocate the scalar tendency fields
    for (FieldMap::iterator it=st.begin(); it!=st.end(); ++it)
        delete it->second;

    // deallocate the diagnostic scalars
    for (FieldMap::iterator it=sd.begin(); it!=sd.end(); ++it)
        delete it->second;

    // deallocate the tmp fields
    for (FieldMap::iterator it=atmp.begin(); it!=atmp.end(); ++it)
        delete it->second;

    // delete the arrays
    delete[] rhoref;
    delete[] rhorefh;
    delete[] umodel;
    delete[] vmodel;

#ifdef USECUDA
    clear_device();
#endif
}

void Fields::init()
{
    // set the convenience pointers
    stats = model->stats;

    int nerror = 0;

    // ALLOCATE ALL THE FIELDS
    // allocate the prognostic velocity fields
    for (FieldMap::iterator it=mp.begin(); it!=mp.end(); ++it)
        nerror += it->second->init();

    // allocate the velocity tendency fields
    for (FieldMap::iterator it=mt.begin(); it!=mt.end(); ++it)
        nerror += it->second->init();

    // allocate the prognostic scalar fields
    for (FieldMap::iterator it=sp.begin(); it!=sp.end(); ++it)
        nerror += it->second->init();

    // allocate the scalar tendency fields
    for (FieldMap::iterator it=st.begin(); it!=st.end(); ++it)
        nerror += it->second->init();

    // allocate the diagnostic scalars
    for (FieldMap::iterator it=sd.begin(); it!=sd.end(); ++it)
        nerror += it->second->init();

    // now that all classes have been able to set the minimum number of tmp fields, initialize them
    for (int i=1; i<=n_tmp_fields; ++i)
    {
        // BvS: the cast to long long is unfortunately necessary for Intel compilers
        // which don't seem to have the full c++11 implementation
        std::string name = "tmp" + std::to_string(static_cast<long long>(i));
        init_tmp_field(name, "", "");
    }

    // allocate the tmp fields
    for (FieldMap::iterator it=atmp.begin(); it!=atmp.end(); ++it)
        nerror += it->second->init();

    if (nerror > 0)
        throw 1;

    // allocate the base density profiles
    rhoref  = new double[grid->kcells];
    rhorefh = new double[grid->kcells];

    // \TODO Define a reference density. Needs to be replaced once anelastic is there
    // BvS: Always init rhoref at 1 for situation with e.g. thermo=0? For anelastic, overwrite it.
    for (int k=0; k<grid->kcells; ++k)
    {
        rhoref[k] = 1.;
        rhorefh[k] = 1.; 
    }

    // allocate help arrays for statistics;
    umodel = new double[grid->kcells];
    vmodel = new double[grid->kcells];

    // Initialize at zero
    for (int k=0; k<grid->kcells; ++k)
    {
        umodel[k] = 0.;
        vmodel[k] = 0.; 
    }

    // Get global cross-list from cross.cxx
    std::vector<std::string> *crosslist_global = model->cross->get_crosslist(); 

    // Check different type of crosses and put them in their respective lists 
    for (FieldMap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
    {
        check_added_cross(it->first, "",        crosslist_global, &crosssimple);
        check_added_cross(it->first, "lngrad",  crosslist_global, &crosslngrad);
        check_added_cross(it->first, "bot",     crosslist_global, &crossbot);
        check_added_cross(it->first, "top",     crosslist_global, &crosstop);
        check_added_cross(it->first, "fluxbot", crosslist_global, &crossfluxbot);
        check_added_cross(it->first, "fluxtop", crosslist_global, &crossfluxtop);
    }

    for (FieldMap::const_iterator it=sd.begin(); it!=sd.end(); ++it)
    {
        check_added_cross(it->first, "",        crosslist_global, &crosssimple);
        check_added_cross(it->first, "lngrad",  crosslist_global, &crosslngrad);
    }

    // Get global dump-list from cross.cxx
    std::vector<std::string> *dumplist_global = model->dump->get_dumplist(); 

    // Check if fields in dumplist are diagnostic fields, if not delete them and print warning
    std::vector<std::string>::iterator dumpvar=dumplist_global->begin();
    while (dumpvar != dumplist_global->end())
    {
        if (sd.count(*dumpvar))
        {
            // Remove variable from global list, put in local list
            dumplist.push_back(*dumpvar);
            dumplist_global->erase(dumpvar); // erase() returns iterator of next element..
        }
        else
            ++dumpvar;
    }
}

void Fields::check_added_cross(std::string var, std::string type, std::vector<std::string> *crosslist, std::vector<std::string> *typelist)
{
    std::vector<std::string>::iterator position;

    position = std::find(crosslist->begin(), crosslist->end(), var + type);
    if (position != crosslist->end()) 
    {
        // don't allow lngrad in 2nd order mode
        if (!(type == "lngrad" && grid->swspatialorder == "2"))
        {
            typelist->push_back(var);
            crosslist->erase(position);
        }
    }
}

#ifndef USECUDA
void Fields::exec()
{
    // calculate the means for the prognostic scalars
    if (calc_mean_profs)
    {
        for (FieldMap::iterator it=sp.begin(); it!=sp.end(); ++it)
            grid->calc_mean(it->second->datamean, it->second->data, grid->kcells);
    }
}
#endif

void Fields::get_mask(Field3d *mfield, Field3d *mfieldh, Mask *m)
{
    if (m->name == "wplus")
        calc_mask_wplus(mfield->data, mfieldh->data, mfieldh->databot, 
                        stats->nmask, stats->nmaskh, &stats->nmaskbot, w->data);
    else if (m->name == "wmin")                                                  
        calc_mask_wmin(mfield->data, mfieldh->data, mfieldh->databot,
                       stats->nmask, stats->nmaskh, &stats->nmaskbot, w->data);
}

void Fields::calc_mask_wplus(double* restrict mask, double* restrict maskh, double* restrict maskbot,
        int* restrict nmask, int* restrict nmaskh, int* restrict nmaskbot,
        double* restrict w)
{
    int ijk,ij,jj,kk,kstart;

    jj = grid->icells;
    kk = grid->ijcells;
    kstart = grid->kstart;

    int ntmp;

    for (int k=grid->kstart; k<grid->kend; k++)
    {
        nmask[k] = 0;
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                ijk = i + j*jj + k*kk;
                ntmp = (w[ijk] + w[ijk+kk]) > 0.;
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
                ijk = i + j*jj + k*kk;
                ntmp = w[ijk] > 0.;
                nmaskh[k] += ntmp;
                maskh[ijk] = (double)ntmp;
            }
    }

    // Set the mask for surface projected quantities
    // In this case: velocity at surface, so zero
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            ij  = i + j*jj;
            ijk = i + j*jj + kstart*kk;
            maskbot[ij] = maskh[ijk];
        }

    grid->boundary_cyclic(mask);
    grid->boundary_cyclic(maskh);
    grid->boundary_cyclic_2d(maskbot);

    master->sum(nmask , grid->kcells);
    master->sum(nmaskh, grid->kcells);
    *nmaskbot = nmaskh[grid->kstart];
}

void Fields::calc_mask_wmin(double* restrict mask, double* restrict maskh, double* restrict maskbot,
        int* restrict nmask, int* restrict nmaskh, int* restrict nmaskbot,
        double* restrict w)
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
                const int ntmp = (w[ijk] + w[ijk+kk]) <= 0.;
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
                const int ntmp = w[ijk] <= 0.;
                nmaskh[k] += ntmp;
                maskh[ijk] = (double)ntmp;
            }
    }

    // Set the mask for surface projected quantities
    // In this case: velocity at surface, so zero
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

void Fields::exec_stats(Mask *m)
{
    // define locations
    const int uloc[] = {1,0,0};
    const int vloc[] = {0,1,0};
    const int wloc[] = {0,0,1};
    const int sloc[] = {0,0,0};

    const int uwloc[] = {1,0,1};
    const int vwloc[] = {0,1,1};

    const double NoOffset = 0.;

    // save the area coverage of the mask
    stats->calc_area(m->profs["area" ].data, sloc, stats->nmask );
    stats->calc_area(m->profs["areah"].data, wloc, stats->nmaskh);

    // start with the stats on the w location, to make the wmean known for the flux calculations
    stats->calc_mean(m->profs["w"].data, w->data, NoOffset, wloc, atmp["tmp4"]->data, stats->nmaskh);
    for (int n=2; n<5; ++n)
    {
        std::stringstream ss;
        ss << n;
        std::string sn = ss.str();
        stats->calc_moment(w->data, m->profs["w"].data, m->profs["w"+sn].data, n, wloc,
                          atmp["tmp4"]->data, stats->nmaskh);
    }

    // calculate the stats on the u location
    // interpolate the mask horizontally onto the u coordinate
    grid->interpolate_2nd(atmp["tmp1"]->data, atmp["tmp3"]->data, sloc, uloc);
    stats->calc_mean(m->profs["u"].data, u->data, grid->utrans, uloc, atmp["tmp1"]->data, stats->nmask);
    stats->calc_mean(umodel            , u->data, NoOffset   , uloc, atmp["tmp1"]->data, stats->nmask);
    for (int n=2; n<5; ++n)
    {
        std::stringstream ss;
        ss << n;
        std::string sn = ss.str();
        stats->calc_moment(u->data, umodel, m->profs["u"+sn].data, n, uloc,
                          atmp["tmp1"]->data, stats->nmask);
    }

    // interpolate the mask on half level horizontally onto the u coordinate
    grid->interpolate_2nd(atmp["tmp1"]->data, atmp["tmp4"]->data, wloc, uwloc);
    if (grid->swspatialorder == "2")
    {
        stats->calc_grad_2nd(u->data, m->profs["ugrad"].data, grid->dzhi, uloc,
                            atmp["tmp1"]->data, stats->nmaskh);
        stats->calc_flux_2nd(u->data, umodel, w->data, m->profs["w"].data,
                            m->profs["uw"].data, atmp["tmp2"]->data, uloc,
                            atmp["tmp1"]->data, stats->nmaskh);
        if (model->diff->get_switch() == "smag2")
            stats->calc_diff_2nd(u->data, w->data, sd["evisc"]->data,
                                m->profs["udiff"].data, grid->dzhi,
                                u->datafluxbot, u->datafluxtop, 1., uloc,
                                atmp["tmp1"]->data, stats->nmaskh);
        else
            stats->calc_diff_2nd(u->data, m->profs["udiff"].data, grid->dzhi, visc, uloc,
                                atmp["tmp1"]->data, stats->nmaskh);

    }
    else if (grid->swspatialorder == "4")
    {
        stats->calc_grad_4th(u->data, m->profs["ugrad"].data, grid->dzhi4, uloc,
                            atmp["tmp1"]->data, stats->nmaskh);
        stats->calc_flux_4th(u->data, w->data, m->profs["uw"].data, atmp["tmp2"]->data, uloc,
                            atmp["tmp1"]->data, stats->nmaskh);
        stats->calc_diff_4th(u->data, m->profs["udiff"].data, grid->dzhi4, visc, uloc,
                            atmp["tmp1"]->data, stats->nmaskh);
    }

    // calculate the stats on the v location
    grid->interpolate_2nd(atmp["tmp1"]->data, atmp["tmp3"]->data, sloc, vloc);
    stats->calc_mean(m->profs["v"].data, v->data, grid->vtrans, vloc, atmp["tmp1"]->data, stats->nmask);
    stats->calc_mean(vmodel            , v->data, NoOffset   , vloc, atmp["tmp1"]->data, stats->nmask);
    for (int n=2; n<5; ++n)
    {
        std::stringstream ss;
        ss << n;
        std::string sn = ss.str();
        stats->calc_moment(v->data, vmodel, m->profs["v"+sn].data, n, vloc,
                          atmp["tmp1"]->data, stats->nmask);
    }

    // interpolate the mask on half level horizontally onto the u coordinate
    grid->interpolate_2nd(atmp["tmp1"]->data, atmp["tmp4"]->data, wloc, vwloc);
    if (grid->swspatialorder == "2")
    {
        stats->calc_grad_2nd(v->data, m->profs["vgrad"].data, grid->dzhi, vloc,
                            atmp["tmp1"]->data, stats->nmaskh);
        stats->calc_flux_2nd(v->data, vmodel, w->data, m->profs["w"].data,
                            m->profs["vw"].data, atmp["tmp2"]->data, vloc,
                            atmp["tmp1"]->data, stats->nmaskh);
        if (model->diff->get_switch() == "smag2")
            stats->calc_diff_2nd(v->data, w->data, sd["evisc"]->data,
                                m->profs["vdiff"].data, grid->dzhi,
                                v->datafluxbot, v->datafluxtop, 1., vloc,
                                atmp["tmp1"]->data, stats->nmaskh);
        else
            stats->calc_diff_2nd(v->data, m->profs["vdiff"].data, grid->dzhi, visc, vloc,
                                atmp["tmp1"]->data, stats->nmaskh);

    }
    else if (grid->swspatialorder == "4")
    {
        stats->calc_grad_4th(v->data, m->profs["vgrad"].data, grid->dzhi4, vloc,
                            atmp["tmp1"]->data, stats->nmaskh);
        stats->calc_flux_4th(v->data, w->data, m->profs["vw"].data, atmp["tmp2"]->data, vloc,
                            atmp["tmp1"]->data, stats->nmaskh);
        stats->calc_diff_4th(v->data, m->profs["vdiff"].data, grid->dzhi4, visc, vloc,
                            atmp["tmp1"]->data, stats->nmaskh);
    }

    // calculate stats for the prognostic scalars
    Diff_smag_2 *diffptr = static_cast<Diff_smag_2 *>(model->diff);
    for (FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
    {
        stats->calc_mean(m->profs[it->first].data, it->second->data, NoOffset, sloc, atmp["tmp3"]->data, stats->nmask);
        for (int n=2; n<5; ++n)
        {
            std::stringstream ss;
            ss << n;
            std::string sn = ss.str();
            stats->calc_moment(it->second->data, m->profs[it->first].data, m->profs[it->first+sn].data, n, sloc,
                    atmp["tmp3"]->data, stats->nmask);
        }
        if (grid->swspatialorder == "2")
        {
            stats->calc_grad_2nd(it->second->data, m->profs[it->first+"grad"].data, grid->dzhi, sloc,
                                atmp["tmp4"]->data, stats->nmaskh);
            stats->calc_flux_2nd(it->second->data, m->profs[it->first].data, w->data, m->profs["w"].data,
                                m->profs[it->first+"w"].data, atmp["tmp1"]->data, sloc,
                                atmp["tmp4"]->data, stats->nmaskh);
            if (model->diff->get_switch() == "smag2")
                stats->calc_diff_2nd(it->second->data, w->data, sd["evisc"]->data,
                                    m->profs[it->first+"diff"].data, grid->dzhi,
                                    it->second->datafluxbot, it->second->datafluxtop, diffptr->tPr, sloc,
                                    atmp["tmp4"]->data, stats->nmaskh);
            else
                stats->calc_diff_2nd(it->second->data, m->profs[it->first+"diff"].data, grid->dzhi, it->second->visc, sloc,
                                    atmp["tmp4"]->data, stats->nmaskh);
        }
        else if (grid->swspatialorder == "4")
        {
            stats->calc_grad_4th(it->second->data, m->profs[it->first+"grad"].data, grid->dzhi4, sloc,
                                atmp["tmp4"]->data, stats->nmaskh);
            stats->calc_flux_4th(it->second->data, w->data, m->profs[it->first+"w"].data, atmp["tmp1"]->data, sloc,
                                atmp["tmp4"]->data, stats->nmaskh);
            stats->calc_diff_4th(it->second->data, m->profs[it->first+"diff"].data, grid->dzhi4, it->second->visc, sloc,
                                atmp["tmp4"]->data, stats->nmaskh);
        }
    }

    // Calculate pressure statistics
    stats->calc_mean(m->profs["p"].data, sd["p"]->data, NoOffset, sloc, atmp["tmp3"]->data, stats->nmask);
    stats->calc_moment(sd["p"]->data, m->profs["p"].data, m->profs["p2"].data, 2, sloc,
                      atmp["tmp1"]->data, stats->nmask);
    if (grid->swspatialorder == "2")
    {
        stats->calc_grad_2nd(sd["p"]->data, m->profs["pgrad"].data, grid->dzhi, sloc,
                             atmp["tmp4"]->data, stats->nmaskh);
        stats->calc_flux_2nd(sd["p"]->data, m->profs["p"].data, w->data, m->profs["w"].data,
                            m->profs["pw"].data, atmp["tmp1"]->data, sloc,
                            atmp["tmp4"]->data, stats->nmaskh);
    }
    else if (grid->swspatialorder == "4")
    {
        stats->calc_grad_4th(sd["p"]->data, m->profs["pgrad"].data, grid->dzhi4, sloc,
                             atmp["tmp4"]->data, stats->nmaskh);
        stats->calc_flux_4th(sd["p"]->data, w->data, m->profs["pw"].data, atmp["tmp1"]->data, sloc,
                             atmp["tmp4"]->data, stats->nmaskh);
    }

    // calculate the total fluxes
    stats->add_fluxes(m->profs["uflux"].data, m->profs["uw"].data, m->profs["udiff"].data);
    stats->add_fluxes(m->profs["vflux"].data, m->profs["vw"].data, m->profs["vdiff"].data);
    for (FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
        stats->add_fluxes(m->profs[it->first+"flux"].data, m->profs[it->first+"w"].data, m->profs[it->first+"diff"].data);

    if (model->diff->get_switch() == "smag2")
        stats->calc_mean(m->profs["evisc"].data, sd["evisc"]->data, NoOffset, sloc, atmp["tmp3"]->data, stats->nmask);
}

void Fields::set_calc_mean_profs(bool sw)
{
    calc_mean_profs = sw;
}

void Fields::set_minimum_tmp_fields(int n)
{
    n_tmp_fields = std::max(n_tmp_fields, n);
}

void Fields::init_momentum_field(Field3d*& fld, Field3d*& fldt, std::string fldname, std::string longname, std::string unit)
{
    if (mp.find(fldname)!=mp.end())
    {
        master->print_error("\"%s\" already exists\n", fldname.c_str());
        throw 1;
    }

    // add a new prognostic momentum variable
    mp[fldname] = new Field3d(grid, master, fldname, longname, unit);

    // add a new tendency for momentum variable
    std::string fldtname  = fldname + "t";
    std::string tunit     = unit + "s-1";
    std::string tlongname = "Tendency of " + longname;
    mt[fldname] = new Field3d(grid, master, fldtname, tlongname, tunit);

    // TODO remove these from the model?
    fld  = mp[fldname];
    fldt = mt[fldname];

    // add the prognostic variable and its tendency to the collection
    // of all fields and tendencies
    a [fldname] = mp[fldname];
    ap[fldname] = mp[fldname];
    at[fldname] = mt[fldname];
}

void Fields::init_prognostic_field(std::string fldname, std::string longname, std::string unit)
{
    if (sp.find(fldname)!=sp.end())
    {
        master->print_error("\"%s\" already exists\n", fldname.c_str());
        throw 1;
    }

    // add a new scalar variable
    sp[fldname] = new Field3d(grid, master, fldname,longname, unit);

    // add a new tendency for scalar variable
    std::string fldtname  = fldname + "t";
    std::string tlongname = "Tendency of " + longname;
    std::string tunit     = unit + "s-1";
    st[fldname] = new Field3d(grid, master, fldtname,tlongname, tunit);

    // add the prognostic variable and its tendency to the collection
    // of all fields and tendencies
    a [fldname] = sp[fldname];
    ap[fldname] = sp[fldname];
    at[fldname] = st[fldname];
}

void Fields::init_diagnostic_field(std::string fldname,std::string longname, std::string unit)
{
    if (sd.find(fldname)!=sd.end())
    {
        master->print_error("\"%s\" already exists\n", fldname.c_str());
        throw 1;
    }

    sd[fldname] = new Field3d(grid, master, fldname, longname, unit);
    a [fldname] = sd[fldname];
}

void Fields::init_tmp_field(std::string fldname,std::string longname, std::string unit)
{
    if (atmp.find(fldname)!=atmp.end())
    {
        master->print_error("\"%s\" already exists\n", fldname.c_str());
        throw 1;
    }

    atmp[fldname] = new Field3d(grid, master, fldname, longname, unit);
}

void Fields::create(Input *inputin)
{
    int nerror = 0;

    // Randomize the momentum
    nerror += randomize(inputin, "u", u->data);
    nerror += randomize(inputin, "w", w->data);
    // Only add perturbation to v in case of a 3d run.
    if (grid->jtot > 1)
        nerror += randomize(inputin, "v", v->data);

    // Randomize the scalars
    for (FieldMap::iterator it=sp.begin(); it!=sp.end(); ++it)
        nerror += randomize(inputin, it->first, it->second->data);

    // Add Vortices
    nerror += add_vortex_pair(inputin);

    // Add the mean profiles to the fields
    nerror += add_mean_prof(inputin, "u", mp["u"]->data, grid->utrans);
    nerror += add_mean_prof(inputin, "v", mp["v"]->data, grid->vtrans);

    for (FieldMap::iterator it=sp.begin(); it!=sp.end(); ++it)
        nerror += add_mean_prof(inputin, it->first, it->second->data, 0.);

    // set w equal to zero at the boundaries, just to be sure
    int lbot = grid->kstart*grid->ijcells;
    int ltop = grid->kend  *grid->ijcells;
    for (int l=0; l<grid->ijcells; ++l)
    {
        w->data[lbot+l] = 0.;
        w->data[ltop+l] = 0.;
    }

    if (nerror)
        throw 1;
}

int Fields::randomize(Input* inputin, std::string fld, double* restrict data)
{
    int nerror = 0;

    // Set mpiid as random seed to avoid having the same field at all procs
    int static seed = 0;

    if (!seed)
    {
        nerror += inputin->get_item(&seed, "fields", "rndseed", "", 0);
        seed += master->mpiid + 2;
        std::srand(seed);
    }

    const int jj = grid->icells;
    const int kk = grid->ijcells;

    // Look up the specific randomizer variables.
    nerror += inputin->get_item(&rndamp, "fields", "rndamp", fld, 0.);
    nerror += inputin->get_item(&rndz  , "fields", "rndz"  , fld, 0.);
    nerror += inputin->get_item(&rndexp, "fields", "rndexp", fld, 0.);

    // Find the location of the randomizer height.
    int kendrnd = grid->kstart;
    while (grid->z[kendrnd] < rndz)
        ++kendrnd;

    if (rndz > grid->zsize)
    {
        master->print_error("randomizer height rndz (%f) higher than domain top (%f)\n", rndz, grid->zsize);
        return 1;
    }

    // Issue a warning if the randomization depth is larger than zero, but less than the first model level.
    if (kendrnd == grid->kstart && rndz > 0.)
        master->print_warning("randomization depth is less than the height of the first model level\n");

    for (int k=grid->kstart; k<kendrnd; ++k)
    {
        const double rndfac = std::pow((rndz-grid->z [k])/rndz, rndexp);
        for (int j=grid->jstart; j<grid->jend; ++j)
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                data[ijk] = rndfac * rndamp * ((double) std::rand() / (double) RAND_MAX - 0.5);
            }
    }

    return nerror;
}

int Fields::add_vortex_pair(Input* inputin)
{
    int nerror = 0;

    // optional parameters
    nerror += inputin->get_item(&vortexnpair, "fields", "vortexnpair", "", 0    );
    nerror += inputin->get_item(&vortexamp  , "fields", "vortexamp"  , "", 1.e-3);
    nerror += inputin->get_item(&vortexaxis , "fields", "vortexaxis" , "", "y"  );

    // add a double vortex to the initial conditions
    const double pi = std::acos((double)-1.);

    const int jj = grid->icells;
    const int kk = grid->ijcells;

    if (vortexnpair > 0)
    {
        if (vortexaxis == "y")
            for (int k=grid->kstart; k<grid->kend; ++k)
                for (int j=grid->jstart; j<grid->jend; ++j)
                    for (int i=grid->istart; i<grid->iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        u->data[ijk] +=  vortexamp*std::sin(vortexnpair*2.*pi*(grid->xh[i])/grid->xsize)*std::cos(pi*grid->z [k]/grid->zsize);
                        w->data[ijk] += -vortexamp*std::cos(vortexnpair*2.*pi*(grid->x [i])/grid->xsize)*std::sin(pi*grid->zh[k]/grid->zsize);
                    }
        else if (vortexaxis == "x")
            for (int k=grid->kstart; k<grid->kend; ++k)
                for (int j=grid->jstart; j<grid->jend; ++j)
                    for (int i=grid->istart; i<grid->iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        v->data[ijk] +=  vortexamp*std::sin(vortexnpair*2.*pi*(grid->yh[j])/grid->ysize)*std::cos(pi*grid->z [k]/grid->zsize);
                        w->data[ijk] += -vortexamp*std::cos(vortexnpair*2.*pi*(grid->y [j])/grid->ysize)*std::sin(pi*grid->zh[k]/grid->zsize);
                    }
    }

    return nerror;
}

int Fields::add_mean_prof(Input* inputin, std::string fld, double* restrict data, double offset)
{
    double proftemp[grid->kmax];

    const int jj = grid->icells;
    const int kk = grid->ijcells;

    if (inputin->get_prof(proftemp, fld, grid->kmax))
        return 1;

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                data[ijk] += proftemp[k-grid->kstart] - offset;
            }

    return 0;
}

void Fields::load(int n)
{
    const double NoOffset = 0.;

    int nerror = 0;

    for (FieldMap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
    {
        // the offset is kept at zero, otherwise bitwise identical restarts is not possible
        char filename[256];
        std::sprintf(filename, "%s.%07d", it->second->name.c_str(), n);
        master->print_message("Loading \"%s\" ... ", filename);
        if (grid->load_field3d(it->second->data, atmp["tmp1"]->data, atmp["tmp2"]->data, filename, NoOffset))
        {
            master->print_message("FAILED\n");
            ++nerror;
        }
        else
        {
            master->print_message("OK\n");
        }  
    }

    if (nerror)
        throw 1;
}

void Fields::create_stats()
{
    int nerror = 0;

    // add the profiles to te statistics
    if (stats->get_switch() == "1")
    {
        // add variables to the statistics
        stats->add_prof(u->name, u->longname, u->unit, "z" );
        stats->add_prof(v->name, v->longname, v->unit, "z" );
        stats->add_prof(w->name, w->longname, w->unit, "zh" );

        for (FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
            stats->add_prof(it->first,it->second->longname, it->second->unit, "z");

        stats->add_prof(sd["p"]->name, sd["p"]->longname, sd["p"]->unit, "z");
        std::string sn("2");
        stats->add_prof(sd["p"]->name + sn,"Moment "+ sn + " of the " + sd["p"]->longname,"(" + sd["p"]->unit + ")"+sn, "z" );
        stats->add_prof(sd["p"]->name +"w", "Turbulent flux of the " + sd["p"]->longname, sd["p"]->unit + " m s-1", "zh");
        stats->add_prof(sd["p"]->name +"grad", "Gradient of the " + sd["p"]->longname, sd["p"]->unit + " m-1", "zh");

        // CvH, shouldn't this call be in the diffusion class?
        if (model->diff->get_switch() == "smag2")
            stats->add_prof(sd["evisc"]->name, sd["evisc"]->longname, sd["evisc"]->unit, "z");

        // moments
        for (int n=2; n<5; ++n)
        {
            std::stringstream ss;
            ss << n;
            std::string sn = ss.str();
            stats->add_prof(u->name + sn,"Moment "+ sn + " of the " + u->longname,"(" + u->unit + ")"+sn, "z" );
            stats->add_prof(v->name + sn,"Moment "+ sn + " of the " + v->longname,"(" + v->unit + ")"+sn, "z" );
            stats->add_prof(w->name + sn,"Moment "+ sn + " of the " + w->longname,"(" + w->unit + ")"+sn, "zh" );
            for (FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
                stats->add_prof(it->first + sn,"Moment "+ sn + " of the " + it->second->longname,"(" + it->second->unit + ")"+sn, "z" );
        }

        // gradients
        stats->add_prof(u->name + "grad", "Gradient of the " + u->longname,"s-1","zh");
        stats->add_prof(v->name + "grad", "Gradient of the " + v->longname,"s-1","zh");
        for (FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
            stats->add_prof(it->first+"grad", "Gradient of the " + it->second->longname, it->second->unit + " m-1", "zh");

        // turbulent fluxes
        stats->add_prof("uw", "Turbulent flux of the " + u->longname, "m2 s-2", "zh");
        stats->add_prof("vw", "Turbulent flux of the " + v->longname, "m2 s-2", "zh");
        for (FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
            stats->add_prof(it->first+"w", "Turbulent flux of the " + it->second->longname, it->second->unit + " m s-1", "zh");

        // Diffusive fluxes
        stats->add_prof("udiff", "Diffusive flux of the " + u->longname, "m2 s-2", "zh");
        stats->add_prof("vdiff", "Diffusive flux of the " + v->longname, "m2 s-2", "zh");
        for (FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
            stats->add_prof(it->first+"diff", "Diffusive flux of the " + it->second->longname, it->second->unit + " m s-1", "zh");

        //Total fluxes
        stats->add_prof("uflux", "Total flux of the " + u->longname, "m2 s-2", "zh");
        stats->add_prof("vflux", "Total flux of the " + v->longname, "m2 s-2", "zh");
        for (FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
            stats->add_prof(it->first+"flux", "Total flux of the " + it->second->longname, it->second->unit + " m s-1", "zh");
    }

    if (nerror)
        throw 1;
}

void Fields::save(int n)
{
    const double NoOffset = 0.;

    int nerror = 0;
    for (FieldMap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", it->second->name.c_str(), n);
        master->print_message("Saving \"%s\" ... ", filename);

        // the offset is kept at zero, because otherwise bitwise identical restarts is not possible
        if (grid->save_field3d(it->second->data, atmp["tmp1"]->data, atmp["tmp2"]->data, filename, NoOffset))
        {
            master->print_message("FAILED\n");
            ++nerror;
        }  
        else
        {
            master->print_message("OK\n");
        }
    }

    if (nerror)
        throw 1;
}

#ifndef USECUDA
double Fields::check_momentum()
{
    return calc_momentum_2nd(u->data, v->data, w->data, grid->dz);
}
#endif

#ifndef USECUDA
double Fields::check_tke()
{
    return calc_tke_2nd(u->data, v->data, w->data, grid->dz);
}
#endif

#ifndef USECUDA
double Fields::check_mass()
{
    // CvH for now, do the mass check on the first scalar... Do we want to change this?
    FieldMap::const_iterator itProg=sp.begin();
    if (sp.begin() != sp.end())
        return calc_mass(itProg->second->data, grid->dz);
    else
        return 0.;
}
#endif

double Fields::calc_mass(double* restrict s, double* restrict dz)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    double mass = 0;

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                mass += s[ijk]*dz[k];
            }

    grid->get_sum(&mass);

    mass /= (grid->itot*grid->jtot*grid->zsize);

    return mass;
}

double Fields::calc_momentum_2nd(double* restrict u, double* restrict v, double* restrict w, double* restrict dz)
{
    using Finite_difference::O2::interp2;

    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    double momentum = 0;

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                momentum += (interp2(u[ijk], u[ijk+ii]) + interp2(v[ijk], v[ijk+jj]) + interp2(w[ijk], w[ijk+kk]))*dz[k];
            }

    grid->get_sum(&momentum);

    momentum /= (grid->itot*grid->jtot*grid->zsize);

    return momentum;
}

double Fields::calc_tke_2nd(double* restrict u, double* restrict v, double* restrict w, double* restrict dz)
{
    using Finite_difference::O2::interp2;

    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    double tke = 0;

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                tke += ( interp2(u[ijk]*u[ijk], u[ijk+ii]*u[ijk+ii]) 
                       + interp2(v[ijk]*v[ijk], v[ijk+jj]*v[ijk+jj]) 
                       + interp2(w[ijk]*w[ijk], w[ijk+kk]*w[ijk+kk]))*dz[k];
            }

    grid->get_sum(&tke);

    tke /= (grid->itot*grid->jtot*grid->zsize);
    tke *= 0.5;

    return tke;
}

void Fields::exec_cross()
{
    int nerror = 0;

    Cross* cross = model->cross;

    for (std::vector<std::string>::const_iterator it=crosssimple.begin(); it<crosssimple.end(); ++it)
        nerror += cross->cross_simple(a[*it]->data, atmp["tmp1"]->data, a[*it]->name);

    for (std::vector<std::string>::const_iterator it=crosslngrad.begin(); it<crosslngrad.end(); ++it)
        nerror += cross->cross_lngrad(a[*it]->data, atmp["tmp1"]->data, atmp["tmp2"]->data, grid->dzi4, a[*it]->name + "lngrad");

    for (std::vector<std::string>::const_iterator it=crossfluxbot.begin(); it<crossfluxbot.end(); ++it)
        nerror += cross->cross_plane(a[*it]->datafluxbot, atmp["tmp1"]->data, a[*it]->name + "fluxbot");

    for (std::vector<std::string>::const_iterator it=crossfluxtop.begin(); it<crossfluxtop.end(); ++it)
        nerror += cross->cross_plane(a[*it]->datafluxtop, atmp["tmp1"]->data, a[*it]->name + "fluxtop");

    for (std::vector<std::string>::const_iterator it=crossbot.begin(); it<crossbot.end(); ++it)
        nerror += cross->cross_plane(a[*it]->databot, atmp["tmp1"]->data, a[*it]->name + "bot");

    for (std::vector<std::string>::const_iterator it=crosstop.begin(); it<crosstop.end(); ++it)
        nerror += cross->cross_plane(a[*it]->datatop, atmp["tmp1"]->data, a[*it]->name + "top");

    if (nerror)
        throw 1;
}

void Fields::exec_dump()
{
    for (std::vector<std::string>::const_iterator it=dumplist.begin(); it<dumplist.end(); ++it)
        model->dump->save_dump(sd[*it]->data, atmp["tmp1"]->data, *it);
}
