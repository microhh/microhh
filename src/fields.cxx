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

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "field3d.h"
#include "input.h"
#include "data_block.h"
#include "defines.h"
#include "finite_difference.h"
#include "stats.h"
#include "column.h"
#include "cross.h"
#include "dump.h"
#include "diff.h"

namespace
{
    enum class Mask_type {Wplus, Wmin};

    template<typename TF, Mask_type mode>
    void calc_mask_w(TF* const restrict mask, TF* const restrict maskh, TF* const restrict maskbot,
                     int* const restrict nmask, int* const restrict nmaskh, const TF* const restrict w,
                     const int istart, const int jstart, const int kstart,
                     const int iend,   const int jend,   const int kend,
                     const int icells, const int ijcells)
    {
        int ntmp = 0;

        // Calculate mask at full levels
        for (int k=kstart; k<kend; k++)
        {
            nmask[k] = 0;
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;

                    if (mode == Mask_type::Wplus)
                        ntmp = (w[ijk] + w[ijk+ijcells]) >  0.;
                    else
                        ntmp = (w[ijk] + w[ijk+ijcells]) <= 0.;

                    nmask[k] += ntmp;
                    mask[ijk] = static_cast<TF>(ntmp);
                }
        }

        // Calculate mask at half levels
        for (int k=kstart; k<kend+1; k++)
        {
            nmaskh[k] = 0;
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;

                    if (mode == Mask_type::Wplus)
                        ntmp = w[ijk] >  0.;
                    else
                        ntmp = w[ijk] <= 0.;

                    nmaskh[k] += ntmp;
                    maskh[ijk] = static_cast<TF>(ntmp);
                }
        }

        // Set the mask for surface projected quantities
        // In this case: velocity at surface, so zero
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij  = i + j*icells;
                const int ijk = i + j*icells + kstart*ijcells;

                maskbot[ij] = maskh[ijk];
            }
    }
}


template<typename TF>
Fields<TF>::Fields(Master& masterin, Grid<TF>& gridin, Input& input) :
    master(masterin),
    grid(gridin),
    boundary_cyclic(master, grid),
    field3d_io(master, grid)
{
    calc_mean_profs = false;

    // Initialize GPU pointers
    // rhoref_g  = 0;
    // rhorefh_g = 0;

    // obligatory parameters
    visc = input.get_item<TF>("fields", "visc", "");

    // Initialize the passive scalars
    std::vector<std::string> slist = input.get_list<std::string>("fields", "slist", "", std::vector<std::string>());
    for (auto& s : slist)
    {
        init_prognostic_field(s, s, "-");
        sp.at(s)->visc = input.get_item<TF>("fields", "svisc", s);
    }

    // initialize the basic set of fields
    init_momentum_field("u", "U velocity", "m s-1");
    init_momentum_field("v", "V velocity", "m s-1");
    init_momentum_field("w", "Vertical velocity", "m s-1");
    init_diagnostic_field("p", "Pressure", "Pa");

    // Set a default of 4 temporary fields. Other classes can increase this number
    // before the init phase, where they are initialized in Fields::init()
    n_tmp_fields = 4;

    // Specify the masks that fields can provide / calculate
    available_masks.insert(available_masks.end(), {"wplus", "wmin"});

    // Remove the data from the input that is not used in run mode, to avoid warnings.
    /*
    if (master.mode == "run")
    {
        input.flag_as_used("fields", "rndamp");
        input.flag_as_used("fields", "rndexp");
        input.flag_as_used("fields", "rndseed");
        input.flag_as_used("fields", "rndz");

        input.flag_as_used("fields", "vortexnpair");
        input.flag_as_used("fields", "vortexamp"  );
        input.flag_as_used("fields", "vortexaxis" );
    }
    */
}

template<typename TF>
Fields<TF>::~Fields()
{
}

template<typename TF>
void Fields<TF>::init(Dump<TF>& dump, Cross<TF>& cross)
{
    int nerror = 0;
    // ALLOCATE ALL THE FIELDS
    // allocate the prognostic velocity fields
    for (auto& it : mp)
        nerror += it.second->init();

    // allocate the velocity tendency fields
    for (auto& it : mt)
        nerror += it.second->init();

    // allocate the prognostic scalar fields
    for (auto& it : sp)
        nerror += it.second->init();

    // allocate the scalar tendency fields
    for (auto& it : st)
        nerror += it.second->init();

    // allocate the diagnostic scalars
    for (auto& it : sd)
        nerror += it.second->init();

    // now that all classes have been able to set the minimum number of tmp fields, initialize them
    for (int i=0; i<n_tmp_fields; ++i)
        init_tmp_field();

    // allocate the tmp fields
    for (auto& tmp : atmp)
        nerror += tmp->init();

    if (nerror > 0)
        throw 1;

    // Get the grid data.
    const Grid_data<TF>& gd = grid.get_grid_data();

    rhoref .resize(gd.kcells);
    rhorefh.resize(gd.kcells);

    // \TODO Define a reference density. Needs to be replaced once anelastic is there
    // BvS: Always init rhoref at 1 for situation with e.g. thermo=0? For anelastic, overwrite it.
    std::fill(rhoref .begin(), rhoref .end(), 1.);
    std::fill(rhorefh.begin(), rhorefh.end(), 1.);

    // Create help arrays for statistics.
    umodel.resize(gd.kcells);
    vmodel.resize(gd.kcells);

    // Init the toolbox classes.
    boundary_cyclic.init();
    field3d_io.init();

    // Set up output classes
    create_dump(dump);
    create_cross(cross);


}

template<typename TF>
void Fields<TF>::create_dump(Dump<TF>& dump)
{
    // add the profiles to the columns
    if (dump.get_switch())
    {
        // Get global dump-list from cross.cxx
        std::vector<std::string> *dumplist_global = dump.get_dumplist();

        // Check if fields in dumplist are diagnostic fields, if not delete them and print warning
        std::vector<std::string>::iterator dumpvar=dumplist_global->begin();
        while (dumpvar != dumplist_global->end())
        {
            if (a.count(*dumpvar))
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

template<typename TF>
void Fields<TF>::create_cross(Cross<TF>& cross)
{

    if (cross.get_switch())
    {

        // Get global cross-list from cross.cxx
        std::vector<std::string> *crosslist_global = cross.get_crosslist();

        // Check different type of crosses and put them in their respective lists
        for (auto& it : ap)
        {
            check_added_cross(it.first, "",        crosslist_global, &cross_simple);
            check_added_cross(it.first, "lngrad",  crosslist_global, &cross_lngrad);
            check_added_cross(it.first, "bot",     crosslist_global, &cross_bot);
            check_added_cross(it.first, "top",     crosslist_global, &cross_top);
            check_added_cross(it.first, "fluxbot", crosslist_global, &cross_fluxbot);
            check_added_cross(it.first, "fluxtop", crosslist_global, &cross_fluxtop);
        }

        for (auto& it : sd)
        {
            check_added_cross(it.first, "",        crosslist_global, &cross_simple);
            check_added_cross(it.first, "lngrad",  crosslist_global, &cross_lngrad);
        }


   }

}

template<typename TF>
void Fields<TF>::check_added_cross(std::string var, std::string type, std::vector<std::string> *crosslist, std::vector<std::string> *typelist)
{
    std::vector<std::string>::iterator position;

    position = std::find(crosslist->begin(), crosslist->end(), var + type);
    if (position != crosslist->end())
    {
        // don't allow lngrad in 2nd order mode
        if (!(type == "lngrad" && grid.swspatialorder == "2"))
        {
            typelist->push_back(var);
            crosslist->erase(position);
        }
    }
}


template<typename TF>
std::shared_ptr<Field3d<TF>> Fields<TF>::get_tmp()
{
    std::shared_ptr<Field3d<TF>> tmp;

    // In case of insufficient tmp fields, allocate a new one.
    if (atmp.empty())
    {
        init_tmp_field();
        tmp = atmp.back();
        tmp->init();
    }
    else
        tmp = atmp.back();

    atmp.pop_back();
    return tmp;
}

template<typename TF>
void Fields<TF>::release_tmp(std::shared_ptr<Field3d<TF>>& tmp)
{
    atmp.push_back(std::move(tmp));
}

#ifdef USECUDA
template<typename TF>
std::shared_ptr<Field3d<TF>> Fields<TF>::get_tmp_g()
{
    std::shared_ptr<Field3d<TF>> tmp;

    // In case of insufficient tmp fields, allocate a new one.
    if (atmp_g.empty())
    {
        init_tmp_field_g();
        tmp = atmp_g.back();
        tmp->init_device();
    }
    else
        tmp = atmp_g.back();

    atmp_g.pop_back();
    return tmp;
}

template<typename TF>
void Fields<TF>::release_tmp_g(std::shared_ptr<Field3d<TF>>& tmp)
{
    atmp_g.push_back(std::move(tmp));
}
#endif

template<typename TF>
void Fields<TF>::get_mask(Field3d<TF>& mfield, Field3d<TF>& mfieldh, Stats<TF>& stats, std::string mask_name)
{
    auto& gd = grid.get_grid_data();

    if (mask_name == "wplus")
        calc_mask_w<TF, Mask_type::Wplus>(mfield.fld.data(), mfieldh.fld.data(), mfieldh.fld_bot.data(),
                                          stats.nmask.data(), stats.nmaskh.data(), mp["w"]->fld.data(),
                                          gd.istart, gd.jstart, gd.kstart, gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);
    else if (mask_name == "wmin")
        calc_mask_w<TF, Mask_type::Wmin >(mfield.fld.data(), mfieldh.fld.data(), mfieldh.fld_bot.data(),
                                          stats.nmask.data(), stats.nmaskh.data(), mp["w"]->fld.data(),
                                          gd.istart, gd.jstart, gd.kstart, gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);

    boundary_cyclic.exec(mfield.fld.data());
    boundary_cyclic.exec(mfieldh.fld.data());
    boundary_cyclic.exec_2d(mfieldh.fld_bot.data());

    master.sum(stats.nmask.data() , gd.kcells);
    master.sum(stats.nmaskh.data(), gd.kcells);
    stats.nmaskbot = stats.nmaskh[gd.kstart];
}

template<typename TF>
void Fields<TF>::exec_stats(Stats<TF>& stats, std::string mask_name, Field3d<TF>& mask_field, Field3d<TF>& mask_fieldh,
        const Diff<TF>& diff)
{
    auto& gd = grid.get_grid_data();

    Mask<TF>& m = stats.masks[mask_name];

    // Define locations
    const int uloc[] = {1,0,0};
    const int vloc[] = {0,1,0};
    const int wloc[] = {0,0,1};
    const int sloc[] = {0,0,0};

    const int uwloc[] = {1,0,1};
    const int vwloc[] = {0,1,1};

    const TF no_offset = 0.;

    // Save the area coverage of the mask
    stats.calc_area(m.profs["area" ].data.data(), sloc, stats.nmask .data());
    stats.calc_area(m.profs["areah"].data.data(), wloc, stats.nmaskh.data());

    // Start with the stats on the w location, to make the wmean known for the flux calculations
    stats.calc_mean(m.profs["w"].data.data(), mp["w"]->fld.data(), no_offset, mask_fieldh.fld.data(), stats.nmaskh.data());
    for (int n=2; n<5; ++n)
    {
        std::string sn = std::to_string(n);
        stats.calc_moment(mp["w"]->fld.data(), m.profs["w"].data.data(), m.profs["w"+sn].data.data(), n, mask_fieldh.fld.data(), stats.nmaskh.data());
    }

    // Calculate the stats on the u location
    // Interpolate the mask horizontally onto the u coordinate.
    auto mask_on_u = get_tmp();
    grid.interpolate_2nd(mask_on_u->fld.data(), mask_field.fld.data(), sloc, uloc);

    stats.calc_mean(m.profs["u"].data.data(), mp["u"]->fld.data(), grid.utrans, mask_on_u->fld.data(), stats.nmask.data());
    stats.calc_mean(umodel.data(), mp["u"]->fld.data(), no_offset, mask_on_u->fld.data(), stats.nmask.data());

    for (int n=2; n<5; ++n)
    {
        std::string sn = std::to_string(n);
        stats.calc_moment(mp["u"]->fld.data(), umodel.data(), m.profs["u"+sn].data.data(), n, mask_on_u->fld.data(), stats.nmask.data());
    }
    release_tmp(mask_on_u);

    // Interpolate the mask on half level horizontally onto the u coordinate
    auto maskh_on_u = get_tmp();
    auto tmp = get_tmp();
    grid.interpolate_2nd(maskh_on_u->fld.data(), mask_fieldh.fld.data(), wloc, uwloc);

    if (grid.swspatialorder == "2")
    {
        stats.calc_grad_2nd(mp["u"]->fld.data(), m.profs["ugrad"].data.data(), gd.dzhi.data(), maskh_on_u->fld.data(), stats.nmaskh.data());
        stats.calc_flux_2nd(mp["u"]->fld.data(), umodel.data(), mp["w"]->fld.data(), m.profs["w"].data.data(),
                            m.profs["uw"].data.data(), tmp->fld.data(), uloc, maskh_on_u->fld.data(), stats.nmaskh.data());

        if (diff.get_switch() == Diffusion_type::Diff_smag2)
        {
            stats.calc_diff_2nd(mp["u"]->fld.data(), mp["w"]->fld.data(), sd["evisc"]->fld.data(),
                                m.profs["udiff"].data.data(), gd.dzhi.data(),
                                mp["u"]->flux_bot.data(), mp["u"]->flux_top.data(), 1., uloc,
                                maskh_on_u->fld.data(), stats.nmaskh.data());
        }
        else
        {
            stats.calc_diff_2nd(mp["u"]->fld.data(), m.profs["udiff"].data.data(), gd.dzhi.data(), visc, uloc, maskh_on_u->fld.data(), stats.nmaskh.data());
        }

    }
    else if (grid.swspatialorder == "4")
    {
        master.print_message("4th order statistics not yet implemented in fields...\n");
        //stats.calc_grad_4th(u->data, m.profs["ugrad"].data, grid.dzhi4, uloc,
        //                    atmp["tmp1"]->data, stats.nmaskh);
        //stats.calc_flux_4th(u->data, w->data, m.profs["uw"].data, atmp["tmp2"]->data, uloc,
        //                    atmp["tmp1"]->data, stats.nmaskh);
        //stats.calc_diff_4th(u->data, m.profs["udiff"].data, grid.dzhi4, visc, uloc,
        //                    atmp["tmp1"]->data, stats.nmaskh);
    }
    release_tmp(maskh_on_u);

    // Calculate the statistics on the v location
    auto mask_on_v = get_tmp();

    grid.interpolate_2nd(mask_on_v->fld.data(), mask_field.fld.data(), sloc, vloc);
    stats.calc_mean(m.profs["v"].data.data(), mp["v"]->fld.data(), grid.vtrans, mask_on_v->fld.data(), stats.nmask.data());
    stats.calc_mean(vmodel.data(), mp["v"]->fld.data(), no_offset, mask_on_v->fld.data(), stats.nmask.data());

    for (int n=2; n<5; ++n)
    {
        std::string sn = std::to_string(n);
        stats.calc_moment(mp["v"]->fld.data(), vmodel.data(), m.profs["v"+sn].data.data(), n, mask_on_v->fld.data(), stats.nmask.data());
    }
    release_tmp(mask_on_v);


    // interpolate the mask on half level horizontally onto the u coordinate
    auto maskh_on_v = get_tmp();
    grid.interpolate_2nd(maskh_on_v->fld.data(), mask_fieldh.fld.data(), wloc, vwloc);

    if (grid.swspatialorder == "2")
    {
        stats.calc_grad_2nd(mp["v"]->fld.data(), m.profs["vgrad"].data.data(), gd.dzhi.data(), maskh_on_v->fld.data(), stats.nmaskh.data());
        stats.calc_flux_2nd(mp["v"]->fld.data(), vmodel.data(), mp["w"]->fld.data(), m.profs["w"].data.data(),
                            m.profs["vw"].data.data(), tmp->fld.data(), vloc, maskh_on_v->fld.data(), stats.nmaskh.data());

        if (diff.get_switch() == Diffusion_type::Diff_smag2)
        {
            stats.calc_diff_2nd(mp["v"]->fld.data(), mp["w"]->fld.data(), sd["evisc"]->fld.data(),
                                m.profs["vdiff"].data.data(), gd.dzhi.data(),
                                mp["v"]->flux_bot.data(), mp["v"]->flux_top.data(), 1., vloc,
                                maskh_on_v->fld.data(), stats.nmaskh.data());
        }
        else
        {
            stats.calc_diff_2nd(mp["v"]->fld.data(), m.profs["vdiff"].data.data(), gd.dzhi.data(), visc, vloc, maskh_on_v->fld.data(), stats.nmaskh.data());
        }

    }
    else if (grid.swspatialorder == "4")
    {
        //stats.calc_grad_4th(v->data, m.profs["vgrad"].data, grid.dzhi4, vloc,
        //                    atmp["tmp1"]->data, stats.nmaskh);
        //stats.calc_flux_4th(v->data, w->data, m.profs["vw"].data, atmp["tmp2"]->data, vloc,
        //                    atmp["tmp1"]->data, stats.nmaskh);
        //stats.calc_diff_4th(v->data, m.profs["vdiff"].data, grid.dzhi4, visc, vloc,
        //                    atmp["tmp1"]->data, stats.nmaskh);
    }
    release_tmp(maskh_on_v);

    // calculate stats for the prognostic scalars
    //Diff_smag_2 *diffptr = static_cast<Diff_smag_2 *>(model->diff);
    for (auto& it : sp)
    {
        stats.calc_mean(m.profs[it.first].data.data(), it.second->fld.data(), no_offset, mask_field.fld.data(), stats.nmask.data());

        for (int n=2; n<5; ++n)
        {
            std::string sn = std::to_string(n);
            stats.calc_moment(it.second->fld.data(), m.profs[it.first].data.data(), m.profs[it.first+sn].data.data(), n, mask_field.fld.data(), stats.nmask.data());
        }

        if (grid.swspatialorder == "2")
        {
            stats.calc_grad_2nd(it.second->fld.data(), m.profs[it.first+"grad"].data.data(), gd.dzhi.data(), mask_fieldh.fld.data(), stats.nmaskh.data());
            stats.calc_flux_2nd(it.second->fld.data(), m.profs[it.first].data.data(), mp["w"]->fld.data(), m.profs["w"].data.data(),
                                m.profs[it.first+"w"].data.data(), tmp->fld.data(), sloc, mask_fieldh.fld.data(), stats.nmaskh.data());

            if (diff.get_switch() == Diffusion_type::Diff_smag2)
            {
                stats.calc_diff_2nd(it.second->fld.data(), mp["w"]->fld.data(), sd["evisc"]->fld.data(),
                                    m.profs[it.first+"diff"].data.data(), gd.dzhi.data(),
                                    it.second->flux_bot.data(), it.second->flux_top.data(), diff.tPr, sloc,
                                    mask_fieldh.fld.data(), stats.nmaskh.data());
            }
            else
            {
                stats.calc_diff_2nd(it.second->fld.data(), m.profs[it.first+"diff"].data.data(), gd.dzhi.data(),
                                    it.second->visc, sloc, mask_fieldh.fld.data(), stats.nmaskh.data());
            }
        }
        else if (grid.swspatialorder == "4")
        {
            //stats.calc_grad_4th(it.second->data, m.profs[it.first+"grad"].data, grid.dzhi4, sloc,
            //                    atmp["tmp4"]->data, stats.nmaskh);
            //stats.calc_flux_4th(it.second->data, w->data, m.profs[it.first+"w"].data, atmp["tmp1"]->data, sloc,
            //                    atmp["tmp4"]->data, stats.nmaskh);
            //stats.calc_diff_4th(it.second->data, m.profs[it.first+"diff"].data, grid.dzhi4, it.second->visc, sloc,
            //                    atmp["tmp4"]->data, stats.nmaskh);
        }
    }

    // Calculate pressure statistics
    stats.calc_mean(m.profs["p"].data.data(), sd["p"]->fld.data(), no_offset, mask_field.fld.data(), stats.nmask.data());
    stats.calc_moment(sd["p"]->fld.data(), m.profs["p"].data.data(), m.profs["p2"].data.data(), 2, mask_field.fld.data(), stats.nmask.data());

    if (grid.swspatialorder == "2")
    {
        stats.calc_grad_2nd(sd["p"]->fld.data(), m.profs["pgrad"].data.data(), gd.dzhi.data(), mask_fieldh.fld.data(), stats.nmaskh.data());
        stats.calc_flux_2nd(sd["p"]->fld.data(), m.profs["p"].data.data(), mp["w"]->fld.data(), m.profs["w"].data.data(),
                            m.profs["pw"].data.data(), tmp->fld.data(), sloc, mask_fieldh.fld.data(), stats.nmaskh.data());
    }
    else if (grid.swspatialorder == "4")
    {
        //stats.calc_grad_4th(sd["p"]->data, m.profs["pgrad"].data, grid.dzhi4, sloc,
        //                     atmp["tmp4"]->data, stats.nmaskh);
        //stats.calc_flux_4th(sd["p"]->data, w->data, m.profs["pw"].data, atmp["tmp1"]->data, sloc,
        //                     atmp["tmp4"]->data, stats.nmaskh);
    }

    // Calculate the total fluxes
    stats.add_fluxes(m.profs["uflux"].data.data(), m.profs["uw"].data.data(), m.profs["udiff"].data.data());
    stats.add_fluxes(m.profs["vflux"].data.data(), m.profs["vw"].data.data(), m.profs["vdiff"].data.data());
    for (auto& it : sp)
        stats.add_fluxes(m.profs[it.first+"flux"].data.data(), m.profs[it.first+"w"].data.data(), m.profs[it.first+"diff"].data.data());

    if (diff.get_switch() == Diffusion_type::Diff_smag2)
        stats.calc_mean(m.profs["evisc"].data.data(), sd["evisc"]->fld.data(), no_offset, mask_field.fld.data(), stats.nmask.data());

    release_tmp(tmp);
}

/*
void Fields::set_calc_mean_profs(bool sw)
{
    calc_mean_profs = sw;
}

void Fields::set_minimum_tmp_fields(int n)
{
    n_tmp_fields = std::max(n_tmp_fields, n);
}
*/

template<typename TF>
void Fields<TF>::init_momentum_field(std::string fldname, std::string longname, std::string unit)
{
    if (mp.find(fldname) != mp.end())
    {
        master.print_error("\"%s\" already exists\n", fldname.c_str());
        throw 1;
    }

    // Add a new prognostic momentum variable.
    mp[fldname] = std::make_shared<Field3d<TF>>(master, grid, fldname, longname, unit);

    // Add a new tendency for momentum variable.
    std::string fldtname  = fldname + "t";
    std::string tunit     = unit + "s-1";
    std::string tlongname = "Tendency of " + longname;
    mt[fldname] = std::make_shared<Field3d<TF>>(master, grid, fldtname, tlongname, tunit);

    // Add the prognostic variable and its tendency to the collection
    // of all fields and tendencies.
    a [fldname] = mp[fldname];
    ap[fldname] = mp[fldname];
    at[fldname] = mt[fldname];
}

template<typename TF>
void Fields<TF>::init_prognostic_field(std::string fldname, std::string longname, std::string unit)
{
    if (sp.find(fldname)!=sp.end())
    {
        master.print_error("\"%s\" already exists\n", fldname.c_str());
        throw 1;
    }

    // add a new scalar variable
    sp[fldname] = std::make_shared<Field3d<TF>>(master, grid, fldname, longname, unit);

    // add a new tendency for scalar variable
    std::string fldtname  = fldname + "t";
    std::string tlongname = "Tendency of " + longname;
    std::string tunit     = unit + "s-1";
    st[fldname] = std::make_shared<Field3d<TF>>(master, grid, fldtname, tlongname, tunit);

    // add the prognostic variable and its tendency to the collection
    // of all fields and tendencies
    a [fldname] = sp[fldname];
    ap[fldname] = sp[fldname];
    at[fldname] = st[fldname];
}

template<typename TF>
void Fields<TF>::init_diagnostic_field(std::string fldname,std::string longname, std::string unit)
{
    if (sd.find(fldname)!=sd.end())
    {
        master.print_error("\"%s\" already exists\n", fldname.c_str());
        throw 1;
    }

    sd[fldname] = std::make_shared<Field3d<TF>>(master, grid, fldname, longname, unit);
    a [fldname] = sd[fldname];
}

template<typename TF>
void Fields<TF>::init_tmp_field()
{
    static int ntmp = 0;
    ++ntmp;
    std::string fldname = "tmp" + std::to_string(ntmp);
    std::string longname = "";
    std::string unit = "";

    std::string message = "Allocating temporary field: " + fldname;
    master.print_message(message);
    atmp.push_back(std::make_shared<Field3d<TF>>(master, grid, fldname, longname, unit));
}

#ifdef USECUDA
template<typename TF>
void Fields<TF>::init_tmp_field_g()
{
    static int ntmp = 0;
    ++ntmp;
    std::string fldname = "tmp_gpu" + std::to_string(ntmp);
    std::string longname = "";
    std::string unit = "";

    std::string message = "Allocating temporary field: " + fldname;
    master.print_message(message);
    atmp_g.push_back(std::make_shared<Field3d<TF>>(master, grid, fldname, longname, unit));
}
#endif


template<typename TF>
void Fields<TF>::create(Input& inputin, Data_block& profs)
{
    // Randomize the momentum
    randomize(inputin, "u", mp.at("u")->fld.data());
    randomize(inputin, "w", mp.at("w")->fld.data());

    // Only add perturbation to v in case of a 3d run.
    const Grid_data<TF>& gd = grid.get_grid_data();
    if (gd.jtot > 1)
        randomize(inputin, "v", mp.at("v")->fld.data());

    // Randomize the scalars
    for (auto& it : sp)
        randomize(inputin, it.first, it.second->fld.data());

    // Add Vortices
    add_vortex_pair(inputin);

    // Add the mean profiles to the fields
    add_mean_profs(profs);

    /*
    nerror += add_mean_prof(inputin, "u", mp["u"]->data, grid.utrans);
    nerror += add_mean_prof(inputin, "v", mp["v"]->data, grid.vtrans);

    for (auto& it : sp)
        nerror += add_mean_prof(inputin, it.first, it.second->data, 0.);
    */

    // Make sure the boundaries of w are zero. Non-zero initialization is not recoverable.
    int lbot = gd.kstart*gd.ijcells;
    int ltop = gd.kend  *gd.ijcells;

    for (int l=0; l<gd.ijcells; ++l)
    {
        mp.at("w")->fld[lbot+l] = 0.;
        mp.at("w")->fld[ltop+l] = 0.;
    }
}

template<typename TF>
void Fields<TF>::randomize(Input& input, std::string fld, TF* const restrict data)
{
    // Set mpiid as random seed to avoid having the same field at all procs
    int static seed = 0;

    if (!seed)
    {
        seed = input.get_item<int>("fields", "rndseed", "", 0);
        seed += master.get_mpiid() + 2;
        std::srand(seed);
    }

    const Grid_data<TF>& gd = grid.get_grid_data();

    const int jj = gd.icells;
    const int kk = gd.ijcells;

    // Look up the specific randomizer variables.
    rndamp = input.get_item<TF>("fields", "rndamp", fld, 0.);
    rndz   = input.get_item<TF>("fields", "rndz"  , fld, 0.);
    rndexp = input.get_item<TF>("fields", "rndexp", fld, 0.);

    if (rndz > gd.zsize)
    {
        master.print_error("randomizer height rndz (%f) higher than domain top (%f)\n", rndz, gd.zsize);
        throw std::runtime_error("Randomizer error");
    }

    // Find the location of the randomizer height.
    int kendrnd = gd.kstart;
    while (gd.z[kendrnd] < rndz)
        ++kendrnd;

    // Issue a warning if the randomization depth is larger than zero, but less than the first model level.
    if (kendrnd == gd.kstart && rndz > 0.)
        master.print_warning("randomization depth is less than the height of the first model level\n");

    for (int k=gd.kstart; k<kendrnd; ++k)
    {
        const TF rndfac = std::pow((rndz-gd.z[k])/rndz, rndexp);
        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                data[ijk] = rndfac * rndamp * ((TF) std::rand() / (TF) RAND_MAX - 0.5);
            }
    }
}

namespace
{
    template<typename TF>
    void add_mean_prof_to_field(TF* restrict const data,
                                const TF* restrict const dataprof,
                                const TF offset,
                                const int istart, const int iend,
                                const int jstart, const int jend,
                                const int kstart, const int kend,
                                const int jj, const int kk)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    data[ijk] += dataprof[k-kstart] - offset;
                }
    }
}

template<typename TF>
void Fields<TF>::add_mean_profs(Data_block& profs)
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    std::vector<TF> prof(gd.ktot);

    profs.get_vector(prof, "u", gd.ktot, 0, 0);
    add_mean_prof_to_field<TF>(mp["u"]->fld.data(), prof.data(), grid.utrans,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    profs.get_vector(prof, "v", gd.ktot, 0, 0);
    add_mean_prof_to_field<TF>(mp["v"]->fld.data(), prof.data(), grid.vtrans,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    for (auto& f : sp)
    {
        profs.get_vector(prof, f.first, gd.ktot, 0, 0);
        add_mean_prof_to_field<TF>(f.second->fld.data(), prof.data(), 0.,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }
}


template<typename TF>
void Fields<TF>::add_vortex_pair(Input& inputin)
{
    auto& gd = grid.get_grid_data();

    // Optional parameters.
    vortexnpair = inputin.get_item<int>        ("fields", "vortexnpair", "", 0    );
    vortexamp   = inputin.get_item<TF>         ("fields", "vortexamp"  , "", 1.e-3);
    vortexaxis  = inputin.get_item<std::string>("fields", "vortexaxis" , "", "y"  );

    // Add a double vortex to the initial conditions.
    const double pi = std::acos((double)-1.);

    if (vortexnpair > 0)
    {
        if (vortexaxis == "y")
            for (int k=gd.kstart; k<gd.kend; ++k)
                for (int j=gd.jstart; j<gd.jend; ++j)
                    for (int i=gd.istart; i<gd.iend; ++i)
                    {
                        const int ijk = i + j*gd.icells + k*gd.ijcells;
                        mp["u"]->fld[ijk] +=  vortexamp*std::sin(vortexnpair*2.*pi*(gd.xh[i])/gd.xsize)*std::cos(pi*gd.z [k]/gd.zsize);
                        mp["w"]->fld[ijk] += -vortexamp*std::cos(vortexnpair*2.*pi*(gd.x [i])/gd.xsize)*std::sin(pi*gd.zh[k]/gd.zsize);
                    }
        else if (vortexaxis == "x")
            for (int k=gd.kstart; k<gd.kend; ++k)
                for (int j=gd.jstart; j<gd.jend; ++j)
                    for (int i=gd.istart; i<gd.iend; ++i)
                    {
                        const int ijk = i + j*gd.icells + k*gd.ijcells;
                        mp["v"]->fld[ijk] +=  vortexamp*std::sin(vortexnpair*2.*pi*(gd.yh[j])/gd.ysize)*std::cos(pi*gd.z [k]/gd.zsize);
                        mp["w"]->fld[ijk] += -vortexamp*std::cos(vortexnpair*2.*pi*(gd.y [j])/gd.ysize)*std::sin(pi*gd.zh[k]/gd.zsize);
                    }
    }
}

//int Fields::add_mean_prof(Input* inputin, std::string fld, double* restrict data, double offset)
//{
//    double proftemp[grid.kmax];
//
//    const int jj = grid.icells;
//    const int kk = grid.ijcells;
//
//    if (input.get_prof(proftemp, fld, grid.kmax))
//        return 1;
//
//    for (int k=grid.kstart; k<grid.kend; ++k)
//        for (int j=grid.jstart; j<grid.jend; ++j)
//            for (int i=grid.istart; i<grid.iend; ++i)
//            {
//                const int ijk = i + j*jj + k*kk;
//                data[ijk] += proftemp[k-grid.kstart] - offset;
//            }
//
//    return 0;
//}

template <typename TF>
void Fields<TF>::create_stats(Stats<TF>& stats, const Diff<TF>& diff)
{
    // Add the profiles to te statistics
    if (stats.get_switch())
    {
        // Mean velocity compontents
        stats.add_prof(ap["u"]->name, ap["u"]->longname, ap["u"]->unit, "z" );
        stats.add_prof(ap["v"]->name, ap["v"]->longname, ap["v"]->unit, "z" );
        stats.add_prof(ap["w"]->name, ap["w"]->longname, ap["w"]->unit, "zh");

        // Mean prognostic scalars
        for (auto& it : sp)
            stats.add_prof(it.first, it.second->longname, it.second->unit, "z");

        // Pressure with its variance, fluxes and gradients
        stats.add_prof(sd["p"]->name, sd["p"]->longname, sd["p"]->unit, "z");
        std::string sn("2");
        stats.add_prof(sd["p"]->name + sn,    "Moment "+ sn + " of the " + sd["p"]->longname,"(" + sd["p"]->unit + ")"+sn, "z" );
        stats.add_prof(sd["p"]->name +"w",    "Turbulent flux of the "   + sd["p"]->longname, sd["p"]->unit + " m s-1", "zh");
        stats.add_prof(sd["p"]->name +"grad", "Gradient of the "         + sd["p"]->longname, sd["p"]->unit + " m-1", "zh");

//        // CvH, shouldn't this call be in the diffusion class? BvS: yes ;-)
        if (diff.get_switch() == Diffusion_type::Diff_smag2)
            stats.add_prof(sd["evisc"]->name, sd["evisc"]->longname, sd["evisc"]->unit, "z");

        // Add the second up to fourth moments of the velocity and scalars
        for (int n=2; n<5; ++n)
        {
            std::string sn = std::to_string(n);
            //std::stringstream ss;
            //ss << n;
            //std::string sn = ss.str();
            stats.add_prof(ap["u"]->name + sn, "Moment "+ sn + " of the " + ap["u"]->longname,"(" + ap["u"]->unit + ")"+sn, "z" );
            stats.add_prof(ap["v"]->name + sn, "Moment "+ sn + " of the " + ap["v"]->longname,"(" + ap["v"]->unit + ")"+sn, "z" );
            stats.add_prof(ap["w"]->name + sn, "Moment "+ sn + " of the " + ap["w"]->longname,"(" + ap["w"]->unit + ")"+sn, "zh" );
            for (auto& it : sp)
                stats.add_prof(it.first + sn, "Moment "+ sn + " of the " + it.second->longname,"(" + it.second->unit + ")"+sn, "z" );
        }

        // Gradients
        stats.add_prof("ugrad", "Gradient of the " + ap["u"]->longname, "s-1", "zh");
        stats.add_prof("vgrad", "Gradient of the " + ap["v"]->longname, "s-1", "zh");
        for (auto& it : sp)
            stats.add_prof(it.first+"grad", "Gradient of the " + it.second->longname, it.second->unit + " m-1", "zh");

        // Turbulent fluxes
        stats.add_prof("uw", "Turbulent flux of the " + ap["u"]->longname, "m2 s-2", "zh");
        stats.add_prof("vw", "Turbulent flux of the " + ap["v"]->longname, "m2 s-2", "zh");
        for (auto& it : sp)
            stats.add_prof(it.first+"w", "Turbulent flux of the " + it.second->longname, it.second->unit + " m s-1", "zh");

        // Diffusive fluxes
        stats.add_prof("udiff", "Diffusive flux of the " + ap["u"]->longname, "m2 s-2", "zh");
        stats.add_prof("vdiff", "Diffusive flux of the " + ap["v"]->longname, "m2 s-2", "zh");
        for (auto& it : sp)
            stats.add_prof(it.first+"diff", "Diffusive flux of the " + it.second->longname, it.second->unit + " m s-1", "zh");

        // Total fluxes
        stats.add_prof("uflux", "Total flux of the " + ap["u"]->longname, "m2 s-2", "zh");
        stats.add_prof("vflux", "Total flux of the " + ap["v"]->longname, "m2 s-2", "zh");
        for (auto& it : sp)
            stats.add_prof(it.first+"flux", "Total flux of the " + it.second->longname, it.second->unit + " m s-1", "zh");
    }
}

template<typename TF>
void Fields<TF>::create_column(Column<TF>& column)
{

    // add the profiles to the columns
    if (column.get_switch())
    {
        // add variables to the statistics
        column.add_prof(ap["u"]->name, ap["u"]->longname, ap["u"]->unit, "z" );
        column.add_prof(ap["v"]->name, ap["v"]->longname, ap["v"]->unit, "z" );
        column.add_prof(ap["w"]->name, ap["w"]->longname, ap["w"]->unit, "zh" );

        for (auto& it : sp)
            column.add_prof(it.first,it.second->longname, it.second->unit, "z");

        column.add_prof(sd["p"]->name, sd["p"]->longname, sd["p"]->unit, "z");
    }

}

template<typename TF>
void Fields<TF>::save(int n)
{
    const TF no_offset = 0.;

    auto tmp1 = get_tmp();
    auto tmp2 = get_tmp();

    int nerror = 0;
    for (auto& f : ap)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", f.second->name.c_str(), n);
        master.print_message("Saving \"%s\" ... ", filename);

        // The offset is kept at zero, because otherwise bitwise identical restarts are not possible.
        if (field3d_io.save_field3d(f.second->fld.data(), tmp1->fld.data(), tmp2->fld.data(),
                    filename, no_offset))
        {
            master.print_message("FAILED\n");
            ++nerror;
        }
        else
        {
            master.print_message("OK\n");
        }
    }

    release_tmp(tmp1);
    release_tmp(tmp2);

    if (nerror)
        throw 1;
}

template<typename TF>
void Fields<TF>::load(int n)
{
    const TF no_offset = 0.;

    auto tmp1 = get_tmp();
    auto tmp2 = get_tmp();

    int nerror = 0;

    for (auto& f : ap)
    {
        // The offset is kept at zero, otherwise bitwise identical restarts is not possible.
        char filename[256];
        std::sprintf(filename, "%s.%07d", f.second->name.c_str(), n);
        master.print_message("Loading \"%s\" ... ", filename);

        if (field3d_io.load_field3d(f.second->fld.data(), tmp1->fld.data(), tmp2->fld.data(),
                    filename, no_offset))
        {
            master.print_message("FAILED\n");
            ++nerror;
        }
        else
        {
            master.print_message("OK\n");
        }
    }

    release_tmp(tmp1);
    release_tmp(tmp2);

    if (nerror)
        throw 1;
}

/*
#ifndef USECUDA
double Fields::check_momentum()
{
    return calc_momentum_2nd(u->data, v->data, w->data, grid.dz);
}
#endif

#ifndef USECUDA
double Fields::check_tke()
{
    return calc_tke_2nd(u->data, v->data, w->data, grid.dz);
}
#endif

#ifndef USECUDA
double Fields::check_mass()
{
    // CvH for now, do the mass check on the first scalar... Do we want to change this?
    Field_map::const_iterator itProg=sp.begin();
    if (sp.begin() != sp.end())
        return calc_mass(itProg->second->data, grid.dz);
    else
        return 0.;
}
#endif

double Fields::calc_mass(double* restrict s, double* restrict dz)
{
    const int jj = grid.icells;
    const int kk = grid.ijcells;

    double mass = 0;

    for (int k=grid.kstart; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                mass += s[ijk]*dz[k];
            }

    grid.get_sum(&mass);

    mass /= (grid.itot*grid.jtot*grid.zsize);

    return mass;
}

double Fields::calc_momentum_2nd(double* restrict u, double* restrict v, double* restrict w, double* restrict dz)
{
    using Finite_difference::O2::interp2;

    const int ii = 1;
    const int jj = grid.icells;
    const int kk = grid.ijcells;

    double momentum = 0;

    for (int k=grid.kstart; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                momentum += (interp2(u[ijk], u[ijk+ii]) + interp2(v[ijk], v[ijk+jj]) + interp2(w[ijk], w[ijk+kk]))*dz[k];
            }

    grid.get_sum(&momentum);

    momentum /= (grid.itot*grid.jtot*grid.zsize);

    return momentum;
}

double Fields::calc_tke_2nd(double* restrict u, double* restrict v, double* restrict w, double* restrict dz)
{
    using Finite_difference::O2::interp2;

    const int ii = 1;
    const int jj = grid.icells;
    const int kk = grid.ijcells;

    double tke = 0;

    for (int k=grid.kstart; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                tke += ( interp2(u[ijk]*u[ijk], u[ijk+ii]*u[ijk+ii])
                       + interp2(v[ijk]*v[ijk], v[ijk+jj]*v[ijk+jj])
                       + interp2(w[ijk]*w[ijk], w[ijk+kk]*w[ijk+kk]))*dz[k];
            }

    grid.get_sum(&tke);

    tke /= (grid.itot*grid.jtot*grid.zsize);
    tke *= 0.5;

    return tke;
}
*/
template<typename TF>
void Fields<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    for (auto& it : cross_simple)
        cross.cross_simple(a.at(it)->fld.data(), a.at(it)->name, iotime);

    for (auto& it : cross_lngrad)
        cross.cross_lngrad(a.at(it)->fld.data(), a.at(it)->name+"lngrad", iotime);

    for (auto& it : cross_fluxbot)
        cross.cross_plane(a.at(it)->flux_bot.data(), a.at(it)->name+"fluxbot", iotime);

    for (auto& it : cross_fluxtop)
        cross.cross_plane(a.at(it)->flux_top.data(), a.at(it)->name+"fluxtop", iotime);

    for (auto& it : cross_bot)
        cross.cross_plane(a.at(it)->fld_bot.data(), a.at(it)->name+"bot", iotime);

    for (auto& it : cross_top)
        cross.cross_plane(a.at(it)->fld_top.data(), a.at(it)->name+"top", iotime);
}

template<typename TF>
void Fields<TF>::exec_dump(Dump<TF>& dump, unsigned long iotime)
{
    for (auto& it : dumplist)
        dump.save_dump(a.at(it)->fld.data(), a.at(it)->name, iotime);
}


template<typename TF>
void Fields<TF>::exec_column(Column<TF>& column)
{
    const double no_offset = 0.;

    column.calc_column("u",mp["u"]->fld.data(), grid.utrans);
    column.calc_column("v",mp["v"]->fld.data(), grid.vtrans);
    column.calc_column("w",mp["w"]->fld.data(), no_offset);
    for (auto& it : sp)
    {
        column.calc_column(it.first, it.second->fld.data(), no_offset);
    }
    column.calc_column("p", sd["p"]->fld.data(), no_offset);
}

template<typename TF>
bool Fields<TF>::has_mask(std::string mask_name)
{
    if (std::find(available_masks.begin(), available_masks.end(), mask_name) != available_masks.end())
        return true;
    else
        return false;
}

template class Fields<double>;
template class Fields<float>;
