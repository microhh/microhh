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

#include <algorithm>
#include <stdexcept>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "input.h"
#include "dump.h"
#include "finite_difference.h"
#include "budget3d.h"

namespace
{
    template<typename TF>
    void accumulate_tendency(
            TF* const restrict acc,
            const TF* const restrict tend,
            const TF sub_dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    acc[ijk] += sub_dt * tend[ijk];
                }
    }


    template<typename TF>
    void scale_field(TF* const restrict fld, const TF fac, const int ncells)
    {
        #pragma ivdep
        for (int n=0; n<ncells; ++n)
            fld[n] *= fac;
    }


    template<typename TF>
    void calc_flux_u(
            TF* const restrict flux,
            const TF* const restrict u,
            const TF* const restrict s,
            const TF* const restrict evisc,
            const TF evisc_fac,
            const TF dxi,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using namespace Finite_difference::O2;

        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    flux[ijk] = u[ijk]*interp2(s[ijk-ii], s[ijk])
                              - evisc_fac*interp2(evisc[ijk-ii], evisc[ijk])*(s[ijk]-s[ijk-ii])*dxi;
                }
    }


    template<typename TF>
    void calc_flux_v(
            TF* const restrict flux,
            const TF* const restrict v,
            const TF* const restrict s,
            const TF* const restrict evisc,
            const TF evisc_fac,
            const TF dyi,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using namespace Finite_difference::O2;

        const int jj = icells;
        const int kk = ijcells;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    flux[ijk] = v[ijk]*interp2(s[ijk-jj], s[ijk])
                              - evisc_fac*interp2(evisc[ijk-jj], evisc[ijk])*(s[ijk]-s[ijk-jj])*dyi;
                }
    }


    template<typename TF>
    void calc_flux_w(
            TF* const restrict flux,
            const TF* const restrict w,
            const TF* const restrict s,
            const TF* const restrict evisc,
            const TF evisc_fac,
            const TF* const restrict dzhi,
            const TF* const restrict flux_bot,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using namespace Finite_difference::O2;

        const int jj = icells;
        const int kk = ijcells;

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    flux[ijk] = w[ijk]*interp2(s[ijk-kk], s[ijk])
                              - evisc_fac*interp2(evisc[ijk-kk], evisc[ijk])*(s[ijk]-s[ijk-kk])*dzhi[k];
                }

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                flux[i + j*jj + kstart*kk] = flux_bot[i + j*jj];
                flux[i + j*jj + kend  *kk] = TF(0.);
            }
    }


    template<typename TF>
    void calc_flux_uw(
            TF* const restrict flux,
            const TF* const restrict u,
            const TF* const restrict w,
            const TF* const restrict evisc,
            const TF dxi,
            const TF* const restrict dzhi,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using namespace Finite_difference::O2;

        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const TF evisc_uw = TF(0.25)*(evisc[ijk-ii-kk] + evisc[ijk-ii] + evisc[ijk-kk] + evisc[ijk]);

                    flux[ijk] = interp2(u[ijk-kk], u[ijk])*interp2(w[ijk-ii], w[ijk])
                              - evisc_uw*( (w[ijk]-w[ijk-ii])*dxi + (u[ijk]-u[ijk-kk])*dzhi[k] );
                }

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                flux[i + j*jj + kstart*kk] = TF(0.);
                flux[i + j*jj + kend  *kk] = TF(0.);
            }
    }


    template<typename TF>
    void calc_flux_vw(
            TF* const restrict flux,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict evisc,
            const TF dyi,
            const TF* const restrict dzhi,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using namespace Finite_difference::O2;

        const int jj = icells;
        const int kk = ijcells;

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const TF evisc_vw = TF(0.25)*(evisc[ijk-jj-kk] + evisc[ijk-jj] + evisc[ijk-kk] + evisc[ijk]);

                    flux[ijk] = interp2(v[ijk-kk], v[ijk])*interp2(w[ijk-jj], w[ijk])
                              - evisc_vw*( (w[ijk]-w[ijk-jj])*dyi + (v[ijk]-v[ijk-kk])*dzhi[k] );
                }

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                flux[i + j*jj + kstart*kk] = TF(0.);
                flux[i + j*jj + kend  *kk] = TF(0.);
            }
    }
}

template<typename TF>
Budget3d<TF>::Budget3d(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Thermo<TF>& thermoin, Dump<TF>& dump, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), thermo(thermoin),
    elapsed_time(0.)
{
    if (fields.sd.find("eviscs") != fields.sd.end())
    {
        evisc_name = "eviscs";
        evisc_fac = TF(1.);
    }
    else if (fields.sd.find("evisc") != fields.sd.end())
    {
        evisc_name = "evisc";
        const TF tPr = inputin.get_item<TF>("diff", "tPr", "", TF(1./3.));
        evisc_fac = TF(1.)/tPr;
    }

    std::vector<std::string>& dumplist = dump.get_dumplist();

    for (auto it = dumplist.begin(); it != dumplist.end(); )
    {
        const std::string name = *it;

        if (name.size() > 5 && name.compare(name.size()-5, 5, "_tend") == 0
                && fields.at.find(name.substr(0, name.size()-5)) != fields.at.end())
        {
            const std::string varname = name.substr(0, name.size()-5);
            auto& tend = *fields.at.at(varname);

            fields.init_diagnostic_field(
                    name, "Accumulated mean tendency of " + varname,
                    tend.unit, "budget3d", tend.loc);

            tendencylist.push_back(varname);
            it = dumplist.erase(it);
        }
        else if (name == "uw" || name == "vw")
        {
            if (fields.sd.find("evisc") == fields.sd.end())
                throw std::runtime_error("[budget3d] dumping \"" + name + "\" requires an eddy viscosity field (\"evisc\")");

            fluxlist.push_back(name);
            it = dumplist.erase(it);
        }
        else if (name.size() > 1 && (name.front() == 'u' || name.front() == 'v' || name.front() == 'w')
                && (fields.sp.find(name.substr(1)) != fields.sp.end() || thermo.check_field_exists(name.substr(1))))
        {
            if (evisc_name.empty())
                throw std::runtime_error("[budget3d] dumping \"" + name + "\" requires an eddy viscosity field (\"evisc\" or \"eviscs\")");

            fluxlist.push_back(name);
            it = dumplist.erase(it);
        }
        else
            ++it;
    }

    swbudget3d = !tendencylist.empty() || !fluxlist.empty();
}

template<typename TF>
Budget3d<TF>::~Budget3d()
{
}

#ifndef USECUDA
template<typename TF>
void Budget3d<TF>::exec(const double sub_dt, const double dt, const bool in_substep)
{
    if (!swbudget3d)
        return;

    auto& gd = grid.get_grid_data();
    const TF sub_dt_tf = static_cast<TF>(sub_dt);

    for (auto& varname : tendencylist)
    {
        auto& acc  = *fields.sd.at(varname + "_tend");
        auto& tend = *fields.at.at(varname);

        accumulate_tendency(
                acc.fld.data(),
                tend.fld.data(),
                sub_dt_tf,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }

    // Sum of cB over the RK substeps is not 1; count the full time step once.
    if (!in_substep)
        elapsed_time += dt;
}
#endif

template<typename TF>
void Budget3d<TF>::prepare_dump()
{
    if (!swbudget3d || elapsed_time <= 0.)
        return;

    // At this point `fld` on the host is up to date: the generic dump/cross/stats
    // block in Model::exec() always calls fields->backward_device() first.
    auto& gd = grid.get_grid_data();
    const TF dti = TF(1.) / static_cast<TF>(elapsed_time);

    for (auto& varname : tendencylist)
    {
        auto& acc = *fields.sd.at(varname + "_tend");
        scale_field(acc.fld.data(), dti, gd.ncells);
    }
}

template<typename TF>
void Budget3d<TF>::exec_dump(Dump<TF>& dump, unsigned long iotime)
{
    if (!swbudget3d)
        return;

    for (auto& varname : tendencylist)
    {
        auto& acc = *fields.sd.at(varname + "_tend");
        dump.save_dump(acc.fld.data(), varname + "_tend", iotime);
    }

    if (fluxlist.empty())
        return;

    auto& gd = grid.get_grid_data();

    auto& u = *fields.mp.at("u");
    auto& v = *fields.mp.at("v");
    auto& w = *fields.mp.at("w");
    auto& evisc = *fields.sd.at(evisc_name);

    auto flux = fields.get_tmp();
    auto s_tmp = fields.get_tmp();

    for (auto& name : fluxlist)
    {
        if (name == "uw")
            calc_flux_uw(
                    flux->fld.data(), u.fld.data(), w.fld.data(), fields.sd.at("evisc")->fld.data(),
                    gd.dxi, gd.dzhi.data(), gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        else if (name == "vw")
            calc_flux_vw(
                    flux->fld.data(), v.fld.data(), w.fld.data(), fields.sd.at("evisc")->fld.data(),
                    gd.dyi, gd.dzhi.data(), gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        else
        {
            const char dir = name[0];
            const std::string varname = name.substr(1);

            TF* s;
            const TF* flux_bot;
            if (fields.sp.find(varname) != fields.sp.end())
            {
                auto& sp = *fields.sp.at(varname);
                s = sp.fld.data();
                flux_bot = sp.flux_bot.data();
            }
            else
            {
                thermo.get_thermo_field(*s_tmp, varname, true, true);
                s = s_tmp->fld.data();
                std::fill(s_tmp->flux_bot.begin(), s_tmp->flux_bot.end(), TF(0.));
                flux_bot = s_tmp->flux_bot.data();
            }

            if (dir == 'u')
                calc_flux_u(
                        flux->fld.data(), u.fld.data(), s, evisc.fld.data(), evisc_fac,
                        gd.dxi, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
            else if (dir == 'v')
                calc_flux_v(
                        flux->fld.data(), v.fld.data(), s, evisc.fld.data(), evisc_fac,
                        gd.dyi, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
            else
                calc_flux_w(
                        flux->fld.data(), w.fld.data(), s, evisc.fld.data(), evisc_fac,
                        gd.dzhi.data(), flux_bot,
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        }

        dump.save_dump(flux->fld.data(), name, iotime);
    }

    fields.release_tmp(flux);
    fields.release_tmp(s_tmp);
}

template<typename TF>
void Budget3d<TF>::reset()
{
    if (!swbudget3d)
        return;

    for (auto& varname : tendencylist)
    {
        auto& acc = *fields.sd.at(varname + "_tend");
        std::fill(acc.fld.begin(), acc.fld.end(), TF(0));

        #ifdef USECUDA
        fields.forward_field_device_3d(acc.fld_g, acc.fld.data());
        #endif
    }

    elapsed_time = 0.;
}


#ifdef FLOAT_SINGLE
template class Budget3d<float>;
#else
template class Budget3d<double>;
#endif
