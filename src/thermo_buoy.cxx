/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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
#include <cmath>
#include <iostream>
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "thermo_buoy.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "master.h"


namespace
{
    using Finite_difference::O2::interp2;
    using Finite_difference::O4::interp4c;

    template<typename TF>
    void calc_buoyancy(TF* restrict b, TF* restrict bin, const int ncells)
    {
        for (int n=0; n<ncells; ++n)
            b[n] = bin[n];
    }

    template<typename TF>
    void calc_N2(TF* const restrict N2, const TF* const restrict b, const TF bg_n2, const TF* const restrict dzi,
                 const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                 const int icells, const int ijcells, const int kcells)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    N2[ijk] = TF(0.5)*(b[ijk+ijcells] - b[ijk-ijcells])*dzi[k] + bg_n2;
                }
    }

    template<typename TF>
    void calc_buoyancy_bot(
            TF* restrict b  , TF* restrict bbot,
            TF* restrict bin, TF* restrict binbot,
            const int icells, const int jcells, const int ijcells, const int kstart)
    {
        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = i + j*icells + kstart*ijcells;
                bbot[ij] = binbot[ij];
                b[ijk]   = bin[ijk];
            }
    }

    template<typename TF>
    void calc_buoyancy_fluxbot(
            TF* restrict bfluxbot, TF* restrict binfluxbot, const int icells, const int jcells)
    {
        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ij  = i + j*icells;
                bfluxbot[ij] = binfluxbot[ij];
            }
    }

    template<typename TF>
    void calc_buoyancy_tend_2nd(
            TF* restrict wt, TF* restrict b,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    wt[ijk] += interp2(b[ijk-ijcells], b[ijk]);
                }
    }

    template<typename TF>
    void calc_buoyancy_tend_u_2nd(
            TF* restrict ut, TF* restrict b, const TF alpha,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const TF sinalpha = std::sin(alpha);

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    ut[ijk] += sinalpha * interp2(b[ijk-1], b[ijk]);
                }
    }

    template<typename TF>
    void calc_buoyancy_tend_w_2nd(
            TF* restrict wt, TF* restrict b, const TF alpha,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const TF cosalpha = std::cos(alpha);

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    wt[ijk] += cosalpha * interp2(b[ijk-ijcells], b[ijk]);
                }
    }

    template<typename TF>
    void calc_buoyancy_tend_b_2nd(
            TF* restrict bt, TF* restrict u, TF* restrict w, const TF alpha, const TF n2, const TF utrans,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const TF sinalpha = std::sin(alpha);
        const TF cosalpha = std::cos(alpha);

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    bt[ijk] -= n2 * ( sinalpha * (interp2(u[ijk], u[ijk+1]) + utrans)
                                    + cosalpha *  interp2(w[ijk], w[ijk+ijcells]) );
                }
    }

    template<typename TF>
    void calc_buoyancy_tend_4th(
            TF* restrict wt, TF* restrict b,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int ijcells)
    {
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk1;
                    wt[ijk] += interp4c(b[ijk-kk2], b[ijk-kk1], b[ijk], b[ijk+kk1]);
                }
    }

    template<typename TF>
    void calc_buoyancy_tend_u_4th(
            TF* restrict ut, TF* restrict b, const TF alpha,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const TF sinalpha = std::sin(alpha);

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    ut[ijk] += sinalpha * interp4c(b[ijk-2], b[ijk-1], b[ijk], b[ijk+1]);
                }
    }

    template<typename TF>
    void calc_buoyancy_tend_w_4th(
            TF* restrict wt, TF* restrict b, const TF alpha,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj  = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        const TF cosalpha = std::cos(alpha);

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk1;
                    wt[ijk] += cosalpha * interp4c(b[ijk-kk2], b[ijk-kk1], b[ijk], b[ijk+kk1]);
                }
    }
    template<typename TF>
    void calc_buoyancy_tend_b_4th(
            TF* restrict bt, TF* restrict u, TF* restrict w, const TF alpha, const TF n2, const TF utrans,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int jj  = 1*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        const TF sinalpha = std::sin(alpha);
        const TF cosalpha = std::cos(alpha);

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk1;
                    bt[ijk] -= n2 * ( sinalpha * (interp4c(u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2]) + utrans)
                                    + cosalpha *  interp4c(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]) );
                }
    }

    template<typename TF>
    void calc_baroclinic_2nd(
            TF* const restrict bt, const TF* const restrict v,
            const TF dbdy_ls,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    bt[ijk] -= dbdy_ls * interp2(v[ijk], v[ijk+jj]);
                }
    }

    template<typename TF>
    void calc_baroclinic_4th(
            TF* const restrict bt, const TF* const restrict v,
            const TF dbdy_ls,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    bt[ijk] -= dbdy_ls * interp4c(v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2]);
                }
    }

    template<typename TF>
    int calc_zi(const TF* const restrict fldmean, const int kstart, const int kend, const int plusminus)
    {
        TF maxgrad = 0.;
        TF grad = 0.;
        int kinv = kstart;
        for (int k=kstart+1; k<kend; ++k)
        {
            grad = plusminus * (fldmean[k] - fldmean[k-1]);
            if (grad > maxgrad)
            {
                maxgrad = grad;
                kinv = k;
            }
        }
        return kinv;
    }
}

template<typename TF>
Thermo_buoy<TF>::Thermo_buoy(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
Thermo<TF>(masterin, gridin, fieldsin, inputin)
{
    auto& gd = grid.get_grid_data();
    swthermo = "buoy";

    const std::string group_name = "thermo";

    fields.init_prognostic_field("b", "Buoyancy", "m s-2", group_name, gd.sloc);

    bs.alpha = inputin.get_item<TF>("thermo", "alpha", "", 0.);
    bs.n2 = inputin.get_item<TF>("thermo", "N2"   , "", 0.);
    fields.sp.at("b")->visc = inputin.get_item<TF>("fields", "svisc", "th");

    bs.has_slope = std::abs(bs.alpha) > 0.;
    bs.has_N2 = std::abs(bs.n2) > 0.;

    if (bs.has_slope || bs.has_N2)
        master.print_message("Slope-enabled thermodynamics is activated\n");

    swbaroclinic = inputin.get_item<bool>("thermo", "swbaroclinic", "", false);
    if (swbaroclinic)
        dbdy_ls = inputin.get_item<TF>("thermo", "dbdy_ls", "");
}

template<typename TF>
Thermo_buoy<TF>::~Thermo_buoy()
{
}

template<typename TF>
void Thermo_buoy<TF>::create(
        Input& inputin, Netcdf_handle& input_nc, Stats<TF>& stats,
        Column<TF>& column, Cross<TF>& cross, Dump<TF>& dump, Timeloop<TF>& timeloop)
{
    stats.add_tendency(*fields.mt.at("w"), "zh", tend_name, tend_longname);
}

#ifndef USECUDA
template<typename TF>
void Thermo_buoy<TF>::exec(const double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();
    if (grid.get_spatial_order() == Grid_order::Second)
    {
        if (bs.has_slope || bs.has_N2)
        {
            calc_buoyancy_tend_u_2nd(fields.mt.at("u")->fld.data(), fields.sp.at("b")->fld.data(), bs.alpha, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
            calc_buoyancy_tend_w_2nd(fields.mt.at("w")->fld.data(), fields.sp.at("b")->fld.data(), bs.alpha, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
            calc_buoyancy_tend_b_2nd(fields.st.at("b")->fld.data(), fields.mp.at("u")->fld.data(), fields.mp.at("w")->fld.data(), bs.alpha, bs.n2, gd.utrans, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        }
        else
        {
            calc_buoyancy_tend_2nd(fields.mt.at("w")->fld.data(), fields.sp.at("b")->fld.data(), gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        }

        if (swbaroclinic)
            calc_baroclinic_2nd(
                    fields.st.at("b")->fld.data(), fields.mp.at("v")->fld.data(),
                    dbdy_ls,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
    }
    else if (grid.get_spatial_order() == Grid_order::Fourth)
    {
        if (bs.has_slope || bs.has_N2)
        {
            calc_buoyancy_tend_u_4th(fields.mt.at("u")->fld.data(), fields.sp.at("b")->fld.data(), bs.alpha, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
            calc_buoyancy_tend_w_4th(fields.mt.at("w")->fld.data(), fields.sp.at("b")->fld.data(), bs.alpha, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
            calc_buoyancy_tend_b_4th(fields.st.at("b")->fld.data(), fields.mp.at("u")->fld.data(), fields.mp.at("w")->fld.data(),  bs.alpha, bs.n2, gd.utrans, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        }
        else
        {
            calc_buoyancy_tend_4th(fields.mt.at("w")->fld.data(), fields.sp.at("b")->fld.data(), gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        }

        if (swbaroclinic)
            calc_baroclinic_4th(
                    fields.st.at("b")->fld.data(), fields.mp.at("v")->fld.data(),
                    dbdy_ls,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
    }
    stats.calc_tend(*fields.mt.at("w"), tend_name);

}
#endif

template<typename TF>
unsigned long Thermo_buoy<TF>::get_time_limit(unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}

template<typename TF>
void Thermo_buoy<TF>::get_thermo_field(
        Field3d<TF>& b, const std::string& name, const bool cyclic, const bool is_stat)
{
    auto& gd = grid.get_grid_data();

    if (name == "b")
        calc_buoyancy(b.fld.data(), fields.sp.at("b")->fld.data(),gd.ncells);
    else if (name == "N2")
        calc_N2(b.fld.data(), fields.sp.at("b")->fld.data(), bs.n2, gd.dzi.data(),gd.istart, gd.iend,
                gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells, gd.kcells);
    else
        throw 1;

    // Note: calc_buoyancy already handles the lateral ghost cells
}

template<typename TF>
void Thermo_buoy<TF>::get_prog_vars(std::vector<std::string>& list)
{
    list.push_back("b");
}

template<typename TF>
void Thermo_buoy<TF>::get_buoyancy_fluxbot(std::vector<TF>& bfluxbot, bool is_stat)
{
    auto& gd = grid.get_grid_data();

    calc_buoyancy_fluxbot(
        bfluxbot.data(),
        fields.sp.at("b")->flux_bot.data(),
        gd.icells, gd.jcells);
}

template<typename TF>
void Thermo_buoy<TF>::get_buoyancy_surf(std::vector<TF>& b, std::vector<TF>& bbot, bool is_stat)
{
    auto& gd = grid.get_grid_data();

    calc_buoyancy_bot(
        b.data(), bbot.data(),
        fields.sp.at("b")->fld.data(),
        fields.sp.at("b")->fld_bot.data(),
        gd.icells, gd.jcells,
        gd.ijcells, gd.kstart);

    //calc_buoyancy_fluxbot(b.flux_bot.data(), fields.sp.at("b")->flux_bot.data(), gd.icells, gd.jcells);
}

template<typename TF>
TF Thermo_buoy<TF>::get_buoyancy_diffusivity()
{
    return fields.sp.at("b")->visc;
}

template<typename TF>
int Thermo_buoy<TF>::get_bl_depth()
{
    // Use the buoyancy gradient to find the BL depth
    auto& gd = grid.get_grid_data();
    return calc_zi(fields.sp.at("b")->fld_mean.data(), gd.kstart, gd.kend, 1);
}

template<typename TF>
bool Thermo_buoy<TF>::check_field_exists(std::string name)
{
    if (name == "b")
        return true;
    else
        return false;
}

template class Thermo_buoy<double>;
template class Thermo_buoy<float>;
