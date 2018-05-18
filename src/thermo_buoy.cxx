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
#include <cmath>
#include <iostream>
#include "grid.h"
#include "fields.h"
#include "thermo_buoy.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "master.h"

using Finite_difference::O2::interp2;
using Finite_difference::O4::interp4;


namespace
{
    template<typename TF>
    void calc_buoyancy(TF* restrict b, TF* restrict bin, const int ncells)
    {
        for (int n=0; n<ncells; ++n)
            b[n] = bin[n];
    }

    template<typename TF>
    void calc_buoyancy_bot(TF* restrict b  , TF* restrict bbot,
                                        TF* restrict bin, TF* restrict binbot, const int icells, const int jcells, const int ijcells, const int kstart)
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
    void calc_buoyancy_fluxbot(TF* restrict bfluxbot, TF* restrict binfluxbot, const int icells, const int jcells)
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
    void calc_buoyancy_tend_2nd(TF* restrict wt, TF* restrict b, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
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
    void calc_buoyancy_tend_u_2nd(TF* restrict ut, TF* restrict b, const TF alpha, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
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
    void calc_buoyancy_tend_w_2nd(TF* restrict wt, TF* restrict b, const TF alpha, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
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
    void calc_buoyancy_tend_b_2nd(TF* restrict bt, TF* restrict u, TF* restrict w, const TF alpha, const TF n2, const TF utrans, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
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
    void calc_buoyancy_tend_4th(TF* restrict wt, TF* restrict b, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int jj, const int ijcells)
    {
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk1;
                    wt[ijk] += interp4(b[ijk-kk2], b[ijk-kk1], b[ijk], b[ijk+kk1]);
                }
    }

    template<typename TF>
    void calc_buoyancy_tend_u_4th(TF* restrict ut, TF* restrict b, const TF alpha, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
    {

        const TF sinalpha = std::sin(alpha);

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    ut[ijk] += sinalpha * interp4(b[ijk-2], b[ijk-1], b[ijk], b[ijk+1]);
                }
    }

    template<typename TF>
    void calc_buoyancy_tend_w_4th(TF* restrict wt, TF* restrict b, const TF alpha, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
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
                    wt[ijk] += cosalpha * interp4(b[ijk-kk2], b[ijk-kk1], b[ijk], b[ijk+kk1]);
                }
    }
    template<typename TF>
    void calc_buoyancy_tend_b_4th(TF* restrict bt, TF* restrict u, TF* restrict w, const TF alpha, const TF n2, const TF utrans, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
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
                    bt[ijk] -= n2 * ( sinalpha * (interp4(u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2]) + utrans)
                                    + cosalpha *  interp4(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]) );
                }
    }
}

template<typename TF>
Thermo_buoy<TF>::Thermo_buoy(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
Thermo<TF>(masterin, gridin, fieldsin, inputin)
{
    swthermo = "buoy";

    fields.init_prognostic_field("b", "Buoyancy", "m s-2");

    bs.alpha = inputin.get_item<TF>("thermo", "alpha", "", 0.);
    bs.n2 = inputin.get_item<TF>("thermo", "N2"   , "", 0.);
    fields.sp.at("b")->visc = inputin.get_item<TF>("fields", "svisc", "th");

	bs.has_slope = std::abs(bs.alpha) > 0.;
	bs.has_N2 = std::abs(bs.n2) > 0.;

    if (bs.has_slope || bs.has_N2)
        master.print_message("Slope-enabled thermodynamics is activated\n");

}

template<typename TF>
Thermo_buoy<TF>::~Thermo_buoy()
{
}

#ifndef USECUDA
template<typename TF>
void Thermo_buoy<TF>::exec()
{
    auto& gd = grid.get_grid_data();
    if (grid.swspatialorder == "2")
    {
	    if (bs.has_slope || bs.has_N2)
	    {
            calc_buoyancy_tend_u_2nd(fields.mt.at("u")->fld.data(), fields.sp.at("b")->fld.data(), bs.alpha, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
            calc_buoyancy_tend_w_2nd(fields.mt.at("w")->fld.data(), fields.sp.at("b")->fld.data(), bs.alpha, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
            calc_buoyancy_tend_b_2nd(fields.st.at("b")->fld.data(), fields.mp.at("u")->fld.data(), fields.mp.at("w")->fld.data(), bs.alpha, bs.n2, grid.utrans, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        }
        else
        {
	        calc_buoyancy_tend_2nd(fields.mt.at("w")->fld.data(), fields.sp.at("b")->fld.data(), gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        }
    }
    else if (grid.swspatialorder == "4")
    {
	    if (bs.has_slope || bs.has_N2)
	    {
		    calc_buoyancy_tend_u_4th(fields.mt.at("u")->fld.data(), fields.sp.at("b")->fld.data(), bs.alpha, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
            calc_buoyancy_tend_w_4th(fields.mt.at("w")->fld.data(), fields.sp.at("b")->fld.data(), bs.alpha, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
            calc_buoyancy_tend_b_4th(fields.st.at("b")->fld.data(), fields.mp.at("u")->fld.data(), fields.mp.at("w")->fld.data(),  bs.alpha, bs.n2, grid.utrans, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
	    }
	    else
	    {
		    calc_buoyancy_tend_4th(fields.mt.at("w")->fld.data(), fields.sp.at("b")->fld.data(), gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
	    }
    }
}
#endif

template<typename TF>
unsigned long Thermo_buoy<TF>::get_time_limit(unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}

template<typename TF>
void Thermo_buoy<TF>::get_thermo_field(Field3d<TF>& b, std::string name, bool cyclic, bool is_stat)
{
    auto& gd = grid.get_grid_data();
    calc_buoyancy(b.fld.data(), fields.sp.at("b")->fld.data(), gd.ncells);

    // Note: calc_buoyancy already handles the lateral ghost cells
}

template<typename TF>
void Thermo_buoy<TF>::get_prog_vars(std::vector<std::string>& list)
{
    list.push_back("b");
}

template<typename TF>
void Thermo_buoy<TF>::get_buoyancy_fluxbot(Field3d<TF>& b, bool is_stat)
{
    auto& gd = grid.get_grid_data();
    calc_buoyancy_fluxbot(b.flux_bot.data(), fields.sp.at("b")->flux_bot.data(), gd.icells, gd.jcells);
}

template<typename TF>
void Thermo_buoy<TF>::get_buoyancy_surf(Field3d<TF>& b, bool is_stat)
{
    auto& gd = grid.get_grid_data();
    calc_buoyancy_bot(b.fld.data()         , b.fld_bot.data(),
                      fields.sp.at("b")->fld.data(), fields.sp.at("b")->fld_bot.data(), gd.icells, gd.jcells, gd.ijcells, gd.kstart);
    calc_buoyancy_fluxbot(b.flux_bot.data(), fields.sp.at("b")->flux_bot.data(), gd.icells, gd.jcells);
}

template<typename TF>
TF Thermo_buoy<TF>::get_buoyancy_diffusivity()
{
    return fields.sp["b"]->visc;
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
