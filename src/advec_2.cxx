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

#include <algorithm>
#include <cmath>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "advec_2.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"

namespace
{
    using namespace Finite_difference::O2;
}

template<typename TF>
Advec_2<TF>::Advec_2(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Advec<TF>(masterin, gridin, fieldsin, inputin)
{
}

template<typename TF>
Advec_2<TF>::~Advec_2() {}

namespace
{
    template<typename TF>
    TF calc_cfl(const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzi, const TF dx, const TF dy,
            const TF dt, Master& master,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii = 1;

        const TF dxi = 1./dx;
        const TF dyi = 1./dy;

        TF cfl = 0;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    cfl = std::max(cfl, std::abs(interp2(u[ijk], u[ijk+ii]))*dxi + std::abs(interp2(v[ijk], v[ijk+jj]))*dyi + std::abs(interp2(w[ijk], w[ijk+kk]))*dzi[k]);
                }

        master.max(&cfl, 1);

        cfl = cfl*dt;

        return cfl;
    }

    template<typename TF>
    void advec_u(TF* const restrict ut,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzi, const TF dx, const TF dy,
            const TF* const restrict rhoref, const TF* const restrict rhorefh,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii = 1;

        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    ut[ijk] +=
                             - ( interp2(u[ijk   ], u[ijk+ii]) * interp2(u[ijk   ], u[ijk+ii])
                               - interp2(u[ijk-ii], u[ijk   ]) * interp2(u[ijk-ii], u[ijk   ]) ) * dxi

                             - ( interp2(v[ijk-ii+jj], v[ijk+jj]) * interp2(u[ijk   ], u[ijk+jj])
                               - interp2(v[ijk-ii   ], v[ijk   ]) * interp2(u[ijk-jj], u[ijk   ]) ) * dyi

                             - ( rhorefh[k+1] * interp2(w[ijk-ii+kk], w[ijk+kk]) * interp2(u[ijk   ], u[ijk+kk])
                               - rhorefh[k  ] * interp2(w[ijk-ii   ], w[ijk   ]) * interp2(u[ijk-kk], u[ijk   ]) ) / rhoref[k] * dzi[k];
                }
    }

    template<typename TF>
    void advec_v(TF* const restrict vt,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzi, const TF dx, const TF dy,
            const TF* const restrict rhoref, const TF* const restrict rhorefh,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii = 1;

        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    vt[ijk] +=
                             - ( interp2(u[ijk+ii-jj], u[ijk+ii]) * interp2(v[ijk   ], v[ijk+ii])
                               - interp2(u[ijk   -jj], u[ijk   ]) * interp2(v[ijk-ii], v[ijk   ]) ) * dxi

                             - ( interp2(v[ijk   ], v[ijk+jj]) * interp2(v[ijk   ], v[ijk+jj])
                               - interp2(v[ijk-jj], v[ijk   ]) * interp2(v[ijk-jj], v[ijk   ]) ) * dyi

                             - ( rhorefh[k+1] * interp2(w[ijk-jj+kk], w[ijk+kk]) * interp2(v[ijk   ], v[ijk+kk])
                               - rhorefh[k  ] * interp2(w[ijk-jj   ], w[ijk   ]) * interp2(v[ijk-kk], v[ijk   ]) ) / rhoref[k] * dzi[k];
                }
    }

    template<typename TF>
    void advec_w(TF* const restrict wt,
            const TF* const restrict u, const TF* const restrict v, TF* const restrict w,
            const TF* const restrict dzhi, const TF dx, const TF dy,
            const TF* const restrict rhoref, const TF* const restrict rhorefh,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii = 1;

        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    wt[ijk] +=
                             - ( interp2(u[ijk+ii-kk], u[ijk+ii]) * interp2(w[ijk   ], w[ijk+ii])
                               - interp2(u[ijk   -kk], u[ijk   ]) * interp2(w[ijk-ii], w[ijk   ]) ) * dxi

                             - ( interp2(v[ijk+jj-kk], v[ijk+jj]) * interp2(w[ijk   ], w[ijk+jj])
                               - interp2(v[ijk   -kk], v[ijk   ]) * interp2(w[ijk-jj], w[ijk   ]) ) * dyi

                             - ( rhoref[k  ] * interp2(w[ijk   ], w[ijk+kk]) * interp2(w[ijk   ], w[ijk+kk])
                               - rhoref[k-1] * interp2(w[ijk-kk], w[ijk   ]) * interp2(w[ijk-kk], w[ijk   ]) ) / rhorefh[k] * dzhi[k];
                }
    }

    template<typename TF>
    void advec_s(TF* const restrict st, const TF* const restrict s,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzi, const TF dx, const TF dy,
            const TF* const restrict rhoref, const TF* const restrict rhorefh,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii = 1;

        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] +=
                             - ( u[ijk+ii] * interp2(s[ijk   ], s[ijk+ii])
                               - u[ijk   ] * interp2(s[ijk-ii], s[ijk   ]) ) * dxi

                             - ( v[ijk+jj] * interp2(s[ijk   ], s[ijk+jj])
                               - v[ijk   ] * interp2(s[ijk-jj], s[ijk   ]) ) * dyi

                             - ( rhorefh[k+1] * w[ijk+kk] * interp2(s[ijk   ], s[ijk+kk])
                               - rhorefh[k  ] * w[ijk   ] * interp2(s[ijk-kk], s[ijk   ]) ) / rhoref[k] * dzi[k];
                }
    }

    template<typename TF>
    void advec_flux_u(
            TF* const restrict st, const TF* const restrict s, const TF* const restrict w,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii = 1;

        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] = interp2(w[ijk-ii], w[ijk]) * interp2(s[ijk-kk], s[ijk]);
                }
    }

    template<typename TF>
    void advec_flux_v(
            TF* const restrict st, const TF* const restrict s, const TF* const restrict w,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] = interp2(w[ijk-jj], w[ijk]) * interp2(s[ijk-kk], s[ijk]);
                }
    }

    template<typename TF>
    void advec_flux_s(
            TF* const restrict st, const TF* const restrict s, const TF* const restrict w,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] = w[ijk] * interp2(s[ijk-kk], s[ijk]);
                }
    }
}

template<typename TF>
void Advec_2<TF>::create(Stats<TF>& stats)
{
    stats.add_tendency(*fields.mt.at("u"), "z", tend_name, tend_longname);
    stats.add_tendency(*fields.mt.at("v"), "z", tend_name, tend_longname);
    stats.add_tendency(*fields.mt.at("w"), "zh", tend_name, tend_longname);
    for (auto it : fields.st)
        stats.add_tendency(*it.second, "z", tend_name, tend_longname);
}

#ifndef USECUDA
template<typename TF>
double Advec_2<TF>::get_cfl(double dt)
{
    auto& gd = grid.get_grid_data();
    TF cfl = calc_cfl<TF>(fields.mp.at("u")->fld.data(),fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                          gd.dzi.data(), gd.dx, gd.dy,
                          dt, master,
                          gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                          gd.icells, gd.ijcells);

    // CFL is kept in double precision for time stepping accuracy
    return static_cast<double>(cfl);
}

template<typename TF>
unsigned long Advec_2<TF>::get_time_limit(unsigned long idt, double dt)
{
    // Calculate cfl and prevent zero divisons.
    auto& gd = grid.get_grid_data();

    double cfl = calc_cfl<TF>(fields.mp.at("u")->fld.data(),fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dx, gd.dy,
            dt, master,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    cfl = std::max(cflmin, cfl);
    return idt * cflmax / cfl;
}


template<typename TF>
void Advec_2<TF>::exec(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();
    advec_u(fields.mt.at("u")->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dx, gd.dy,
            fields.rhoref.data(), fields.rhorefh.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    advec_v(fields.mt.at("v")->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dx, gd.dy,
            fields.rhoref.data(), fields.rhorefh.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    advec_w(fields.mt.at("w")->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            gd.dzhi.data(), gd.dx, gd.dy,
            fields.rhoref.data(), fields.rhorefh.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    for (auto& it : fields.st)
        advec_s(it.second->fld.data(), fields.sp.at(it.first)->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                gd.dzi.data(), gd.dx, gd.dy,
                fields.rhoref.data(), fields.rhorefh.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);
    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);
}
#endif

template<typename TF>
void Advec_2<TF>::get_advec_flux(Field3d<TF>& advec_flux, const Field3d<TF>& fld)
{
    auto& gd = grid.get_grid_data();

    if (fld.loc == gd.uloc)
    {
        advec_flux_u(
                advec_flux.fld.data(), fld.fld.data(), fields.mp.at("w")->fld.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }
    else if (fld.loc == gd.vloc)
    {
        advec_flux_v(
                advec_flux.fld.data(), fld.fld.data(), fields.mp.at("w")->fld.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }
    else if (fld.loc == gd.sloc)
    {
        advec_flux_s(
                advec_flux.fld.data(), fld.fld.data(), fields.mp.at("w")->fld.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }
    else
        throw std::runtime_error("Advec_2 cannot deliver flux field at that location");
}

template class Advec_2<double>;
template class Advec_2<float>;
