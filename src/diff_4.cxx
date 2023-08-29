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
#include <algorithm>

#include "grid.h"
#include "fields.h"
#include "master.h"
#include "stats.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"

#include "diff_4.h"

namespace
{
    using namespace Finite_difference::O4;

    template<typename TF, bool dim3>
    void diff_c(TF* restrict at, const TF* restrict a, const TF visc,
                const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                const int jj, const int kk, const TF dx, const TF dy, const TF* restrict dzi4, const TF* restrict dzhi4)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int jj3 = 3*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;
        const int kk3 = 3*kk;

        const TF dxidxi = 1./(dx*dx);
        const TF dyidyi = 1./(dy*dy);

        // bottom boundary
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj1 + kstart*kk1;
                at[ijk] += visc * (cdg3<TF>*a[ijk-ii3] + cdg2<TF>*a[ijk-ii2] + cdg1<TF>*a[ijk-ii1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+ii1] + cdg2<TF>*a[ijk+ii2] + cdg3<TF>*a[ijk+ii3])*dxidxi;
                if (dim3)
                    at[ijk] += visc * (cdg3<TF>*a[ijk-jj3] + cdg2<TF>*a[ijk-jj2] + cdg1<TF>*a[ijk-jj1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+jj1] + cdg2<TF>*a[ijk+jj2] + cdg3<TF>*a[ijk+jj3])*dyidyi;
                at[ijk] += visc * ( cg0<TF>*(bg0<TF>*a[ijk-kk2] + bg1<TF>*a[ijk-kk1] + bg2<TF>*a[ijk    ] + bg3<TF>*a[ijk+kk1]) * dzhi4[kstart-1]
                                  + cg1<TF>*(cg0<TF>*a[ijk-kk2] + cg1<TF>*a[ijk-kk1] + cg2<TF>*a[ijk    ] + cg3<TF>*a[ijk+kk1]) * dzhi4[kstart  ]
                                  + cg2<TF>*(cg0<TF>*a[ijk-kk1] + cg1<TF>*a[ijk    ] + cg2<TF>*a[ijk+kk1] + cg3<TF>*a[ijk+kk2]) * dzhi4[kstart+1]
                                  + cg3<TF>*(cg0<TF>*a[ijk    ] + cg1<TF>*a[ijk+kk1] + cg2<TF>*a[ijk+kk2] + cg3<TF>*a[ijk+kk3]) * dzhi4[kstart+2] )
                                * dzi4[kstart];
            }

        for (int k=kstart+1; k<kend-1; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    at[ijk] += visc * (cdg3<TF>*a[ijk-ii3] + cdg2<TF>*a[ijk-ii2] + cdg1<TF>*a[ijk-ii1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+ii1] + cdg2<TF>*a[ijk+ii2] + cdg3<TF>*a[ijk+ii3])*dxidxi;
                    if (dim3)
                        at[ijk] += visc * (cdg3<TF>*a[ijk-jj3] + cdg2<TF>*a[ijk-jj2] + cdg1<TF>*a[ijk-jj1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+jj1] + cdg2<TF>*a[ijk+jj2] + cdg3<TF>*a[ijk+jj3])*dyidyi;
                    at[ijk] += visc * ( cg0<TF>*(cg0<TF>*a[ijk-kk3] + cg1<TF>*a[ijk-kk2] + cg2<TF>*a[ijk-kk1] + cg3<TF>*a[ijk    ]) * dzhi4[k-1]
                                      + cg1<TF>*(cg0<TF>*a[ijk-kk2] + cg1<TF>*a[ijk-kk1] + cg2<TF>*a[ijk    ] + cg3<TF>*a[ijk+kk1]) * dzhi4[k  ]
                                      + cg2<TF>*(cg0<TF>*a[ijk-kk1] + cg1<TF>*a[ijk    ] + cg2<TF>*a[ijk+kk1] + cg3<TF>*a[ijk+kk2]) * dzhi4[k+1]
                                      + cg3<TF>*(cg0<TF>*a[ijk    ] + cg1<TF>*a[ijk+kk1] + cg2<TF>*a[ijk+kk2] + cg3<TF>*a[ijk+kk3]) * dzhi4[k+2] )
                                    * dzi4[k];
                }

        // top boundary
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj1 + (kend-1)*kk1;
                at[ijk] += visc * (cdg3<TF>*a[ijk-ii3] + cdg2<TF>*a[ijk-ii2] + cdg1<TF>*a[ijk-ii1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+ii1] + cdg2<TF>*a[ijk+ii2] + cdg3<TF>*a[ijk+ii3])*dxidxi;
                if (dim3)
                    at[ijk] += visc * (cdg3<TF>*a[ijk-jj3] + cdg2<TF>*a[ijk-jj2] + cdg1<TF>*a[ijk-jj1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+jj1] + cdg2<TF>*a[ijk+jj2] + cdg3<TF>*a[ijk+jj3])*dyidyi;
                at[ijk] += visc * ( cg0<TF>*(cg0<TF>*a[ijk-kk3] + cg1<TF>*a[ijk-kk2] + cg2<TF>*a[ijk-kk1] + cg3<TF>*a[ijk    ]) * dzhi4[kend-2]
                                  + cg1<TF>*(cg0<TF>*a[ijk-kk2] + cg1<TF>*a[ijk-kk1] + cg2<TF>*a[ijk    ] + cg3<TF>*a[ijk+kk1]) * dzhi4[kend-1]
                                  + cg2<TF>*(cg0<TF>*a[ijk-kk1] + cg1<TF>*a[ijk    ] + cg2<TF>*a[ijk+kk1] + cg3<TF>*a[ijk+kk2]) * dzhi4[kend  ]
                                  + cg3<TF>*(tg0<TF>*a[ijk-kk1] + tg1<TF>*a[ijk    ] + tg2<TF>*a[ijk+kk1] + tg3<TF>*a[ijk+kk2]) * dzhi4[kend+1] )
                                * dzi4[kend-1];
            }
    }

    template<typename TF, bool dim3>
    void diff_w(TF* restrict at, const TF* restrict a, const TF visc,
                const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                const int jj, const int kk, const TF dx, const TF dy, const TF* restrict dzi4, const TF* restrict dzhi4)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int jj3 = 3*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;
        const int kk3 = 3*kk;

        const TF dxidxi = 1/(dx*dx);
        const TF dyidyi = 1/(dy*dy);

        // bottom boundary
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj1 + (kstart+1)*kk1;
                at[ijk] += visc * (cdg3<TF>*a[ijk-ii3] + cdg2<TF>*a[ijk-ii2] + cdg1<TF>*a[ijk-ii1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+ii1] + cdg2<TF>*a[ijk+ii2] + cdg3<TF>*a[ijk+ii3])*dxidxi;
                if (dim3)
                    at[ijk] += visc * (cdg3<TF>*a[ijk-jj3] + cdg2<TF>*a[ijk-jj2] + cdg1<TF>*a[ijk-jj1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+jj1] + cdg2<TF>*a[ijk+jj2] + cdg3<TF>*a[ijk+jj3])*dyidyi;
                at[ijk] += visc * ( cg0<TF>*(bg0<TF>*a[ijk-kk2] + bg1<TF>*a[ijk-kk1] + bg2<TF>*a[ijk    ] + bg3<TF>*a[ijk+kk1]) * dzi4[kstart-1]
                                  + cg1<TF>*(cg0<TF>*a[ijk-kk2] + cg1<TF>*a[ijk-kk1] + cg2<TF>*a[ijk    ] + cg3<TF>*a[ijk+kk1]) * dzi4[kstart  ]
                                  + cg2<TF>*(cg0<TF>*a[ijk-kk1] + cg1<TF>*a[ijk    ] + cg2<TF>*a[ijk+kk1] + cg3<TF>*a[ijk+kk2]) * dzi4[kstart+1]
                                  + cg3<TF>*(cg0<TF>*a[ijk    ] + cg1<TF>*a[ijk+kk1] + cg2<TF>*a[ijk+kk2] + cg3<TF>*a[ijk+kk3]) * dzi4[kstart+2] )
                                * dzhi4[kstart+1];
            }

        for (int k=kstart+2; k<kend-1; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    at[ijk] += visc * (cdg3<TF>*a[ijk-ii3] + cdg2<TF>*a[ijk-ii2] + cdg1<TF>*a[ijk-ii1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+ii1] + cdg2<TF>*a[ijk+ii2] + cdg3<TF>*a[ijk+ii3])*dxidxi;
                    if (dim3)
                        at[ijk] += visc * (cdg3<TF>*a[ijk-jj3] + cdg2<TF>*a[ijk-jj2] + cdg1<TF>*a[ijk-jj1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+jj1] + cdg2<TF>*a[ijk+jj2] + cdg3<TF>*a[ijk+jj3])*dyidyi;
                    at[ijk] += visc * ( cg0<TF>*(cg0<TF>*a[ijk-kk3] + cg1<TF>*a[ijk-kk2] + cg2<TF>*a[ijk-kk1] + cg3<TF>*a[ijk    ]) * dzi4[k-2]
                                      + cg1<TF>*(cg0<TF>*a[ijk-kk2] + cg1<TF>*a[ijk-kk1] + cg2<TF>*a[ijk    ] + cg3<TF>*a[ijk+kk1]) * dzi4[k-1]
                                      + cg2<TF>*(cg0<TF>*a[ijk-kk1] + cg1<TF>*a[ijk    ] + cg2<TF>*a[ijk+kk1] + cg3<TF>*a[ijk+kk2]) * dzi4[k  ]
                                      + cg3<TF>*(cg0<TF>*a[ijk    ] + cg1<TF>*a[ijk+kk1] + cg2<TF>*a[ijk+kk2] + cg3<TF>*a[ijk+kk3]) * dzi4[k+1] )
                                    * dzhi4[k];
                }

        // top boundary
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj1 + (kend-1)*kk1;
                at[ijk] += visc * (cdg3<TF>*a[ijk-ii3] + cdg2<TF>*a[ijk-ii2] + cdg1<TF>*a[ijk-ii1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+ii1] + cdg2<TF>*a[ijk+ii2] + cdg3<TF>*a[ijk+ii3])*dxidxi;
                if (dim3)
                    at[ijk] += visc * (cdg3<TF>*a[ijk-jj3] + cdg2<TF>*a[ijk-jj2] + cdg1<TF>*a[ijk-jj1] + cdg0<TF>*a[ijk] + cdg1<TF>*a[ijk+jj1] + cdg2<TF>*a[ijk+jj2] + cdg3<TF>*a[ijk+jj3])*dyidyi;
                at[ijk] += visc * ( cg0<TF>*(cg0<TF>*a[ijk-kk3] + cg1<TF>*a[ijk-kk2] + cg2<TF>*a[ijk-kk1] + cg3<TF>*a[ijk    ]) * dzi4[kend-3]
                                  + cg1<TF>*(cg0<TF>*a[ijk-kk2] + cg1<TF>*a[ijk-kk1] + cg2<TF>*a[ijk    ] + cg3<TF>*a[ijk+kk1]) * dzi4[kend-2]
                                  + cg2<TF>*(cg0<TF>*a[ijk-kk1] + cg1<TF>*a[ijk    ] + cg2<TF>*a[ijk+kk1] + cg3<TF>*a[ijk+kk2]) * dzi4[kend-1]
                                  + cg3<TF>*(tg0<TF>*a[ijk-kk1] + tg1<TF>*a[ijk    ] + tg2<TF>*a[ijk+kk1] + tg3<TF>*a[ijk+kk2]) * dzi4[kend  ] )
                                * dzhi4[kend-1];
            }
    }

    template<typename TF>
    void calc_diff_flux(
            TF* const restrict out, const TF* const restrict data, const TF visc, const TF* const dzhi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj  = 1*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        #pragma omp parallel for
        for (int k=kstart; k<kend+1; ++k)
        {
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk1;
                        out[ijk] = - visc*(cg0<TF>*data[ijk-kk2] + cg1<TF>*data[ijk-kk1] + cg2<TF>*data[ijk] + cg3<TF>*data[ijk+kk1])*dzhi[k];
                    }
        }
    }
}

template<typename TF>
Diff_4<TF>::Diff_4(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Boundary<TF>& boundaryin, Input& inputin) :
    Diff<TF>(masterin, gridin, fieldsin, boundaryin, inputin)
{
    dnmax = inputin.get_item<double>("diff", "dnmax", "", 0.4);
}

template <typename TF>
Diff_4<TF>::~Diff_4()
{
}

template <typename TF>
Diffusion_type Diff_4<TF>::get_switch() const
{
    return swdiff;
}

template<typename TF>
unsigned long Diff_4<TF>::get_time_limit(const unsigned long idt, const double dt)
{
    return idt * dnmax / (dt * dnmul);
}

template<typename TF>
double Diff_4<TF>::get_dn(const double dt)
{
    return dnmul*dt;
}

template<typename TF>
void Diff_4<TF>::create(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // Get the maximum viscosity
    TF viscmax = fields.visc;
    for (auto& it : fields.sp)
        viscmax = std::max(it.second->visc, viscmax);

    // Calculate time step multiplier for diffusion number
    dnmul = 0;
    for (int k=gd.kstart; k<gd.kend; ++k)
        dnmul = std::max(dnmul, std::abs(viscmax * (1./(gd.dx*gd.dx) + 1./(gd.dy*gd.dy) + 1./(gd.dz[k]*gd.dz[k]))));

    stats.add_tendency(*fields.mt.at("u"), "z", tend_name, tend_longname);
    stats.add_tendency(*fields.mt.at("v"), "z", tend_name, tend_longname);
    stats.add_tendency(*fields.mt.at("w"), "zh", tend_name, tend_longname);
    for (auto it : fields.st)
        stats.add_tendency(*it.second, "z", tend_name, tend_longname);
}

#ifndef USECUDA
template<typename TF>
void Diff_4<TF>::exec(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // 2D simulation
    if (gd.jtot == 1)
    {
        diff_c<TF,0>(fields.mt.at("u")->fld.data(), fields.mp.at("u")->fld.data(), fields.visc,
                     gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells,
                     gd.dx, gd.dy, gd.dzi4.data(), gd.dzhi4.data());

        diff_c<TF,0>(fields.mt.at("v")->fld.data(), fields.mp.at("v")->fld.data(), fields.visc,
                     gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells,
                     gd.dx, gd.dy, gd.dzi4.data(), gd.dzhi4.data());

        diff_w<TF,0>(fields.mt.at("w")->fld.data(), fields.mp.at("w")->fld.data(), fields.visc,
                     gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells,
                     gd.dx, gd.dy, gd.dzi4.data(), gd.dzhi4.data());

        for (auto& it : fields.st)
            diff_c<TF,0>(it.second->fld.data(), fields.sp.at(it.first)->fld.data(), fields.sp.at(it.first)->visc,
                         gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells,
                         gd.dx, gd.dy, gd.dzi4.data(), gd.dzhi4.data());
    }
    else
    {
        diff_c<TF,1>(fields.mt.at("u")->fld.data(), fields.mp.at("u")->fld.data(), fields.visc,
                     gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells,
                     gd.dx, gd.dy, gd.dzi4.data(), gd.dzhi4.data());

        diff_c<TF,1>(fields.mt.at("v")->fld.data(), fields.mp.at("v")->fld.data(), fields.visc,
                     gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells,
                     gd.dx, gd.dy, gd.dzi4.data(), gd.dzhi4.data());

        diff_w<TF,1>(fields.mt.at("w")->fld.data(), fields.mp.at("w")->fld.data(), fields.visc,
                     gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells,
                     gd.dx, gd.dy, gd.dzi4.data(), gd.dzhi4.data());

        for (auto& it : fields.st)
            diff_c<TF,1>(it.second->fld.data(), fields.sp.at(it.first)->fld.data(), fields.sp.at(it.first)->visc,
                         gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells,
                         gd.dx, gd.dy, gd.dzi4.data(), gd.dzhi4.data());
    }

    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);
    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);
}
#endif


template<typename TF>
void Diff_4<TF>::diff_flux(Field3d<TF>& restrict out, const Field3d<TF>& in)
{
    auto& gd = grid.get_grid_data();
    calc_diff_flux(
            out.fld.data(), in.fld.data(), in.visc, gd.dzhi4.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);
}


#ifdef FLOAT_SINGLE
template class Diff_4<float>;
#else
template class Diff_4<double>;
#endif
