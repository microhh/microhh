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
#include <algorithm>

#include "grid.h"
#include "fields.h"
#include "master.h"
#include "defines.h"
#include "model.h"
#include "constants.h"

#include "diff_2.h"

namespace
{
    template<typename TF>
    void diff_c(TF* restrict at, TF* restrict a, TF visc,
                const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                const int jj, const int kk, const TF dx, const TF dy, const TF* const dzi, const TF* const dzhi)
    {
        const int ii = 1;
        const double dxidxi = 1/(dx*dx);
        const double dyidyi = 1/(dy*dy);
    
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    at[ijk] += visc * (
                            + ( (a[ijk+ii] - a[ijk   ]) 
                              - (a[ijk   ] - a[ijk-ii]) ) * dxidxi 
                            + ( (a[ijk+jj] - a[ijk   ]) 
                              - (a[ijk   ] - a[ijk-jj]) ) * dyidyi
                            + ( (a[ijk+kk] - a[ijk   ]) * dzhi[k+1]
                              - (a[ijk   ] - a[ijk-kk]) * dzhi[k]   ) * dzi[k] );
                }
    }
    
    template<typename TF>
    void diff_w(TF* restrict wt, TF* restrict w, TF visc,
                const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                const int jj, const int kk, const TF dx, const TF dy, const TF* const dzi, const TF* const dzhi)
    {
        const int ii = 1;
        const double dxidxi = 1/(dx*dx);
        const double dyidyi = 1/(dy*dy);
    
        for (int k=kstart+1; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    wt[ijk] += visc * (
                            + ( (w[ijk+ii] - w[ijk   ]) 
                              - (w[ijk   ] - w[ijk-ii]) ) * dxidxi 
                            + ( (w[ijk+jj] - w[ijk   ]) 
                              - (w[ijk   ] - w[ijk-jj]) ) * dyidyi
                            + ( (w[ijk+kk] - w[ijk   ]) * dzi[k]
                              - (w[ijk   ] - w[ijk-kk]) * dzi[k-1] ) * dzhi[k] );
                }
    }
}

template<typename TF>
Diff_2<TF>::Diff_2(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Diff<TF>(masterin, gridin, fieldsin, inputin)
{
    dnmax = inputin.get_item<double>("diff", "dnmax", "", 0.4);
}

template <typename TF>
Diff_2<TF>::~Diff_2()
{
}

template <typename TF>
Diffusion_type Diff_2<TF>::get_switch()
{
    return swdiff;
}

template<typename TF>
unsigned long Diff_2<TF>::get_time_limit(const unsigned long idt, const double dt)
{
    return idt * dnmax / (dt * dnmul);
}

template<typename TF>
double Diff_2<TF>::get_dn(const double dt)
{
    return dnmul*dt;
}

template<typename TF>
void Diff_2<TF>::set_values()
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
}

template<typename TF>
void Diff_2<TF>::exec()
{
    auto& gd = grid.get_grid_data();

    diff_c<TF>(fields.mt.at("u")->data.data(), fields.mp.at("u")->data.data(), fields.visc, 
               gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells,
               gd.dx, gd.dy, gd.dzi.data(), gd.dzhi.data());
    
    diff_c<TF>(fields.mt.at("v")->data.data(), fields.mp.at("v")->data.data(), fields.visc, 
               gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells,
               gd.dx, gd.dy, gd.dzi.data(), gd.dzhi.data());
    
    diff_w<TF>(fields.mt.at("w")->data.data(), fields.mp.at("w")->data.data(), fields.visc, 
               gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells,
               gd.dx, gd.dy, gd.dzi.data(), gd.dzhi.data());
    
    for (auto& it : fields.st)
        diff_c<TF>(it.second->data.data(), fields.sp.at(it.first)->data.data(), fields.sp.at(it.first)->visc, 
                   gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells,
                   gd.dx, gd.dy, gd.dzi.data(), gd.dzhi.data());
}

template class Diff_2<double>;
template class Diff_2<float>;