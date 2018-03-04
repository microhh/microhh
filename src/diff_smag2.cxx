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
#include "constants.h"
#include "monin_obukhov.h"

#include "diff_smag2.h"

namespace
{
    namespace most = Monin_obukhov;

    template <typename TF, bool resolved_wall> 
    void calc_strain2(TF* restrict strain2,
                      TF* restrict u, TF* restrict v, TF* restrict w,
                      TF* restrict ufluxbot, TF* restrict vfluxbot,
                      TF* restrict ustar, TF* restrict obuk,
                      const TF* restrict z, const TF* restrict dzi, const TF* restrict dzhi,
                      const TF dxi, const TF dyi,
                      const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                      const int jj, const int kk)
    {
        const int ii = 1;
        const int k_offset = resolved_wall ? 0 : 1;
    
        // If the wall isn't resolved, calculate du/dz and dv/dz at lowest grid height using MO
        if (!resolved_wall)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
    
                    strain2[ijk] = 2.*(
                                   // du/dx + du/dx
                                   + std::pow((u[ijk+ii]-u[ijk])*dxi, 2)
    
                                   // dv/dy + dv/dy
                                   + std::pow((v[ijk+jj]-v[ijk])*dyi, 2)
    
                                   // dw/dz + dw/dz
                                   + std::pow((w[ijk+kk]-w[ijk])*dzi[kstart], 2)
    
                                   // du/dy + dv/dx
                                   + 0.125*std::pow((u[ijk      ]-u[ijk   -jj])*dyi  + (v[ijk      ]-v[ijk-ii   ])*dxi, 2)
                                   + 0.125*std::pow((u[ijk+ii   ]-u[ijk+ii-jj])*dyi  + (v[ijk+ii   ]-v[ijk      ])*dxi, 2)
                                   + 0.125*std::pow((u[ijk   +jj]-u[ijk      ])*dyi  + (v[ijk   +jj]-v[ijk-ii+jj])*dxi, 2)
                                   + 0.125*std::pow((u[ijk+ii+jj]-u[ijk+ii   ])*dyi  + (v[ijk+ii+jj]-v[ijk   +jj])*dxi, 2)
    
                                   // du/dz
                                   + 0.5*std::pow(-0.5*(ufluxbot[ij]+ufluxbot[ij+ii])/(Constants::kappa*z[kstart]*ustar[ij])*most::phim(z[kstart]/obuk[ij]), 2)
    
                                   // dw/dx
                                   + 0.125*std::pow((w[ijk      ]-w[ijk-ii   ])*dxi, 2)
                                   + 0.125*std::pow((w[ijk+ii   ]-w[ijk      ])*dxi, 2)
                                   + 0.125*std::pow((w[ijk   +kk]-w[ijk-ii+kk])*dxi, 2)
                                   + 0.125*std::pow((w[ijk+ii+kk]-w[ijk   +kk])*dxi, 2)
    
                                   // dv/dz
                                   + 0.5*std::pow(-0.5*(vfluxbot[ij]+vfluxbot[ij+jj])/(Constants::kappa*z[kstart]*ustar[ij])*most::phim(z[kstart]/obuk[ij]), 2)
    
                                   // dw/dy
                                   + 0.125*std::pow((w[ijk      ]-w[ijk-jj   ])*dyi, 2)
                                   + 0.125*std::pow((w[ijk+jj   ]-w[ijk      ])*dyi, 2)
                                   + 0.125*std::pow((w[ijk   +kk]-w[ijk-jj+kk])*dyi, 2)
                                   + 0.125*std::pow((w[ijk+jj+kk]-w[ijk   +kk])*dyi, 2) );
    
                    // add a small number to avoid zero divisions
                    strain2[ijk] += Constants::dsmall;
                }
        }
    
        for (int k=kstart+k_offset; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    strain2[ijk] = 2.*(
                                   // du/dx + du/dx
                                   + std::pow((u[ijk+ii]-u[ijk])*dxi, 2)
    
                                   // dv/dy + dv/dy
                                   + std::pow((v[ijk+jj]-v[ijk])*dyi, 2)
    
                                   // dw/dz + dw/dz
                                   + std::pow((w[ijk+kk]-w[ijk])*dzi[k], 2)
    
                                   // du/dy + dv/dx
                                   + 0.125*std::pow((u[ijk      ]-u[ijk   -jj])*dyi  + (v[ijk      ]-v[ijk-ii   ])*dxi, 2)
                                   + 0.125*std::pow((u[ijk+ii   ]-u[ijk+ii-jj])*dyi  + (v[ijk+ii   ]-v[ijk      ])*dxi, 2)
                                   + 0.125*std::pow((u[ijk   +jj]-u[ijk      ])*dyi  + (v[ijk   +jj]-v[ijk-ii+jj])*dxi, 2)
                                   + 0.125*std::pow((u[ijk+ii+jj]-u[ijk+ii   ])*dyi  + (v[ijk+ii+jj]-v[ijk   +jj])*dxi, 2)
    
                                   // du/dz + dw/dx
                                   + 0.125*std::pow((u[ijk      ]-u[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-ii   ])*dxi, 2)
                                   + 0.125*std::pow((u[ijk+ii   ]-u[ijk+ii-kk])*dzhi[k  ] + (w[ijk+ii   ]-w[ijk      ])*dxi, 2)
                                   + 0.125*std::pow((u[ijk   +kk]-u[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-ii+kk])*dxi, 2)
                                   + 0.125*std::pow((u[ijk+ii+kk]-u[ijk+ii   ])*dzhi[k+1] + (w[ijk+ii+kk]-w[ijk   +kk])*dxi, 2)
    
                                   // dv/dz + dw/dy
                                   + 0.125*std::pow((v[ijk      ]-v[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-jj   ])*dyi, 2)
                                   + 0.125*std::pow((v[ijk+jj   ]-v[ijk+jj-kk])*dzhi[k  ] + (w[ijk+jj   ]-w[ijk      ])*dyi, 2)
                                   + 0.125*std::pow((v[ijk   +kk]-v[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-jj+kk])*dyi, 2)
                                   + 0.125*std::pow((v[ijk+jj+kk]-v[ijk+jj   ])*dzhi[k+1] + (w[ijk+jj+kk]-w[ijk   +kk])*dyi, 2) );
    
                           // Add a small number to avoid zero divisions.
                           strain2[ijk] += Constants::dsmall;
                }
    }

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
Diff_smag2<TF>::Diff_smag2(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Diff<TF>(masterin, gridin, fieldsin, inputin)
{
    dnmax = inputin.get_item<double>("diff", "dnmax", "", 0.4  );
    cs    = inputin.get_item<double>("diff", "cs"   , "", 0.23 );
    tPr   = inputin.get_item<double>("diff", "tPr"  , "", 1./3.);

    fields.init_diagnostic_field("evisc", "Eddy viscosity", "m2 s-1");
}

template <typename TF>
Diff_smag2<TF>::~Diff_smag2()
{
}

template <typename TF>
Diffusion_type Diff_smag2<TF>::get_switch()
{
    return swdiff;
}

template<typename TF>
unsigned long Diff_smag2<TF>::get_time_limit(const unsigned long idt, const double dt)
{
    return idt * dnmax / (dt * dnmul);
}

template<typename TF>
double Diff_smag2<TF>::get_dn(const double dt)
{
    return dnmul*dt;
}

template<typename TF>
void Diff_smag2<TF>::set_values()
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

#ifndef USECUDA
template<typename TF>
void Diff_smag2<TF>::exec()
{
    auto& gd = grid.get_grid_data();
}

template<typename TF>
void Diff_smag2<TF>::exec_viscosity()
{
    auto& gd = grid.get_grid_data();

    // Do a cast because the base boundary class does not have the MOST related variables.
    // Boundary_surface* boundaryptr = static_cast<Boundary_surface*>(model->boundary);

    // Calculate strain rate using MO for velocity gradients lowest level
    // if (model->boundary->get_switch() == "surface")
    //     calc_strain2<false>(fields->sd["evisc"]->data,
    //                         fields->u->data, fields->v->data, fields->w->data,
    //                         fields->u->datafluxbot, fields->v->datafluxbot,
    //                         boundaryptr->ustar, boundaryptr->obuk,
    //                         grid->z, grid->dzi, grid->dzhi);
    // Calculate strain rate using resolved boundaries
    // else
        calc_strain2<TF, true>(fields.sd["evisc"]->fld.data(),
                               fields.mp["u"]->fld.data(), fields.mp["v"]->fld.data(), fields.mp["w"]->fld.data(),
                               fields.mp["u"]->flux_bot.data(), fields.mp["v"]->flux_bot.data(),
                               nullptr, nullptr,
                               gd.z.data(), gd.dzi.data(), gd.dzhi.data(), 1./gd.dx, 1./gd.dy,
                               gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                               gd.icells, gd.ijcells);

    // start with retrieving the stability information
    // if (model->thermo->get_switch() == "0")
    // {
    //     // Calculate eddy viscosity using MO at lowest model level
    //     if (model->boundary->get_switch() == "surface")
    //         calc_evisc_neutral<false>(fields->sd["evisc"]->data,
    //                                   fields->u->data, fields->v->data, fields->w->data,
    //                                   fields->u->datafluxbot, fields->v->datafluxbot,
    //                                   grid->z, grid->dz, boundaryptr->z0m, fields->visc);
    //     // Calculate eddy viscosity assuming resolved walls
    //     else
    ///////        calc_evisc_neutral<true>(fields->sd["evisc"]->data,
    ///////                                 fields->u->data, fields->v->data, fields->w->data,
    ///////                                 fields->u->datafluxbot, fields->v->datafluxbot,
    ///////                                 grid->z, grid->dz, 0, fields->visc); // BvS, for now....
    // assume buoyancy calculation is needed
    // else
    // {
    //     // store the buoyancyflux in tmp1
    //     model->thermo->get_buoyancy_fluxbot(fields->atmp["tmp1"]);
    //     // retrieve the full field in tmp1 and use tmp2 for temporary calculations
    //     model->thermo->get_thermo_field(fields->atmp["tmp1"], fields->atmp["tmp2"], "N2", false);
    //     // model->thermo->getThermoField(fields->sd["tmp1"], fields->sd["tmp2"], "b");

    //     calc_evisc(fields->sd["evisc"]->data,
    //                fields->u->data, fields->v->data, fields->w->data, fields->atmp["tmp1"]->data,
    //                fields->u->datafluxbot, fields->v->datafluxbot, fields->atmp["tmp1"]->datafluxbot,
    //                boundaryptr->ustar, boundaryptr->obuk,
    //                grid->z, grid->dz, grid->dzi,
    //                boundaryptr->z0m);
    // }
}
#endif

template class Diff_smag2<double>;
template class Diff_smag2<float>;
