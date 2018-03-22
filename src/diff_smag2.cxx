/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
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
#include <iostream>

#include "grid.h"
#include "fields.h"
#include "master.h"
#include "defines.h"
#include "constants.h"
#include "monin_obukhov.h"
#include "thermo.h"
#include "boundary.h"

#include "diff_smag2.h"

namespace
{
    namespace most = Monin_obukhov;

    enum class Surface_model {Enabled, Disabled};

    template <typename TF, Surface_model surface_model>
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
        const int k_offset = (surface_model == Surface_model::Disabled) ? 0 : 1;
    
        // If the wall isn't resolved, calculate du/dz and dv/dz at lowest grid height using MO
        if (surface_model == Surface_model::Enabled)
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

    template <typename TF, Surface_model surface_model>
    void calc_evisc_neutral(TF* restrict evisc,
                            TF* restrict u, TF* restrict v, TF* restrict w,
                            TF* restrict ufluxbot, TF* restrict vfluxbot,
                            const TF* restrict z, const TF* restrict dz, const TF z0m, const TF mvisc,
                            const TF dx, const TF dy, const TF cs,
                            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                            const int icells, const int jcells, const int ijcells,
                            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int jj = icells;
        const int kk = ijcells;

        // Wall damping constant.
        const int n = 2;
    
        if (surface_model == Surface_model::Disabled)
        {
            for (int k=kstart; k<kend; ++k)
            {
                const double mlen = pow(cs*std::pow(dx*dy*dz[k], 1./3.), 2);
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        evisc[ijk] = mlen * std::sqrt(evisc[ijk]) + mvisc;
                    }
            }
    
            boundary_cyclic.exec(evisc);

            /*
            // For a resolved wall the viscosity at the wall is needed. For now, assume that the eddy viscosity
            // is zero, so set ghost cell such that the viscosity interpolated to the surface equals the molecular viscosity.
            const int kb = kstart;
            const int kt = kend-1;
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ijkb = i + j*jj + kb*kk;
                    const int ijkt = i + j*jj + kt*kk;
                    evisc[ijkb-kk] = 2 * mvisc - evisc[ijkb];
                    evisc[ijkt+kk] = 2 * mvisc - evisc[ijkt];
                }
                */
        }
        else
        {
            for (int k=kstart; k<kend; ++k)
            {
                // Calculate smagorinsky constant times filter width squared, use wall damping according to Mason's paper.
                const double mlen0 = cs*std::pow(dx*dy*dz[k], 1./3.);
                const double mlen  = std::pow(1./(1./std::pow(mlen0, n) + 1./(std::pow(Constants::kappa*(z[k]+z0m), n))), 1./n);
                const double fac   = std::pow(mlen, 2);
    
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        evisc[ijk] = fac * std::sqrt(evisc[ijk]) + mvisc;
                    }
            }
    
            boundary_cyclic.exec(evisc);
        }
    }

    template<typename TF, Surface_model surface_model>
    void calc_evisc(TF* restrict evisc,
                    TF* restrict u, TF* restrict v, TF* restrict w,  TF* restrict N2,
                    TF* restrict ufluxbot, TF* restrict vfluxbot, TF* restrict bfluxbot,
                    TF* restrict ustar, TF* restrict obuk,
                    const TF* restrict z, const TF* restrict dz, const TF* restrict dzi,
                    const TF dx, const TF dy,
                    const TF z0m, const TF mvisc, const TF cs, const TF tPr,
                    const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                    const int icells, const int jcells, const int ijcells,
                    Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int jj = icells;
        const int kk = ijcells;

        if (surface_model == Surface_model::Disabled)
        {
            for (int k=kstart; k<kend; ++k)
            {
                // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
                const TF mlen = cs*std::pow(dx*dy*dz[k], 1./3.);
                const TF fac = std::pow(mlen, 2);
    
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        // Add the buoyancy production to the TKE
                        TF RitPrratio = N2[ijk] / evisc[ijk] / tPr;
                        RitPrratio = std::min(RitPrratio, static_cast<TF>(1.-Constants::dsmall));
                        evisc[ijk] = fac * std::sqrt(evisc[ijk]) * std::sqrt(1.-RitPrratio);
                    }
            }

            /*
            // For a resolved wall the viscosity at the wall is needed. For now, assume that the eddy viscosity
            // is zero, so set ghost cell such that the viscosity interpolated to the surface equals the molecular viscosity.
            const int kb = kstart;
            const int kt = kend-1;
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ijkb = i + j*jj + kb*kk;
                    const int ijkt = i + j*jj + kt*kk;
                    evisc[ijkb-kk] = 2 * mvisc - evisc[ijkb];
                    evisc[ijkt+kk] = 2 * mvisc - evisc[ijkt];
                }
                */
        }
        else
        {
            // Variables for the wall damping.
            const TF n = 2.;
    
            // Bottom boundary, here strain is fully parametrized using MO.
            // Calculate smagorinsky constant times filter width squared, use wall damping according to Mason.
            const TF mlen0 = cs*std::pow(dx*dy*dz[kstart], 1./3.);
            const TF mlen = std::pow(1./(1./std::pow(mlen0, n) + 1./(std::pow(Constants::kappa*(z[kstart]+z0m), n))), 1./n);
            const TF fac = std::pow(mlen, 2);
    
            for (int j=jstart; j<jend; ++j)
            {
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    // TODO use the thermal expansion coefficient from the input later, what to do if there is no buoyancy?
                    // Add the buoyancy production to the TKE
                    TF RitPrratio = -bfluxbot[ij]/(Constants::kappa*z[kstart]*ustar[ij])*most::phih(z[kstart]/obuk[ij]) / evisc[ijk] / tPr;
                    RitPrratio = std::min(RitPrratio, static_cast<TF>(1.-Constants::dsmall));
                    evisc[ijk] = fac * std::sqrt(evisc[ijk]) * std::sqrt(static_cast<TF>(1.)-RitPrratio);
                }
            }
    
            for (int k=kstart+1; k<kend; ++k)
            {
                // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
                const TF mlen0 = cs*std::pow(dx*dy*dz[k], 1./3.);
                const TF mlen = std::pow(1./(1./std::pow(mlen0, n) + 1./(std::pow(Constants::kappa*(z[k]+z0m), n))), 1./n);
                const TF fac = std::pow(mlen, 2);
    
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        // Add the buoyancy production to the TKE
                        TF RitPrratio = N2[ijk] / evisc[ijk] / tPr;
                        RitPrratio = std::min(RitPrratio, static_cast<TF>(1.-Constants::dsmall));
                        evisc[ijk] = fac * std::sqrt(evisc[ijk]) * std::sqrt(1.-RitPrratio);
                    }
            }
        }
   
        boundary_cyclic.exec(evisc);
    }

    template <typename TF, Surface_model surface_model>
    void diff_u(TF* restrict ut,
                const TF* restrict u, const TF* restrict v, const TF* restrict w,
                const TF* restrict dzi, const TF* restrict dzhi, const TF dxi, const TF dyi,
                const TF* restrict evisc,
                const TF* restrict fluxbot, const TF* restrict fluxtop,
                const TF* restrict rhoref, const TF* restrict rhorefh,
                const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                const int jj, const int kk)
    {
        const int ii = 1;
    
        TF eviscn, eviscs, eviscb, evisct;
    
        const int k_offset = (surface_model == Surface_model::Disabled) ? 0 : 1;
       
        if (surface_model == Surface_model::Enabled)
        {
            // bottom boundary
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    eviscn = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
                    eviscs = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
                    evisct = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
                    eviscb = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);
    
                    ut[ijk] +=
                             // du/dx + du/dx
                             + ( evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
                               - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
                             // du/dy + dv/dx
                             + ( eviscn*((u[ijk+jj]-u[ijk   ])*dyi + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                               - eviscs*((u[ijk   ]-u[ijk-jj])*dyi + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                             // du/dz + dw/dx
                             + ( rhorefh[kstart+1] * evisct*((u[ijk+kk]-u[ijk   ])* dzhi[kstart+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
                               + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
                }
    
            // top boundary
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + (kend-1)*kk;
                    eviscn = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
                    eviscs = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
                    evisct = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
                    eviscb = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);
                    ut[ijk] +=
                             // du/dx + du/dx
                             + ( evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
                               - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
                             // du/dy + dv/dx
                             + ( eviscn*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                               - eviscs*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                             // du/dz + dw/dx
                             + (- rhorefh[kend  ] * fluxtop[ij]
                                - rhorefh[kend-1] * eviscb*((u[ijk   ]-u[ijk-kk])* dzhi[kend-1] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhoref[kend-1] * dzi[kend-1];
                }
        }
    
        for (int k=kstart+k_offset; k<kend-k_offset; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    eviscn = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
                    eviscs = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
                    evisct = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
                    eviscb = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);
                    ut[ijk] +=
                             // du/dx + du/dx
                             + ( evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
                               - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
                             // du/dy + dv/dx
                             + ( eviscn*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                               - eviscs*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                             // du/dz + dw/dx
                             + ( rhorefh[k+1] * evisct*((u[ijk+kk]-u[ijk   ])* dzhi[k+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
                               - rhorefh[k  ] * eviscb*((u[ijk   ]-u[ijk-kk])* dzhi[k  ] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhoref[k] * dzi[k];
                }
    }

    template <typename TF, Surface_model surface_model>
    void diff_v(TF* restrict vt,
                const TF* restrict u, const TF* restrict v, const TF* restrict w,
                const TF* restrict dzi, const TF* restrict dzhi, const TF dxi, const TF dyi,
                const TF* restrict evisc,
                TF* restrict fluxbot, TF* restrict fluxtop,
                TF* restrict rhoref, TF* restrict rhorefh,
                const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                const int jj, const int kk)

    {
        const int ii = 1;
    
        TF evisce, eviscw, eviscb, evisct;
    
        const int k_offset = (surface_model == Surface_model::Disabled) ? 0 : 1;
    
        if (surface_model == Surface_model::Enabled)
        {
            // bottom boundary
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    evisce = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
                    eviscw = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
                    evisct = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
                    eviscb = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);
                    vt[ijk] +=
                             // dv/dx + du/dy
                             + ( evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                               - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                             // dv/dy + dv/dy
                             + ( evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
                               - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
                             // dv/dz + dw/dy
                             + ( rhorefh[kstart+1] * evisct*((v[ijk+kk]-v[ijk   ])*dzhi[kstart+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
                               + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
                }
    
            // top boundary
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + (kend-1)*kk;
                    evisce = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
                    eviscw = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
                    evisct = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
                    eviscb = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);
                    vt[ijk] +=
                             // dv/dx + du/dy
                             + ( evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                               - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                             // dv/dy + dv/dy
                             + ( evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
                               - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
                             // dv/dz + dw/dy
                             + (- rhorefh[kend  ] * fluxtop[ij]
                                - rhorefh[kend-1] * eviscb*((v[ijk   ]-v[ijk-kk])*dzhi[kend-1] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) / rhoref[kend-1] * dzi[kend-1];
                }
        }
    
        for (int k=kstart+k_offset; k<kend-k_offset; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    evisce = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
                    eviscw = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
                    evisct = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
                    eviscb = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);
                    vt[ijk] +=
                             // dv/dx + du/dy
                             + ( evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                               - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                             // dv/dy + dv/dy
                             + ( evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
                               - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
                             // dv/dz + dw/dy
                             + ( rhorefh[k+1] * evisct*((v[ijk+kk]-v[ijk   ])*dzhi[k+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
                               - rhorefh[k  ] * eviscb*((v[ijk   ]-v[ijk-kk])*dzhi[k  ] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) / rhoref[k] * dzi[k];
                }
    }
   
    template <typename TF>
    void diff_w(TF* restrict wt,
                const TF* restrict u, const TF* restrict v, const TF* restrict w,
                const TF* restrict dzi, const TF* restrict dzhi, const TF dxi, const TF dyi,
                const TF* restrict evisc,
                const TF* restrict rhoref, const TF* restrict rhorefh,
                const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                const int jj, const int kk)
    {
        const int ii = 1;
    
        TF evisce, eviscw, eviscn, eviscs;
    
        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    evisce = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+ii-kk] + evisc[ijk+ii]);
                    eviscw = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-ii] + evisc[ijk   -kk] + evisc[ijk   ]);
                    eviscn = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+jj-kk] + evisc[ijk+jj]);
                    eviscs = 0.25*(evisc[ijk-jj-kk] + evisc[ijk-jj] + evisc[ijk   -kk] + evisc[ijk   ]);
                    wt[ijk] +=
                             // dw/dx + du/dz
                             + ( evisce*((w[ijk+ii]-w[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-kk])*dzhi[k])
                               - eviscw*((w[ijk   ]-w[ijk-ii])*dxi + (u[ijk   ]-u[ijk+  -kk])*dzhi[k]) ) * dxi
                             // dw/dy + dv/dz
                             + ( eviscn*((w[ijk+jj]-w[ijk   ])*dyi + (v[ijk+jj]-v[ijk+jj-kk])*dzhi[k])
                               - eviscs*((w[ijk   ]-w[ijk-jj])*dyi + (v[ijk   ]-v[ijk+  -kk])*dzhi[k]) ) * dyi
                             // dw/dz + dw/dz
                             + ( rhoref[k  ] * evisc[ijk   ]*(w[ijk+kk]-w[ijk   ])*dzi[k  ]
                               - rhoref[k-1] * evisc[ijk-kk]*(w[ijk   ]-w[ijk-kk])*dzi[k-1] ) / rhorefh[k] * 2.* dzhi[k];
                }
    }
    
    template <typename TF>
    void diff_c(TF* restrict at, const TF* restrict a,
                const TF* restrict dzi, const TF* restrict dzhi, const TF dxidxi, const TF dyidyi,
                const TF* restrict evisc,
                const TF* restrict fluxbot, const TF* restrict fluxtop, 
                const TF* restrict rhoref, const TF* restrict rhorefh, const TF tPr,
                const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                const int jj, const int kk)
    {
        const int ii = 1;
    
        TF evisce, eviscw, eviscn, eviscs, evisct, eviscb;
    
        // bottom boundary
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])/tPr;
                eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])/tPr;
                eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])/tPr;
                eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])/tPr;
                evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])/tPr;
                eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])/tPr;
    
                at[ijk] +=
                         + ( evisce*(a[ijk+ii]-a[ijk   ]) 
                           - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi 
                         + ( eviscn*(a[ijk+jj]-a[ijk   ]) 
                           - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                         + ( rhorefh[kstart+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[kstart+1]
                           + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
            }
    
        for (int k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])/tPr;
                    eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])/tPr;
                    eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])/tPr;
                    eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])/tPr;
                    evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])/tPr;
                    eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])/tPr;
    
                    at[ijk] +=
                             + ( evisce*(a[ijk+ii]-a[ijk   ]) 
                               - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi 
                             + ( eviscn*(a[ijk+jj]-a[ijk   ]) 
                               - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                             + ( rhorefh[k+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[k+1]
                               - rhorefh[k  ] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[k]  ) / rhoref[k] * dzi[k];
                }
    
        // top boundary
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + (kend-1)*kk;
                evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])/tPr;
                eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])/tPr;
                eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])/tPr;
                eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])/tPr;
                evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])/tPr;
                eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])/tPr;
    
                at[ijk] +=
                         + ( evisce*(a[ijk+ii]-a[ijk   ]) 
                           - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi 
                         + ( eviscn*(a[ijk+jj]-a[ijk   ]) 
                           - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                         + (-rhorefh[kend  ] * fluxtop[ij]
                           - rhorefh[kend-1] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[kend-1] ) / rhoref[kend-1] * dzi[kend-1];
            }
    }

    template<typename TF>
    TF calc_dnmul(TF* restrict evisc, const TF* restrict dzi, const TF dxidxi, const TF dyidyi, const TF tPr,
                  const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                  const int jj, const int kk)
    {
        const TF one = 1.;
        const TF tPrfac = std::min(one, tPr);
        TF dnmul = 0;
    
        // get the maximum time step for diffusion
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    dnmul = std::max(dnmul, std::abs(tPrfac*evisc[ijk]*(dxidxi + dyidyi + dzi[k]*dzi[k])));
                }
    
        // get_max(&dnmul);
    
        return dnmul;
    }
} // End namespace.

template<typename TF>
Diff_smag2<TF>::Diff_smag2(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Diff<TF>(masterin, gridin, fieldsin, inputin),
    boundary_cyclic(master, grid)
{
    dnmax = inputin.get_item<double>("diff", "dnmax", "", 0.4  );
    cs    = inputin.get_item<double>("diff", "cs"   , "", 0.23 );
    tPr   = inputin.get_item<double>("diff", "tPr"  , "", 1./3.);

    fields.init_diagnostic_field("evisc", "Eddy viscosity", "m2 s-1");
}

template<typename TF>
Diff_smag2<TF>::~Diff_smag2()
{
}

template<typename TF>
void Diff_smag2<TF>::init()
{
    boundary_cyclic.init();
}

template<typename TF>
Diffusion_type Diff_smag2<TF>::get_switch() const
{
    return swdiff;
}

template<typename TF>
unsigned long Diff_smag2<TF>::get_time_limit(const unsigned long idt, const double dt)
{
    auto& gd = grid.get_grid_data();
    double dnmul = calc_dnmul<TF>(fields.sd["evisc"]->fld.data(), gd.dzi.data(), 1./(gd.dx*gd.dx), 1./(gd.dy*gd.dy), tPr,
                                  gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                                  gd.icells, gd.ijcells);
    master.max(&dnmul, 1);

    // Avoid zero division.
    dnmul = std::max(Constants::dsmall, dnmul);

    return idt * dnmax / (dt * dnmul);
}

template<typename TF>
double Diff_smag2<TF>::get_dn(const double dt)
{
    auto& gd = grid.get_grid_data();
    double dnmul = calc_dnmul<TF>(fields.sd["evisc"]->fld.data(), gd.dzi.data(), 1./(gd.dx*gd.dx), 1./(gd.dy*gd.dy), tPr,
                                  gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                                  gd.icells, gd.ijcells);
    master.max(&dnmul, 1);

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

    diff_u<TF, Surface_model::Enabled>(
            fields.mt["u"]->fld.data(),
            fields.mp["u"]->fld.data(), fields.mp["v"]->fld.data(), fields.mp["w"]->fld.data(),
            gd.dzi.data(), gd.dzhi.data(), 1./gd.dx, 1./gd.dy,
            fields.sd["evisc"]->fld.data(),
            fields.mp["u"]->flux_bot.data(), fields.mp["u"]->flux_top.data(),
            fields.rhoref.data(), fields.rhorefh.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    diff_v<TF, Surface_model::Enabled>(
            fields.mt["v"]->fld.data(),
            fields.mp["u"]->fld.data(), fields.mp["v"]->fld.data(), fields.mp["w"]->fld.data(),
            gd.dzi.data(), gd.dzhi.data(), 1./gd.dx, 1./gd.dy,
            fields.sd["evisc"]->fld.data(),
            fields.mp["v"]->flux_bot.data(), fields.mp["v"]->flux_top.data(),
            fields.rhoref.data(), fields.rhorefh.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    diff_w<TF>(fields.mt["w"]->fld.data(),
               fields.mp["u"]->fld.data(), fields.mp["v"]->fld.data(), fields.mp["w"]->fld.data(),
               gd.dzi.data(), gd.dzhi.data(), 1./gd.dx, 1./gd.dy,
               fields.sd["evisc"]->fld.data(),
               fields.rhoref.data(), fields.rhorefh.data(),
               gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
               gd.icells, gd.ijcells);

    for (auto it : fields.st)
        diff_c<TF>(it.second->fld.data(), fields.sp[it.first]->fld.data(),
                   gd.dzi.data(), gd.dzhi.data(), 1./(gd.dx*gd.dx), 1./(gd.dy*gd.dy),
                   fields.sd["evisc"]->fld.data(),
                   fields.sp[it.first]->flux_bot.data(), fields.sp[it.first]->flux_top.data(),
                   fields.rhoref.data(), fields.rhorefh.data(), tPr,
                   gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                   gd.icells, gd.ijcells);
}

template<typename TF>
void Diff_smag2<TF>::exec_viscosity(Boundary<TF>& boundary, Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();

    // Calculate strain rate using MO for velocity gradients lowest level.
    if (boundary.get_switch() == "surface")
        calc_strain2<TF, Surface_model::Enabled>(
                fields.sd["evisc"]->fld.data(),
                fields.mp["u"]->fld.data(), fields.mp["v"]->fld.data(), fields.mp["w"]->fld.data(),
                fields.mp["u"]->flux_bot.data(), fields.mp["v"]->flux_bot.data(),
                boundary.ustar, boundary.obuk,
                gd.z.data(), gd.dzi.data(), gd.dzhi.data(), 1./gd.dx, 1./gd.dy,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

    // Calculate strain rate using resolved boundaries.
    else
        calc_strain2<TF, Surface_model::Disabled>(
                fields.sd["evisc"]->fld.data(),
                fields.mp["u"]->fld.data(), fields.mp["v"]->fld.data(), fields.mp["w"]->fld.data(),
                fields.mp["u"]->flux_bot.data(), fields.mp["v"]->flux_bot.data(),
                nullptr, nullptr,
                gd.z.data(), gd.dzi.data(), gd.dzhi.data(), 1./gd.dx, 1./gd.dy,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

    // Start with retrieving the stability information
    if (thermo.get_switch() == "0")
    {
         // Calculate eddy viscosity using MO at lowest model level
        if (boundary.get_switch() == "surface")
            calc_evisc_neutral<TF, Surface_model::Enabled>(
                    fields.sd["evisc"]->fld.data(),
                    fields.mp["u"]->fld.data(), fields.mp["v"]->fld.data(), fields.mp["w"]->fld.data(),
                    fields.mp["u"]->flux_bot.data(), fields.mp["v"]->flux_bot.data(),
                    gd.z.data(), gd.dz.data(), boundary.z0m, fields.visc,
                    gd.dx, gd.dy, this->cs,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.jcells, gd.ijcells,
                    boundary_cyclic);

         // Calculate eddy viscosity assuming resolved walls
        else
            calc_evisc_neutral<TF, Surface_model::Disabled>(
                    fields.sd["evisc"]->fld.data(),
                    fields.mp["u"]->fld.data(), fields.mp["v"]->fld.data(), fields.mp["w"]->fld.data(),
                    fields.mp["u"]->flux_bot.data(), fields.mp["v"]->flux_bot.data(),
                    gd.z.data(), gd.dz.data(), boundary.z0m, fields.visc,
                    gd.dx, gd.dy, this->cs,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.jcells, gd.ijcells,
                    boundary_cyclic);
    }
    // assume buoyancy calculation is needed
    else
    {
        // store the buoyancyflux in tmp1
        auto& gd = grid.get_grid_data();
        auto buoy_tmp = fields.get_tmp();
        auto tmp = fields.get_tmp();
        thermo.get_buoyancy_fluxbot(*buoy_tmp);
        thermo.get_thermo_field(*buoy_tmp, "N2", false);

        if (boundary.get_switch() == "surface")
            calc_evisc<TF, Surface_model::Enabled>(
                    fields.sd["evisc"]->fld.data(),
                    fields.mp["u"]->fld.data(), fields.mp["v"]->fld.data(), fields.mp["w"]->fld.data(), buoy_tmp->fld.data(),
                    fields.mp["u"]->flux_bot.data(), fields.mp["v"]->flux_bot.data(), buoy_tmp->flux_bot.data(),
                    boundary.ustar, boundary.obuk,
                    gd.z.data(), gd.dz.data(), gd.dzi.data(),
                    gd.dx, gd.dy,
                    boundary.z0m, fields.visc, this->cs, this->tPr,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.jcells, gd.ijcells,
                    boundary_cyclic);
        else
            calc_evisc<TF, Surface_model::Disabled>(
                    fields.sd["evisc"]->fld.data(),
                    fields.mp["u"]->fld.data(), fields.mp["v"]->fld.data(), fields.mp["w"]->fld.data(), buoy_tmp->fld.data(),
                    fields.mp["u"]->flux_bot.data(), fields.mp["v"]->flux_bot.data(), buoy_tmp->flux_bot.data(),
                    boundary.ustar, boundary.obuk,
                    gd.z.data(), gd.dz.data(), gd.dzi.data(),
                    gd.dx, gd.dy,
                    boundary.z0m, fields.visc, this->cs, this->tPr,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.jcells, gd.ijcells,
                    boundary_cyclic);

        fields.release_tmp(buoy_tmp);
        fields.release_tmp(tmp);
    }
}
#endif

template class Diff_smag2<double>;
template class Diff_smag2<float>;
