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
#include <iostream>

#include "grid.h"
#include "fields.h"
#include "master.h"
#include "defines.h"
#include "constants.h"
#include "monin_obukhov.h"
#include "thermo.h"
#include "boundary.h"
#include "stats.h"
#include "fast_math.h"

#include "diff_deardorff.h"
#include "diffusion_kernels.h"

namespace
{
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;
    namespace dk = Diffusion_kernels;

    // SvdL, 07.06.22, minimum value of the eddy diffusivities (as in DALES). What would be a suitable value here?
    // for now, set this to zero as molecular viscosity is anyway added in the diffusion functions
    template<typename TF> constexpr TF mvisc = 0.;//1e-5;

    template <typename TF>
    void check_for_minval(TF* const restrict a,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells, const int ijcells,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int jj = icells;
        const int kk = ijcells;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // SvdL, 08.07.22: after initialization of sgstke12 field, still do check if values exceed minimum value
                    if ( a[ijk] < Constants::sgstke12_min<TF>)
                        a[ijk] = Constants::sgstke12_min<TF>;
                }

    boundary_cyclic.exec(a);
    }


    // SvdL, 07.06.22: not sure yet how to implement deardorff in resolved walls. For now, only works with surface_model=Enabled
    template <typename TF, Surface_model surface_model>
    void calc_evisc_neutral(
            TF* const restrict evisc,
            const TF* const restrict sgstke12,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* const restrict dzhi,
            const TF* z0m,
            const TF dx, const TF dy, const TF zsize,
            const TF cm, const TF cn, const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells, const int ijcells,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int jj = icells;
        const int kk = ijcells;

        // Wall damping constant.
        constexpr TF n_mason = TF(2.);
        constexpr TF A_vandriest = TF(26.);

        if (surface_model == Surface_model::Disabled)
        {
            //for (int k=kstart; k<kend; ++k)
            //{
            //    // const TF mlen_wall = Constants::kappa<TF>*std::min(z[k], zsize-z[k]);
            //    const TF mlen_smag = cs*std::pow(dx*dy*dz[k], TF(1./3.));

            //    for (int j=jstart; j<jend; ++j)
            //        #pragma ivdep
            //        for (int i=istart; i<iend; ++i)
            //        {
            //            const int ijk_bot = i + j*jj + kstart*kk;
            //            const int ijk_top = i + j*jj + kend*kk;
            //            const TF u_tau_bot = std::pow(
            //                    fm::pow2( visc*(u[ijk_bot] - u[ijk_bot-kk] )*dzhi[kstart] )
            //                  + fm::pow2( visc*(v[ijk_bot] - v[ijk_bot-kk] )*dzhi[kstart] ), TF(0.25) );
            //            const TF u_tau_top = std::pow(
            //                    fm::pow2( visc*(u[ijk_top] - u[ijk_top-kk] )*dzhi[kend] )
            //                  + fm::pow2( visc*(v[ijk_top] - v[ijk_top-kk] )*dzhi[kend] ), TF(0.25) );
            //            const TF fac_bot = TF(1.) - std::exp( -(       z[k] *u_tau_bot) / (A_vandriest*visc) );
            //            const TF fac_top = TF(1.) - std::exp( -((zsize-z[k])*u_tau_top) / (A_vandriest*visc) );
            //            const TF fac = std::min( fac_bot, fac_top );

            //            const int ijk = i + j*jj + k*kk;
            //            evisc[ijk] = fm::pow2(fac * mlen_smag) * std::sqrt(evisc[ijk]);
            //        }
            //}

            //// For a resolved wall the viscosity at the wall is needed. For now, assume that the eddy viscosity
            //// is mirrored around the surface.
            //const int kb = kstart;
            //const int kt = kend-1;
            //for (int j=0; j<jcells; ++j)
            //    #pragma ivdep
            //    for (int i=0; i<icells; ++i)
            //    {
            //        const int ijkb = i + j*jj + kb*kk;
            //        const int ijkt = i + j*jj + kt*kk;
            //        evisc[ijkb-kk] = evisc[ijkb];
            //        evisc[ijkt+kk] = evisc[ijkt];
            //    }
        }
        else
        {
            // SvdL, 07.06.22: this part below should be fine, but has to be tested.
            for (int k=kstart; k<kend; ++k)
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ij = i + j*jj;
                        const int ijk = i + j*jj + k*kk;

                        // Apply Mason's wall correction here, as in DALES
                        const TF fac = std::pow(TF(1.)/(TF(1.)/std::pow(fac, n_mason)
                                  + TF(1.)/(std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);
                        const TF kvisc = cm * fac * sgstke12[ijk];

                        // Also, as in DALES: enforce minimum eddy viscosity (mvisc)
                        evisc[ijk] = std::max(kvisc, mvisc<TF>);
                    }
        }

        boundary_cyclic.exec(evisc);
    }

    // SvdL, 07.06.22: not sure yet how to implement deardorff in resolved walls. For now, only works with surface_model=Enabled
    template<typename TF, Surface_model surface_model>
    void calc_evisc(
            TF* const restrict evisc,
            const TF* const restrict sgstke12,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict N2,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* const restrict dzi,
            const TF* const restrict z0m,
            const TF dx, const TF dy,
            const TF cn, const TF cm,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells, const int ijcells,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int jj = icells;
        const int kk = ijcells;

       // SvdL, 07.06.22: still fix for resolved walls
       if (surface_model == Surface_model::Disabled)
       {
           // for (int k=kstart; k<kend; ++k)
           // {
           //     // calculate smagorinsky constant times filter width squared, do not use wall damping with resolved walls.
           //     const TF mlen = std::pow(dx*dy*dz[k], TF(1./3.));
           //     const TF fac = fm::pow2(mlen);
           //
           //     for (int j=jstart; j<jend; ++j)
           //         #pragma ivdep
           //         for (int i=istart; i<iend; ++i)
           //         {
           //             const int ijk = i + j*jj + k*kk;
           //             // Add the buoyancy production to the TKE
           //             TF RitPrratio = N2[ijk] / evisc[ijk] / tPr;
           //             RitPrratio = std::min(RitPrratio, TF(1.-Constants::dsmall));
           //             evisc[ijk] = fac * std::sqrt(evisc[ijk]) * std::sqrt(TF(1.)-RitPrratio);
           //         }
           // }
           //
           // // For a resolved wall the viscosity at the wall is needed. For now, assume that the eddy viscosity
           // // is mirrored over the surface.
           // const int kb = kstart;
           // const int kt = kend-1;
           // for (int j=0; j<jcells; ++j)
           //     #pragma ivdep
           //     for (int i=0; i<icells; ++i)
           //     {
           //         const int ijkb = i + j*jj + kb*kk;
           //         const int ijkt = i + j*jj + kt*kk;
           //         evisc[ijkb-kk] = evisc[ijkb];
           //         evisc[ijkt+kk] = evisc[ijkt];
           //     }
       }
       // SvdL, 07.06.22: this part below should be fine, but has to be tested.
       else
       {
            // Variables for the wall damping and length scales
            const TF n_mason = 2.;
            TF mlen ;
            TF fac  ;
            TF kvisc;

            /* SvdL: following implementation in DALES, there is no separate treatment for the first level.
            // Bottom boundary, here strain is fully parametrized using MO.
            // Calculate geometric filter width, based on Deardorff (1980)
            const TF mlen0 = std::pow(dx*dy*dz[kstart], TF(1./3.));
            const TF mlen = std::pow(TF(1.)/(TF(1.)/std::pow(mlen0, n) + TF(1.)/(std::pow(Constants::kappa<TF>*(z[kstart]+z0m), n))), TF(1.)/n);
            const TF fac = fm::pow2(mlen);

            for (int j=jstart; j<jend; ++j)
            {
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    // TODO use the thermal expansion coefficient from the input later, what to do if there is no buoyancy?
                    // Add the buoyancy production to the TKE
                    TF RitPrratio = -bfluxbot[ij]/(Constants::kappa<TF>*z[kstart]*ustar[ij])*most::phih(z[kstart]/obuk[ij]) / evisc[ijk] / tPr;
                    RitPrratio = std::min(RitPrratio, TF(1.-Constants::dsmall));
                    evisc[ijk] = fac * std::sqrt(evisc[ijk]) * std::sqrt(TF(1.)-RitPrratio);
                }
            }
            */

            // tentativechange, SvdL: counter loopt hier direct vanaf kstart
            for (int k=kstart; k<kend; ++k) // See above, this counter thus starts at kstart
            {
                // Calculate geometric filter width, based on Deardorff (1980)
                const TF mlen0 = std::pow(dx*dy*dz[k], TF(1./3.));

                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ij = i + j*jj;
                        const int ijk = i + j*jj + k*kk;

                        mlen = mlen0; // reset at every iteration of j,i

                        if( N2[ijk] > 0 ) // only if stably stratified, adapt length scale
                        {
                            mlen = cn * sgstke12[ijk] / std::sqrt(N2[ijk]);
                        }

                        fac  = std::min(mlen0, mlen);
                        // Apply Mason's wall correction here, as in DALES
                        fac = std::pow(TF(1.)/(TF(1.)/std::pow(fac, n_mason) + TF(1.)/
                                    (std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);
                        kvisc = cm * fac * sgstke12[ijk];

                        // Formally, this is only Km (momentum). Kh is calculated when needed.
                        // Also, as in DALES: enforce minimum eddy viscosity (mvisc)
                        evisc[ijk] = std::max(kvisc, mvisc<TF>);
                    }
            }

       }

        boundary_cyclic.exec(evisc);
    }

    // SvdL, 07.06.22: not sure yet how to implement deardorff in resolved walls. For now, only works with surface_model=Enabled
    template<typename TF, Surface_model surface_model>
    void calc_evisc_heat(
            TF* const restrict evisch,
            const TF* const restrict evisc,
            const TF* const restrict sgstke12,
            const TF* const restrict N2,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* const restrict z0m,
            const TF dx, const TF dy,
            const TF cn, const TF ch1, const TF ch2,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells, const int ijcells,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int jj = icells;
        const int kk = ijcells;

       // SvdL, 07.06.22: still fix for resolved walls
       if (surface_model == Surface_model::Disabled)
       {
           // for (int k=kstart; k<kend; ++k)
           // {
           //     // calculate smagorinsky constant times filter width squared, do not use wall damping with resolved walls.
           //     const TF mlen = std::pow(dx*dy*dz[k], TF(1./3.));
           //     const TF fac = fm::pow2(mlen);
           //
           //     for (int j=jstart; j<jend; ++j)
           //         #pragma ivdep
           //         for (int i=istart; i<iend; ++i)
           //         {
           //             const int ijk = i + j*jj + k*kk;
           //             // Add the buoyancy production to the TKE
           //             TF RitPrratio = N2[ijk] / evisc[ijk] / tPr;
           //             RitPrratio = std::min(RitPrratio, TF(1.-Constants::dsmall));
           //             evisc[ijk] = fac * std::sqrt(evisc[ijk]) * std::sqrt(TF(1.)-RitPrratio);
           //         }
           // }
           //
           // // For a resolved wall the viscosity at the wall is needed. For now, assume that the eddy viscosity
           // // is mirrored over the surface.
           // const int kb = kstart;
           // const int kt = kend-1;
           // for (int j=0; j<jcells; ++j)
           //     #pragma ivdep
           //     for (int i=0; i<icells; ++i)
           //     {
           //         const int ijkb = i + j*jj + kb*kk;
           //         const int ijkt = i + j*jj + kt*kk;
           //         evisc[ijkb-kk] = evisc[ijkb];
           //         evisc[ijkt+kk] = evisc[ijkt];
           //     }
       }
       // SvdL, 07.06.22: this part below should be fine, but has to be tested.
       else
       {
            // Variables for the wall damping and length scales
            const TF n_mason = 2.;
            TF mlen ;
            TF fac  ;
            TF kvisc;

            // tentativechange, SvdL: counter loopt hier direct vanaf kstart
            for (int k=kstart; k<kend; ++k) // See above, this counter thus starts at kstart
            {
                // Calculate geometric filter width, based on Deardorff (1980)
                const TF mlen0 = std::pow(dx*dy*dz[k], TF(1./3.));

                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ij = i + j*jj;
                        const int ijk = i + j*jj + k*kk;

                        mlen = mlen0; // reset at every iteration of j,i

                        if( N2[ijk] > 0 ) // only if stably stratified, adapt length scale
                        {
                            mlen = cn * sgstke12[ijk] / std::sqrt(N2[ijk]);
                        }

                        fac  = std::min(mlen0, mlen);
                        // Apply Mason's wall correction here, as in DALES
                        fac = std::pow(TF(1.)/(TF(1.)/std::pow(fac, n_mason) + TF(1.)/
                                    (std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);
                        kvisc = (ch1 + ch2 * fac / mlen0 ) * evisc[ijk];

                        // Also, as in DALES: enforce minimum eddy viscosity (mvisc)
                        evisch[ijk] = std::max(kvisc, mvisc<TF>);
                    }
            }
       }

        boundary_cyclic.exec(evisch);
    }

    // SvdL, 07.06.22: functions specifically for the balance of sgs tke are introduced here.
    template <typename TF>
    void sgstke_shear_tend(
            TF* const restrict at,
            const TF* const restrict a,
            const TF* const restrict evisc,
            const TF* const restrict strain2,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // Calculate shear production of TKE based on Deardorff (1980)
                    // S2 appears to be simply 2*(strain squared). Therefore, use strain2[ijk] ( * 2 which
                    // drops out again after division with 2*sgstke12)
                    // https://github.com/dalesteam/dales/blob/master/src/modsubgrid.f90

                    at[ijk] += evisc[ijk] * strain2[ijk] / a[ijk];
                }

    }

    template <typename TF>
    void sgstke_buoy_tend(
            TF* const restrict at,
            const TF* const restrict a,
            const TF* const restrict evisch,
            const TF* const restrict N2,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // Calculate buoyancy destruction of TKE based on Deardorff (1980)
                    at[ijk] -= evisch[ijk] * N2[ijk] / a[ijk] / TF(2.);
                }

    }

    template <typename TF>
    void sgstke_diss_tend(
            TF* const restrict at,
            const TF* const restrict a,
            const TF* const restrict N2,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* const restrict z0m,
            const TF dx, const TF dy,
            const TF cn,
            const TF ce1, const TF ce2,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
      const TF n_mason = 2.;
      TF mlen ;
      TF fac  ;

        for (int k=kstart; k<kend; ++k) // See above, this counter thus starts at kstart
        {
            // Calculate geometric filter width, based on Deardorff (1980)
            const TF mlen0 = std::pow(dx*dy*dz[k], TF(1./3.));

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    const int ijk = i + j*jj + k*kk;

                    mlen = mlen0; // reset at every iteration of j,i

                    if( N2[ijk] > 0 ) // only if stably stratified, adapt length scale
                    {
                        mlen = cn * a[ijk] / std::sqrt(N2[ijk]);
                    }

                    fac  = std::min(mlen0, mlen);
                    // Apply Mason's wall correction here, as in DALES
                    fac = std::pow(TF(1.)/(TF(1.)/std::pow(fac, n_mason) + TF(1.)/
                                (std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);

                    // SvdL, 07.06.22: quite strange (so check later), because why would l (fac) be altered by the Mason correction, but mlen0 not.
                    // Why wouldn't both have to be altered?
                    // Calculate dissipation of TKE based on Deardorff (1980)
                    at[ijk] -= (ce1 + ce2 * fac / mlen0 ) * std::pow(a[ijk], TF(2.)) / fac / TF(2.);
                }
        }
    }

    template <typename TF>
    void sgstke_diss_tend_neutral(
            TF* const restrict at,
            const TF* const restrict a,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* const restrict z0m,
            const TF dx, const TF dy,
            const TF ce1, const TF ce2,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const TF n_mason = 2.;

        for (int k=kstart; k<kend; ++k) // See above, this counter thus starts at kstart
        {
            // Calculate geometric filter width, based on Deardorff (1980)
            const TF mlen0 = std::pow(dx*dy*dz[k], TF(1./3.));

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    const int ijk = i + j*jj + k*kk;

                    // Apply Mason's wall correction here, as in DALES
                    const TF fac = std::pow(TF(1.)/(TF(1.)/std::pow(mlen0, n_mason) + TF(1.)/
                                (std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);

                    // SvdL, 07.06.22: quite strange (so check later), because why would l (fac) be altered by the Mason correction, but mlen0 not.
                    // Why wouldn't both have to be altered?
                    // Calculate dissipation of TKE based on Deardorff (1980)
                    at[ijk] -= (ce1 + ce2 * fac / mlen0 ) * std::pow(a[ijk], TF(2.)) / fac / TF(2.);
                }
        }
    }

    // SvdL, 07.06.22: why is the molecular viscosity explicitly added everywhere to this? Then minimum on evisc(h) as above should be unnecessary?
    // Additionally, how does it combine with minimum of sgstke12? This makes sure zero divisions are not possible, but additionally already sets a minimum diffusivity?
    template <typename TF, Surface_model surface_model>
    void diff_u(
            TF* const restrict ut,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict dzi,
            const TF* const restrict dzhi,
            const TF dxi, const TF dyi,
            const TF* const restrict evisc,
            const TF* const restrict fluxbot,
            const TF* const restrict fluxtop,
            const TF* const restrict rhoref,
            const TF* const restrict rhorefh,
            const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        constexpr int k_offset = (surface_model == Surface_model::Disabled) ? 0 : 1;

        const int ii = 1;

        if (surface_model == Surface_model::Enabled)
        {
            // bottom boundary
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    const TF evisce = evisc[ijk   ] + visc;
                    const TF eviscw = evisc[ijk-ii] + visc;
                    const TF eviscn = TF(0.25)*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]) + visc;
                    const TF eviscs = TF(0.25)*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]) + visc;
                    const TF evisct = TF(0.25)*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]) + visc;

                    ut[ijk] +=
                             // du/dx + du/dx
                             + ( evisce*(u[ijk+ii]-u[ijk   ])*dxi
                               - eviscw*(u[ijk   ]-u[ijk-ii])*dxi ) * TF(2.)*dxi
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
                    const TF evisce = evisc[ijk   ] + visc;
                    const TF eviscw = evisc[ijk-ii] + visc;
                    const TF eviscn = TF(0.25)*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]) + visc;
                    const TF eviscs = TF(0.25)*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]) + visc;
                    const TF eviscb = TF(0.25)*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]) + visc;

                    ut[ijk] +=
                             // du/dx + du/dx
                             + ( evisce*(u[ijk+ii]-u[ijk   ])*dxi
                               - eviscw*(u[ijk   ]-u[ijk-ii])*dxi ) * TF(2.)*dxi
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
                    const TF evisce = evisc[ijk   ] + visc;
                    const TF eviscw = evisc[ijk-ii] + visc;
                    const TF eviscn = TF(0.25)*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]) + visc;
                    const TF eviscs = TF(0.25)*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]) + visc;
                    const TF evisct = TF(0.25)*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]) + visc;
                    const TF eviscb = TF(0.25)*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]) + visc;
                    ut[ijk] +=
                             // du/dx + du/dx
                             + ( evisce*(u[ijk+ii]-u[ijk   ])*dxi
                               - eviscw*(u[ijk   ]-u[ijk-ii])*dxi ) * TF(2.)*dxi
                             // du/dy + dv/dx
                             + ( eviscn*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                               - eviscs*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                             // du/dz + dw/dx
                             + ( rhorefh[k+1] * evisct*((u[ijk+kk]-u[ijk   ])* dzhi[k+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
                               - rhorefh[k  ] * eviscb*((u[ijk   ]-u[ijk-kk])* dzhi[k  ] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhoref[k] * dzi[k];
                }
    }

    template <typename TF, Surface_model surface_model>
    void diff_v(
            TF* const restrict vt,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict dzi,
            const TF* const restrict dzhi,
            const TF dxi, const TF dyi,
            const TF* const restrict evisc,
            const TF* const restrict fluxbot,
            const TF* const restrict fluxtop,
            const TF* const restrict rhoref,
            const TF* const restrict rhorefh,
            const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        constexpr int k_offset = (surface_model == Surface_model::Disabled) ? 0 : 1;

        const int ii = 1;

        if (surface_model == Surface_model::Enabled)
        {
            // bottom boundary
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    const TF evisce = TF(0.25)*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]) + visc;
                    const TF eviscw = TF(0.25)*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]) + visc;
                    const TF eviscn = evisc[ijk   ] + visc;
                    const TF eviscs = evisc[ijk-jj] + visc;
                    const TF evisct = TF(0.25)*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]) + visc;

                    vt[ijk] +=
                             // dv/dx + du/dy
                             + ( evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                               - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                             // dv/dy + dv/dy
                             + ( eviscn*(v[ijk+jj]-v[ijk   ])*dyi
                               - eviscs*(v[ijk   ]-v[ijk-jj])*dyi ) * TF(2.)*dyi
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
                    const TF evisce = TF(0.25)*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]) + visc;
                    const TF eviscw = TF(0.25)*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]) + visc;
                    const TF eviscn = evisc[ijk   ] + visc;
                    const TF eviscs = evisc[ijk-jj] + visc;
                    const TF eviscb = TF(0.25)*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]) + visc;

                    vt[ijk] +=
                             // dv/dx + du/dy
                             + ( evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                               - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                             // dv/dy + dv/dy
                             + ( eviscn*(v[ijk+jj]-v[ijk   ])*dyi
                               - eviscs*(v[ijk   ]-v[ijk-jj])*dyi ) * TF(2.)*dyi
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
                    const TF evisce = TF(0.25)*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]) + visc;
                    const TF eviscw = TF(0.25)*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]) + visc;
                    const TF eviscn = evisc[ijk   ] + visc;
                    const TF eviscs = evisc[ijk-jj] + visc;
                    const TF evisct = TF(0.25)*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]) + visc;
                    const TF eviscb = TF(0.25)*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]) + visc;
                    vt[ijk] +=
                             // dv/dx + du/dy
                             + ( evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                               - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                             // dv/dy + dv/dy
                             + ( eviscn*(v[ijk+jj]-v[ijk   ])*dyi
                               - eviscs*(v[ijk   ]-v[ijk-jj])*dyi ) * TF(2.)*dyi
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
                const TF visc,
                const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                const int jj, const int kk)
    {
        const int ii = 1;

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const TF evisce = TF(0.25)*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+ii-kk] + evisc[ijk+ii]) + visc;
                    const TF eviscw = TF(0.25)*(evisc[ijk-ii-kk] + evisc[ijk-ii] + evisc[ijk   -kk] + evisc[ijk   ]) + visc;
                    const TF eviscn = TF(0.25)*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+jj-kk] + evisc[ijk+jj]) + visc;
                    const TF eviscs = TF(0.25)*(evisc[ijk-jj-kk] + evisc[ijk-jj] + evisc[ijk   -kk] + evisc[ijk   ]) + visc;
                    const TF evisct = evisc[ijk   ] + visc;
                    const TF eviscb = evisc[ijk-kk] + visc;
                    wt[ijk] +=
                             // dw/dx + du/dz
                             + ( evisce*((w[ijk+ii]-w[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-kk])*dzhi[k])
                               - eviscw*((w[ijk   ]-w[ijk-ii])*dxi + (u[ijk   ]-u[ijk+  -kk])*dzhi[k]) ) * dxi
                             // dw/dy + dv/dz
                             + ( eviscn*((w[ijk+jj]-w[ijk   ])*dyi + (v[ijk+jj]-v[ijk+jj-kk])*dzhi[k])
                               - eviscs*((w[ijk   ]-w[ijk-jj])*dyi + (v[ijk   ]-v[ijk+  -kk])*dzhi[k]) ) * dyi
                             // dw/dz + dw/dz
                             + ( rhoref[k  ] * evisct*(w[ijk+kk]-w[ijk   ])*dzi[k  ]
                               - rhoref[k-1] * eviscb*(w[ijk   ]-w[ijk-kk])*dzi[k-1] ) / rhorefh[k] * TF(2.)*dzhi[k];
                }
    }

    template <typename TF, Surface_model surface_model>
    void diff_c(
            TF* const restrict at,
            const TF* const restrict a,
            const TF* const restrict dzi,
            const TF* const restrict dzhi,
            const TF dxidxi, const TF dyidyi,
            const TF* const restrict evisc,
            const TF* const restrict fluxbot,
            const TF* const restrict fluxtop,
            const TF* const restrict rhoref,
            const TF* const restrict rhorefh,
            const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        constexpr int k_offset = (surface_model == Surface_model::Disabled) ? 0 : 1;

        const int ii = 1;

        if (surface_model == Surface_model::Enabled)
        {
            // bottom boundary
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    const TF evisce = TF(0.5)*(evisc[ijk   ]+evisc[ijk+ii]) + visc;
                    const TF eviscw = TF(0.5)*(evisc[ijk-ii]+evisc[ijk   ]) + visc;
                    const TF eviscn = TF(0.5)*(evisc[ijk   ]+evisc[ijk+jj]) + visc;
                    const TF eviscs = TF(0.5)*(evisc[ijk-jj]+evisc[ijk   ]) + visc;
                    const TF evisct = TF(0.5)*(evisc[ijk   ]+evisc[ijk+kk]) + visc;

                    at[ijk] +=
                             + ( evisce*(a[ijk+ii]-a[ijk   ])
                               - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                             + ( eviscn*(a[ijk+jj]-a[ijk   ])
                               - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                             + ( rhorefh[kstart+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[kstart+1]
                               + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
                }

            // top boundary
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + (kend-1)*kk;
                    const TF evisce = TF(0.5)*(evisc[ijk   ]+evisc[ijk+ii]) + visc;
                    const TF eviscw = TF(0.5)*(evisc[ijk-ii]+evisc[ijk   ]) + visc;
                    const TF eviscn = TF(0.5)*(evisc[ijk   ]+evisc[ijk+jj]) + visc;
                    const TF eviscs = TF(0.5)*(evisc[ijk-jj]+evisc[ijk   ]) + visc;
                    const TF eviscb = TF(0.5)*(evisc[ijk-kk]+evisc[ijk   ]) + visc;

                    at[ijk] +=
                             + ( evisce*(a[ijk+ii]-a[ijk   ])
                               - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                             + ( eviscn*(a[ijk+jj]-a[ijk   ])
                               - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                             + (-rhorefh[kend  ] * fluxtop[ij]
                               - rhorefh[kend-1] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[kend-1] ) / rhoref[kend-1] * dzi[kend-1];
                }
        }

        for (int k=kstart+k_offset; k<kend-k_offset; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const TF evisce = TF(0.5)*(evisc[ijk   ]+evisc[ijk+ii]) + visc;
                    const TF eviscw = TF(0.5)*(evisc[ijk-ii]+evisc[ijk   ]) + visc;
                    const TF eviscn = TF(0.5)*(evisc[ijk   ]+evisc[ijk+jj]) + visc;
                    const TF eviscs = TF(0.5)*(evisc[ijk-jj]+evisc[ijk   ]) + visc;
                    const TF evisct = TF(0.5)*(evisc[ijk   ]+evisc[ijk+kk]) + visc;
                    const TF eviscb = TF(0.5)*(evisc[ijk-kk]+evisc[ijk   ]) + visc;

                    at[ijk] +=
                             + ( evisce*(a[ijk+ii]-a[ijk   ])
                               - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                             + ( eviscn*(a[ijk+jj]-a[ijk   ])
                               - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                             + ( rhorefh[k+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[k+1]
                               - rhorefh[k  ] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[k]  ) / rhoref[k] * dzi[k];
                }
    }

    // SvdL, 07.06.22: Because of the rewriting in sqrt(sgstke), an additional factor two comes in the diffusion (see DALES paper).
    // this factor then cancels with the 0.5 in the averaging of eddy viscosity (see below). Therefore, a separate diffusion function is required for the diffusion of sgstke12
    template <typename TF, Surface_model surface_model>
    void diff_sgstke12(
            TF* const restrict at,
            const TF* const restrict a,
            const TF* const restrict dzi,
            const TF* const restrict dzhi,
            const TF dxidxi, const TF dyidyi,
            const TF* const restrict evisc,
            const TF* const restrict fluxbot,
            const TF* const restrict fluxtop,
            const TF* const restrict rhoref,
            const TF* const restrict rhorefh,
            const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        constexpr int k_offset = (surface_model == Surface_model::Disabled) ? 0 : 1;

        const int ii = 1;

        if (surface_model == Surface_model::Enabled)
        {
            // bottom boundary
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    const TF evisce = (evisc[ijk   ]+evisc[ijk+ii]) + TF(2.)*visc;
                    const TF eviscw = (evisc[ijk-ii]+evisc[ijk   ]) + TF(2.)*visc;
                    const TF eviscn = (evisc[ijk   ]+evisc[ijk+jj]) + TF(2.)*visc;
                    const TF eviscs = (evisc[ijk-jj]+evisc[ijk   ]) + TF(2.)*visc;
                    const TF evisct = (evisc[ijk   ]+evisc[ijk+kk]) + TF(2.)*visc;

                    at[ijk] +=
                             + ( evisce*(a[ijk+ii]-a[ijk   ])
                               - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                             + ( eviscn*(a[ijk+jj]-a[ijk   ])
                               - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                             + ( rhorefh[kstart+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[kstart+1]
                               ) / rhoref[kstart] * dzi[kstart];
                }

            // top boundary
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + (kend-1)*kk;
                    const TF evisce = (evisc[ijk   ]+evisc[ijk+ii]) + TF(2.)*visc;
                    const TF eviscw = (evisc[ijk-ii]+evisc[ijk   ]) + TF(2.)*visc;
                    const TF eviscn = (evisc[ijk   ]+evisc[ijk+jj]) + TF(2.)*visc;
                    const TF eviscs = (evisc[ijk-jj]+evisc[ijk   ]) + TF(2.)*visc;
                    const TF eviscb = (evisc[ijk-kk]+evisc[ijk   ]) + TF(2.)*visc;

                    at[ijk] +=
                             + ( evisce*(a[ijk+ii]-a[ijk   ])
                               - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                             + ( eviscn*(a[ijk+jj]-a[ijk   ])
                               - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                             + (-rhorefh[kend  ] * Constants::sgstke12_min<TF> // SvdL (9.7.21): quick & dirty: just put this to sgstke12_min as DALES. Should be a better way to do this? Still leaves possibility of flux of sgstke12 from top
                               - rhorefh[kend-1] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[kend-1] ) / rhoref[kend-1] * dzi[kend-1];
                             // + (-rhorefh[kend  ] * fluxtop[ij]
                             //   - rhorefh[kend-1] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[kend-1] ) / rhoref[kend-1] * dzi[kend-1];
                }
        }

        for (int k=kstart+k_offset; k<kend-k_offset; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const TF evisce = (evisc[ijk   ]+evisc[ijk+ii]) + TF(2.)*visc;
                    const TF eviscw = (evisc[ijk-ii]+evisc[ijk   ]) + TF(2.)*visc;
                    const TF eviscn = (evisc[ijk   ]+evisc[ijk+jj]) + TF(2.)*visc;
                    const TF eviscs = (evisc[ijk-jj]+evisc[ijk   ]) + TF(2.)*visc;
                    const TF evisct = (evisc[ijk   ]+evisc[ijk+kk]) + TF(2.)*visc;
                    const TF eviscb = (evisc[ijk-kk]+evisc[ijk   ]) + TF(2.)*visc;

                    at[ijk] +=
                             + ( evisce*(a[ijk+ii]-a[ijk   ])
                               - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                             + ( eviscn*(a[ijk+jj]-a[ijk   ])
                               - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                             + ( rhorefh[k+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[k+1]
                               - rhorefh[k  ] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[k]  ) / rhoref[k] * dzi[k];
                }
    }

    // SGS_TODO, SvdL, 07.06.22: maximum dnmul is given by absolute maximum of
    // either evisc and evisch. So would be better to pass both fields and do a
    // double check on maximum instead of passing "wrong" Prandtl number
    template<typename TF>
    TF calc_dnmul(
            TF* const restrict evisc,
            const TF* const restrict dzi,
            const TF dxidxi, const TF dyidyi,
            const TF tPr,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const TF tPrfac = std::min(TF(1.), tPr);
        TF dnmul = 0;

        // get the maximum time step for diffusion
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    dnmul = std::max(dnmul, std::abs(evisc[ijk]/tPrfac*(dxidxi + dyidyi + dzi[k]*dzi[k])));
                }

        return dnmul;
    }

    //** SvdL (9.07.21): without the proper recalculation of the Prandtl number, these fluxes are not correct??
    // or are they, because the correct viscosity for scalars is passed?
    template <typename TF, Surface_model surface_model>
    void calc_diff_flux_c(
            TF* const restrict out,
            const TF* const restrict data,
            const TF* const restrict evisc,
            const TF* const restrict dzhi,
            const TF tPr, const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        constexpr int k_offset = (surface_model == Surface_model::Disabled) ? 0 : 1;

        #pragma omp parallel for
        for (int k=kstart+k_offset; k<(kend+1-k_offset); ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const TF eviscc = 0.5*(evisc[ijk-kk]+evisc[ijk])/tPr + visc;

                    out[ijk] = - eviscc*(data[ijk] - data[ijk-kk])*dzhi[k];
                }
        }
    }

    template <typename TF, Surface_model surface_model>
    void calc_diff_flux_u(
            TF* const restrict out,
            const TF* const restrict data,
            const TF* const restrict w,
            const TF* const evisc,
            const TF dxi, const TF* const dzhi,
            const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        constexpr int k_offset = (surface_model == Surface_model::Disabled) ? 0 : 1;

        const int ii = 1;
        #pragma omp parallel for
        for (int k=kstart+k_offset; k<(kend+1-k_offset); ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const TF eviscu = 0.25*(evisc[ijk-ii-ijcells]+evisc[ijk-ii]+evisc[ijk-ijcells]+evisc[ijk]) + visc;
                    out[ijk] = - eviscu*( (data[ijk]-data[ijk-ijcells])*dzhi[k] + (w[ijk]-w[ijk-ii])*dxi );
                }
        }
    }

    template <typename TF, Surface_model surface_model>
    void calc_diff_flux_v(
            TF* const restrict out,
            const TF* const restrict data,
            const TF* const restrict w,
            const TF* const evisc,
            const TF dyi, const TF* const dzhi,
            const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        constexpr int k_offset = (surface_model == Surface_model::Disabled) ? 0 : 1;

        #pragma omp parallel for
        for (int k=kstart+k_offset; k<(kend+1-k_offset); ++k)
        {
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*icells + k*ijcells;
                        const TF eviscv = 0.25*(evisc[ijk-icells-ijcells]+evisc[ijk-icells]+evisc[ijk-ijcells]+evisc[ijk]) + visc;
                        out[ijk] = - eviscv*( (data[ijk]-data[ijk-ijcells])*dzhi[k] + (w[ijk]-w[ijk-icells])*dyi );
                    }
        }
    }

    template<typename TF>
    void calc_diff_flux_bc(
            TF* const restrict out,
            const TF* const restrict data,
            const int istart, const int iend,
            const int jstart, const int jend, const int k,
            const int icells, const int ijcells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = i + j*icells + k*ijcells;
                out[ijk] = data[ij];
            }
    }

}

template<typename TF>
Diff_deardorff<TF>::Diff_deardorff(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Boundary<TF>& boundaryin, Input& inputin) :
    Diff<TF>(masterin, gridin, fieldsin, boundaryin, inputin),
    boundary_cyclic(master, grid),
    field3d_operators(master, grid, fields)
{
    auto& gd = grid.get_grid_data();
    dnmax = inputin.get_item<TF>("diff", "dnmax", "", 0.4  );

    // Read constants of the Deardorff subgrid tke scheme
    ap    = inputin.get_item<TF>("diff", "ap"   , "", 1.5  );
    cf    = inputin.get_item<TF>("diff", "cf"   , "", 2.5  );
    ce1   = inputin.get_item<TF>("diff", "ce1"  , "", 0.19 );
    ce2   = inputin.get_item<TF>("diff", "ce2"  , "", 0.51 );
    cm    = inputin.get_item<TF>("diff", "cm"   , "", 0.12 );
    ch1   = inputin.get_item<TF>("diff", "ch1"  , "", 1.   );
    ch2   = inputin.get_item<TF>("diff", "ch2"  , "", 2.   );
    cn    = inputin.get_item<TF>("diff", "cn"   , "", 0.76 );

    const std::string group_name = "sgstke";

    // Set the switch between buoy/no buoy once
    const std::string sw_thermo = inputin.get_item<std::string>("thermo", "swthermo", "");
    sw_buoy = (sw_thermo == "0") ? false : true;

    // As in Deardorff (1980) work with square root of sgs-tke:
    fields.init_prognostic_field("sgstke12", "Square root of SGS TKE", "m s-1", group_name, gd.sloc);

    // SGS_TODO, SvdL, 07.06.22: If I remember correctly, exactly this was needed to avoid zero divisions somewhere? ... I'll have to check
    // maybe it was just because of a call to a non-existing variable?
    fields.sp.at("sgstke12")->visc = inputin.get_item<TF>("fields", "svisc", "sgstke12");

    fields.init_diagnostic_field("evisc",  "Eddy viscosity for momentum", "m2 s-1", group_name, gd.sloc);

    if (sw_buoy)
        fields.init_diagnostic_field("eviscs", "Eddy viscosity for scalars", "m2 s-1",  group_name, gd.sloc);

    // Checks on input
    if (grid.get_spatial_order() != Grid_order::Second)
        throw std::runtime_error("Diff_deardorff only runs with second order grids.");

    if (boundary.get_switch() == "default")
        throw std::runtime_error("Diff_deardorff does not support resolved walls.");
}

template<typename TF>
Diff_deardorff<TF>::~Diff_deardorff()
{
}

template<typename TF>
void Diff_deardorff<TF>::init()
{
    boundary_cyclic.init();
}

template<typename TF>
Diffusion_type Diff_deardorff<TF>::get_switch() const
{
    return swdiff;
}

// SGS_TODO, SvdL, 07.06.22: the regular fixed tPr has no meaning in this tke-scheme. Not sure how to best fix it here while maintaining correct fractional timestep.
// Potentially, see comment above (calc_dnmul function): pass both evisc and eviscs to check on maximum on both? More future proof, if Deardorff scheme would later allow for Pr > 1 !

#ifndef USECUDA
template<typename TF>
unsigned long Diff_deardorff<TF>::get_time_limit(const unsigned long idt, const double dt)
{
    auto& gd = grid.get_grid_data();

    const TF tempPr = 1.;
    double dnmul = calc_dnmul<TF>(
            fields.sd.at("evisc")->fld.data(),
            gd.dzi.data(), 1./(gd.dx*gd.dx), 1./(gd.dy*gd.dy), tempPr,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    master.max(&dnmul, 1);

    // Avoid zero division.
    dnmul = std::max(Constants::dsmall, dnmul);

    return idt * dnmax / (dt * dnmul);
}
#endif

#ifndef USECUDA
template<typename TF>
double Diff_deardorff<TF>::get_dn(const double dt)
{
    auto& gd = grid.get_grid_data();

    const TF tempPr = 1.;
    double dnmul = calc_dnmul<TF>(
            fields.sd.at("evisc")->fld.data(),
            gd.dzi.data(), 1./(gd.dx*gd.dx), 1./(gd.dy*gd.dy), tempPr,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    master.max(&dnmul, 1);

    return dnmul*dt;
}
#endif


// SvdL, 07.06.22: What is the maximum viscosity exactly and where is it set?
template<typename TF>
void Diff_deardorff<TF>::create(Stats<TF>& stats)
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

    create_stats(stats);

    // tentativechange, SvdL, 08.06.22: maybe also not the best spot, but check here if initialized sgstke12 is non-zero (> minvalue).
    // if not, force minimum value here at start of simulation. Field should have been filled already.
    // this is actually an ugly position, but alternative to writing additional lines/calls in model.cxx or fields.cxx
    check_for_minval<TF>(fields.sp.at("sgstke12")->fld.data(),
                      gd.istart, gd.iend,
                      gd.jstart, gd.jend,
                      gd.kstart, gd.kend,
                      gd.icells, gd.jcells, gd.ijcells,
                      boundary_cyclic);
}

#ifndef USECUDA
template<typename TF>
void Diff_deardorff<TF>::exec(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    diff_u<TF, Surface_model::Enabled>(
            fields.mt.at("u")->fld.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dzhi.data(),
            1./gd.dx, 1./gd.dy,
            fields.sd.at("evisc")->fld.data(),
            fields.mp.at("u")->flux_bot.data(),
            fields.mp.at("u")->flux_top.data(),
            fields.rhoref.data(), fields.rhorefh.data(),
            fields.visc,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    diff_v<TF, Surface_model::Enabled>(
            fields.mt.at("v")->fld.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dzhi.data(),
            1./gd.dx, 1./gd.dy,
            fields.sd.at("evisc")->fld.data(),
            fields.mp.at("v")->flux_bot.data(),
            fields.mp.at("v")->flux_top.data(),
            fields.rhoref.data(), fields.rhorefh.data(),
            fields.visc,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    diff_w<TF>(
            fields.mt.at("w")->fld.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dzhi.data(),
            1./gd.dx, 1./gd.dy,
            fields.sd.at("evisc")->fld.data(),
            fields.rhoref.data(), fields.rhorefh.data(),
            fields.visc,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    for (auto it : fields.st)
    {
        if( it.first == "sgstke12" ) // sgstke12 diffuses with Km AND separate diffusion function
        {
            diff_sgstke12<TF, Surface_model::Enabled>(
                    it.second->fld.data(),
                    fields.sp.at(it.first)->fld.data(),
                    gd.dzi.data(), gd.dzhi.data(), 1./(gd.dx*gd.dx), 1./(gd.dy*gd.dy),
                    fields.sd.at("evisc")->fld.data(),
                    fields.sp.at(it.first)->flux_bot.data(), fields.sp.at(it.first)->flux_top.data(),
                    fields.rhoref.data(), fields.rhorefh.data(),
                    fields.sp.at(it.first)->visc,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
        }
        else // all other scalars, diffuse with Kh
        {
            if( sw_buoy ) // if no buoyancy, Kh = Km
            {
                diff_c<TF, Surface_model::Enabled>(
                        it.second->fld.data(),
                        fields.sp.at(it.first)->fld.data(),
                        gd.dzi.data(), gd.dzhi.data(),
                        1./(gd.dx*gd.dx), 1./(gd.dy*gd.dy),
                        fields.sd.at("evisc")->fld.data(),
                        fields.sp.at(it.first)->flux_bot.data(),
                        fields.sp.at(it.first)->flux_top.data(),
                        fields.rhoref.data(), fields.rhorefh.data(),
                        fields.sp.at(it.first)->visc,
                        gd.istart, gd.iend,
                        gd.jstart, gd.jend,
                        gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);
            }
            else // assume buoyancy calculation is needed
            {
                diff_c<TF, Surface_model::Enabled>(
                        it.second->fld.data(),
                        fields.sp.at(it.first)->fld.data(),
                        gd.dzi.data(), gd.dzhi.data(),
                        1./(gd.dx*gd.dx), 1./(gd.dy*gd.dy),
                        fields.sd.at("eviscs")->fld.data(),
                        fields.sp.at(it.first)->flux_bot.data(),
                        fields.sp.at(it.first)->flux_top.data(),
                        fields.rhoref.data(), fields.rhorefh.data(),
                        fields.sp.at(it.first)->visc,
                        gd.istart, gd.iend,
                        gd.jstart, gd.jend,
                        gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);
            }
        }
    }

    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);

    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);
}

template<typename TF>
void Diff_deardorff<TF>::exec_viscosity(Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();
    auto str2_tmp = fields.get_tmp();

    // SvdL, 07.06.22: Looking at DALES, the shear_squared just appears to be 2*(strain squared)..
    // That would be a simple fix.. https://github.com/dalesteam/dales/blob/master/src/modsubgrid.f90
    // Which additional fix would be needed again?

    // Calculate strain rate using MO for velocity gradients lowest level.
    const std::vector<TF>& dudz = boundary.get_dudz();
    const std::vector<TF>& dvdz = boundary.get_dvdz();
    const std::vector<TF>& z0m = boundary.get_z0m();

    dk::calc_strain2<TF, Surface_model::Enabled>(
            fields.sd.at("evisc")->fld.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            dudz.data(),
            dvdz.data(),
            gd.z.data(),
            gd.dzi.data(),
            gd.dzhi.data(),
            1./gd.dx, 1./gd.dy,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    // Start with retrieving the stability information
    if ( !sw_buoy )
    {
        // Calculate eddy viscosity using MO at lowest model level
        calc_evisc_neutral<TF, Surface_model::Enabled>(
                fields.sd.at("evisc")->fld.data(),
                fields.sp.at("sgstke12")->fld.data(),
                fields.mp.at("u")->fld.data(),
                fields.mp.at("v")->fld.data(),
                fields.mp.at("w")->fld.data(),
                gd.z.data(), gd.dz.data(),
                gd.dzhi.data(), z0m.data(),
                gd.dx, gd.dy, gd.zsize,
                this->cm, this->cn, fields.visc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.jcells, gd.ijcells,
                boundary_cyclic);

        sgstke_diss_tend_neutral<TF>(
                fields.st.at("sgstke12")->fld.data(),
                fields.sp.at("sgstke12")->fld.data(),
                gd.z.data(),
                gd.dz.data(),
                z0m.data(),
                gd.dx,
                gd.dy,
                this->ce1,
                this->ce2,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }
    else
    {
        // Assume buoyancy calculation is needed
        auto buoy_tmp = fields.get_tmp();
        thermo.get_thermo_field(*buoy_tmp, "N2", false, false);

        calc_evisc<TF, Surface_model::Enabled>(
                fields.sd.at("evisc")->fld.data(),
                fields.sp.at("sgstke12")->fld.data(),
                fields.mp.at("u")->fld.data(),
                fields.mp.at("v")->fld.data(),
                fields.mp.at("w")->fld.data(),
                buoy_tmp->fld.data(),
                gd.z.data(), gd.dz.data(),
                gd.dzi.data(), z0m.data(),
                gd.dx, gd.dy,
                this->cn, this->cm,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.jcells, gd.ijcells,
                boundary_cyclic);

        // SGS_TODO, SvdL, 07.06.22: where is the tmp field here? eviscs is not just a temporary field..
        // Calculate the eddy diffusivity for heat and scalars and store in tmp field
        calc_evisc_heat<TF, Surface_model::Enabled>(
                fields.sd.at("eviscs")->fld.data(),
                fields.sd.at("evisc")->fld.data(),
                fields.sp.at("sgstke12")->fld.data(),
                buoy_tmp->fld.data(),
                gd.z.data(), gd.dz.data(), z0m.data(),
                gd.dx, gd.dy,
                this->cn, this->ch1, this->ch2,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.jcells, gd.ijcells,
                boundary_cyclic);

        // BvS: I left the tendency calculations of sgstke here; feels a bit strange
        // to calculate them in `exec_viscosity`, but otherwise strain^2 has to be
        // recalculated in diff->exec()...
        sgstke_buoy_tend<TF>(
                fields.st.at("sgstke12")->fld.data(),
                fields.sp.at("sgstke12")->fld.data(),
                fields.sd.at("eviscs")->fld.data(),
                buoy_tmp->fld.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        sgstke_diss_tend<TF>(
                fields.st.at("sgstke12")->fld.data(),
                fields.sp.at("sgstke12")->fld.data(),
                buoy_tmp->fld.data(),
                gd.z.data(),
                gd.dz.data(),
                z0m.data(),
                gd.dx, gd.dy,
                this->cn,
                this->ce1,
                this->ce2,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        fields.release_tmp(buoy_tmp);
    }

    // SGS_TODO, SvdL, 07.06.22: indeed check that factor 2 of strain^2 is accounted for, or not?
    sgstke_shear_tend<TF>(
            fields.st.at("sgstke12")->fld.data(),
            fields.sp.at("sgstke12")->fld.data(),
            fields.sd.at("evisc")->fld.data(),
            str2_tmp->fld.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    // Release temporary fields
    fields.release_tmp(str2_tmp);
}
#endif

template<typename TF>
void Diff_deardorff<TF>::create_stats(Stats<TF>& stats)
{
    const std::string group_name_tke = "sgstke";
    const std::string group_name_default = "default";

    // Add variables to the statistics
    if (stats.get_switch())
    {
        // Always add statistics of eddy viscosity for momentum (!)
        stats.add_profs(*fields.sd.at("evisc"), "z", {"mean", "2"}, group_name_default);

        // Always add profiles of shear production and dissipation of sgstke12
        stats.add_prof("sgstke12_shear", "Shear production term in (SGS TKE)^(1/2) budget", "m2 s-3", "z" , group_name_tke);
        stats.add_prof("sgstke12_diss", "Dissipation term in (SGS TKE)^(1/2) budget", "m2 s-3", "z" , group_name_tke);

        if (sw_buoy)
        {
            stats.add_profs(*fields.sd.at("eviscs"), "z", {"mean", "2"}, group_name_default);
            stats.add_prof("sgstke12_buoy", "Buoyancy production term in (SGS TKE)^(1/2) budget", "m2 s-3", "z" , group_name_tke);
        }

        stats.add_tendency(*fields.mt.at("u"), "z",  tend_name, tend_longname);
        stats.add_tendency(*fields.mt.at("v"), "z",  tend_name, tend_longname);
        stats.add_tendency(*fields.mt.at("w"), "zh", tend_name, tend_longname);

        for (auto it : fields.st)
            stats.add_tendency(*it.second, "z", tend_name, tend_longname);
    }
}

template<typename TF>
void Diff_deardorff<TF>::exec_stats(Stats<TF>& stats, Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();

    const TF no_offset = 0.;
    const TF no_threshold = 0.;

    stats.calc_stats("evisc", *fields.sd.at("evisc"), no_offset, no_threshold);
    if (sw_buoy)
        stats.calc_stats("eviscs", *fields.sd.at("eviscs"), no_offset, no_threshold);

    // Calculate budget terms
    auto tmp = fields.get_tmp();
    auto strain2 = fields.get_tmp();

    // Calculate strain rate using MO for velocity gradients lowest level.
    const std::vector<TF>& dudz = boundary.get_dudz();
    const std::vector<TF>& dvdz = boundary.get_dvdz();
    const std::vector<TF>& z0m = boundary.get_z0m();

    dk::calc_strain2<TF, Surface_model::Enabled>(
            fields.sd.at("evisc")->fld.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            dudz.data(),
            dvdz.data(),
            gd.z.data(),
            gd.dzi.data(),
            gd.dzhi.data(),
            1./gd.dx, 1./gd.dy,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    //
    // Shear production
    //
    std::fill(tmp->fld.begin(), tmp->fld.end(), TF(0));

    sgstke_shear_tend<TF>(
            tmp->fld.data(),
            fields.sp.at("sgstke12")->fld.data(),
            fields.sd.at("evisc")->fld.data(),
            strain2->fld.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    stats.calc_stats("sgstke12_shear", *tmp, no_offset, no_threshold);

    fields.release_tmp(strain2);

    if (sw_buoy)
    {
        //
        // Dissipation : non-neutral
        //
        std::fill(tmp->fld.begin(), tmp->fld.end(), TF(0));

        auto N2 = fields.get_tmp();
        thermo.get_thermo_field(*N2, "N2", false, false);

        sgstke_diss_tend<TF>(
                tmp->fld.data(),
                fields.sp.at("sgstke12")->fld.data(),
                N2->fld.data(),
                gd.z.data(),
                gd.dz.data(),
                z0m.data(),
                gd.dx,
                gd.dy,
                this->cn,
                this->ce1,
                this->ce2,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_stats("sgstke12_diss", *tmp, no_offset, no_threshold);

        //
        // Buoyancy production/destruction
        //
        std::fill(tmp->fld.begin(), tmp->fld.end(), TF(0));

        sgstke_buoy_tend<TF>(
                tmp->fld.data(),
                fields.sp.at("sgstke12")->fld.data(),
                fields.sd.at("eviscs")->fld.data(),
                N2->fld.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_stats("sgstke12_buoy", *tmp, no_offset, no_threshold);

        fields.release_tmp(N2);
    }
    else
    {
        std::fill(tmp->fld.begin(), tmp->fld.end(), TF(0));

        sgstke_diss_tend_neutral<TF>(
                tmp->fld.data(),
                fields.sp.at("sgstke12")->fld.data(),
                gd.z.data(),
                gd.dz.data(),
                z0m.data(),
                gd.dx, gd.dy,
                this->ce1,
                this->ce2,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_stats("sgstke12_diss", *tmp, no_offset, no_threshold);
    }

    fields.release_tmp(tmp);
}

template<typename TF>
void Diff_deardorff<TF>::diff_flux(
        Field3d<TF>& restrict out, const Field3d<TF>& restrict fld_in)
{
    auto& gd = grid.get_grid_data();

    // SGS_TODO, SvdL, 07.06.22: this can be done much nicer by just adapting these functions completely to avoid need of tPr..
    // Just a dummy value; scalars have their own evisc, so no Prandtl number needed..
    const TF tPr_one = TF(1.);

    // Calculate the boundary fluxes.
    calc_diff_flux_bc(
            out.fld.data(), fld_in.flux_bot.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.icells, gd.ijcells);

    calc_diff_flux_bc(
            out.fld.data(), fld_in.flux_top.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kend, gd.icells, gd.ijcells);

    // Calculate the interior.
    if (fld_in.loc[0] == 1)
        calc_diff_flux_u<TF, Surface_model::Enabled>(
                out.fld.data(), fld_in.fld.data(),
                fields.mp.at("w")->fld.data(),
                fields.sd.at("evisc")->fld.data(),
                gd.dxi, gd.dzhi.data(),
                fields.visc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

    else if (fld_in.loc[1] == 1)
        calc_diff_flux_v<TF, Surface_model::Enabled>(
                out.fld.data(), fld_in.fld.data(),
                fields.mp.at("w")->fld.data(),
                fields.sd.at("evisc")->fld.data(),
                gd.dyi, gd.dzhi.data(),
                fields.visc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    else
    {
        // SvdL, 08.07.22: if no buoyancy all scalars diffuse with Km, in any case sgstke12 has to diffuse with Km
        std::string varname = fld_in.name;
        if (!sw_buoy || varname == "sgstke12" || varname == "w")
            calc_diff_flux_c<TF, Surface_model::Enabled>(
                    out.fld.data(), fld_in.fld.data(),
                    fields.sd.at("evisc")->fld.data(),
                    gd.dzhi.data(),
                    tPr_one, fld_in.visc,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
        else
            calc_diff_flux_c<TF, Surface_model::Enabled>(
                    out.fld.data(), fld_in.fld.data(),
                    fields.sd.at("eviscs")->fld.data(),
                    gd.dzhi.data(),
                    tPr_one, fld_in.visc,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
    }
}

template class Diff_deardorff<double>;
template class Diff_deardorff<float>;
