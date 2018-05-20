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
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <netcdf>
#include "grid.h"
#include "fields.h"
#include "thermo_moist.h"
#include "diff_smag2.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "data_block.h"
#include "model.h"
#include "stats.h"
#include "master.h"
#include "cross.h"
#include "dump.h"
#include "column.h"
#include "thermo_moist_functions.h"
#include "timeloop.h"
#include "field3d_operators.h"


using Finite_difference::O2::interp2;
using Finite_difference::O4::interp4;
using namespace Constants;
using namespace Thermo_moist_functions;

namespace
{

    // Settings for now here for convenience....
    template<typename TF>
    const TF pi      = std::acos(-1.);
    template<typename TF>
    const TF Nc0     = 70e6;     // Fixed cloud droplet number
    template<typename TF>
    const TF K_t     = 2.5e-2;   // Conductivity of heat [J/(sKm)]
    template<typename TF>
    const TF D_v     = 3.e-5;    // Diffusivity of water vapor [m2/s]
    template<typename TF>
    const TF rho_w   = 1.e3;     // Density water
    template<typename TF>
    const TF rho_0   = 1.225;    // SB06, p48
    template<typename TF>
    const TF pirhow  = pi<TF> * rho_w<TF> / 6.;
    template<typename TF>
    const TF mc_min  = 4.2e-15;  // Min mean mass of cloud droplet
    template<typename TF>
    const TF mc_max  = 2.6e-10;  // Max mean mass of cloud droplet
    template<typename TF>
    const TF mr_min  = mc_max<TF>;   // Min mean mass of precipitation drop
    template<typename TF>
    const TF mr_max  = 3e-6;     // Max mean mass of precipitation drop // as in UCLA-LES
    template<typename TF>
    const TF ql_min  = 1.e-6;    // Min cloud liquid water for which calculations are performed
    template<typename TF>
    const TF qr_min  = 1.e-15;   // Min rain liquid water for which calculations are performed

    template<typename TF>
    struct Struct_sat_adjust
    {
        TF ql;
        TF t;
        TF qs;
    };

    template<typename TF>
    inline Struct_sat_adjust<TF> sat_adjust(const TF thl, const TF qt, const TF p, const TF exn)
    {
        int niter = 0;
        int nitermax = 30;
        TF ql, tl, tnr_old = 1.e9, tnr, qs=0;
        tl = thl * exn;
        Struct_sat_adjust<TF> ans;

        // Calculate if q-qs(Tl) <= 0. If so, return 0. Else continue with saturation adjustment
        ans.ql = 0;
        ans.t = tl;
        ans.qs = qsat(p,tl);
        if(qt-ans.qs <= 0)
            return ans;

        tnr = tl;
        while (std::fabs(tnr-tnr_old)/tnr_old> 1e-5 && niter < nitermax)
        {
            ++niter;
            tnr_old = tnr;
            qs = qsat(p,tnr);
            tnr = tnr - (tnr+(Lv/cp)*qs-tl-(Lv/cp)*qt)/(1+(std::pow(Lv,2)*qs)/ (Rv*cp*std::pow(tnr,2)));
        }

        if (niter == nitermax)
        {
            throw std::runtime_error("Non-converging saturation adjustment.");
        }

        ql = std::max(TF(0.),qt - qs);

        ans.ql = ql;
        ans.t = tnr;
        ans.qs = qs;
        return ans;
    }

    template<typename TF>
    void  calc_base_state(TF* restrict pref,    TF* restrict prefh,
                                       TF* restrict rho,     TF* restrict rhoh,
                                       TF* restrict thv,     TF* restrict thvh,
                                       TF* restrict ex,      TF* restrict exh,
                                       TF* restrict thlmean, TF* restrict qtmean, const TF pbot,
                                   const int kstart, const int kend, const TF* restrict z, const TF* restrict dz, const TF* const dzh)
   {

       const TF thlsurf = interp2(thlmean[kstart-1], thlmean[kstart]);
       const TF qtsurf  = interp2(qtmean[kstart-1],  qtmean[kstart]);

       TF ql;

       // Calculate the values at the surface (half level == kstart)
       prefh[kstart] = pbot;
       exh[kstart]   = exner(prefh[kstart]);
       ql            = sat_adjust(thlsurf, qtsurf, prefh[kstart], exh[kstart]).ql;
       thvh[kstart]  = virtual_temperature(exh[kstart], thlsurf, qtsurf, ql);
       rhoh[kstart]  = pbot / (Rd * exh[kstart] * thvh[kstart]);

       // Calculate the first full level pressure
       pref[kstart]  = prefh[kstart] * std::exp(-grav * z[kstart] / (Rd * exh[kstart] * thvh[kstart]));

       for (int k=kstart+1; k<kend+1; ++k)
       {
           // 1. Calculate remaining values (thv and rho) at full-level[k-1]
           ex[k-1]  = exner(pref[k-1]);
           ql       = sat_adjust(thlmean[k-1], qtmean[k-1], pref[k-1], ex[k-1]).ql;
           thv[k-1] = virtual_temperature(ex[k-1], thlmean[k-1], qtmean[k-1], ql);
           rho[k-1] = pref[k-1] / (Rd * ex[k-1] * thv[k-1]);

           // 2. Calculate pressure at half-level[k]
           prefh[k] = prefh[k-1] * std::exp(-grav * dz[k-1] / (Rd * ex[k-1] * thv[k-1]));
           exh[k]   = exner(prefh[k]);

           // 3. Use interpolated conserved quantities to calculate half-level[k] values
           const TF thli = interp2(thlmean[k-1], thlmean[k]);
           const TF qti  = interp2(qtmean [k-1], qtmean [k]);
           const TF qli  = sat_adjust(thli, qti, prefh[k], exh[k]).ql;

           thvh[k]  = virtual_temperature(exh[k], thli, qti, qli);
           rhoh[k]  = prefh[k] / (Rd * exh[k] * thvh[k]);

           // 4. Calculate pressure at full-level[k]
           pref[k] = pref[k-1] * std::exp(-grav * dzh[k] / (Rd * exh[k] * thvh[k]));
       }

       pref[kstart-1] = TF(2.)*prefh[kstart] - pref[kstart];
   }

   template<typename TF>
   void calc_top_and_bot(TF* restrict thl0, TF* restrict qt0, const TF* const z, const TF* const zh, const TF* const dzhi, const int kstart, const int kend)
   {
       // Calculate surface and model top values thl and qt
       TF thl0s, qt0s, thl0t, qt0t;
       thl0s = thl0[kstart] - z[kstart]*(thl0[kstart+1]-thl0[kstart])*dzhi[kstart+1];
       qt0s  = qt0[kstart]  - z[kstart]*(qt0[kstart+1] -qt0[kstart] )*dzhi[kstart+1];
       thl0t = thl0[kend-1] + (zh[kend]-z[kend-1])*(thl0[kend-1]-thl0[kend-2])*dzhi[kend-1];
       qt0t  = qt0[kend-1]  + (zh[kend]-z[kend-1])*(qt0[kend-1]- qt0[kend-2] )*dzhi[kend-1];

       // Set the ghost cells for the reference temperature and moisture
       thl0[kstart-1]  = TF(2.)*thl0s - thl0[kstart];
       thl0[kend]      = TF(2.)*thl0t - thl0[kend-1];
       qt0[kstart-1]   = TF(2.)*qt0s  - qt0[kstart];
       qt0[kend]       = TF(2.)*qt0t  - qt0[kend-1];
    }

   template<typename TF>
   void calc_buoyancy_tend_2nd(TF* restrict wt, TF* restrict thl, TF* restrict qt,
                                             TF* restrict ph, TF* restrict thlh, TF* restrict qth,
                                             TF* restrict ql, TF* restrict thvrefh, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int jj, const int kk)
   {

       TF tl, exnh;

       for (int k=kstart+1; k<kend; k++)
       {
           exnh = exner(ph[k]);
           for (int j=jstart; j<jend; j++)
               #pragma ivdep
               for (int i=istart; i<iend; i++)
               {
                   const int ijk = i + j*jj + k*kk;
                   const int ij  = i + j*jj;
                   thlh[ij] = interp2(thl[ijk-kk], thl[ijk]);
                   qth[ij]  = interp2(qt[ijk-kk], qt[ijk]);
               }
           for (int j=jstart; j<jend; j++)
               #pragma ivdep
               for (int i=istart; i<iend; i++)
               {
                   const int ij  = i + j*jj;
                   ql[ij] = sat_adjust(thlh[ij], qth[ij], ph[k], exnh).ql;
               }
           for (int j=jstart; j<jend; j++)
               #pragma ivdep
               for (int i=istart; i<iend; i++)
               {
                   const int ijk = i + j*jj + k*kk;
                   const int ij  = i + j*jj;
                   wt[ijk] += buoyancy(exnh, thlh[ij], qth[ij], ql[ij], thvrefh[k]);
               }
       }
   }


   template<typename TF>
   void calc_buoyancy(TF* restrict b, TF* restrict thl, TF* restrict qt,
                                    TF* restrict p, TF* restrict ql, TF* restrict thvref, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kcells, const int jj, const int kk)
   {

       TF tl, ex;

       for (int k=0; k<kcells; k++)
       {
           ex = exner(p[k]);
           if (k>=kstart)
           {
               for (int j=jstart; j<jend; j++)
                   #pragma ivdep
                   for (int i=istart; i<iend; i++)
                   {
                       const int ijk = i + j*jj + k*kk;
                       const int ij  = i + j*jj;
                       ql[ijk] = sat_adjust(thl[ijk], qt[ijk], p[k], ex).ql;
                   }
           }
           else
           {
                for (int j=jstart; j<jend; j++)
                   #pragma ivdep
                   for (int i=istart; i<iend; i++)
                   {
                       const int ijk  = i + j*jj+k*kk;
                       ql[ijk] = 0.;
                   }

           }
           for (int j=jstart; j<jend; j++)
               #pragma ivdep
               for (int i=istart; i<iend; i++)
               {
                   const int ijk = i + j*jj + k*kk;
                   const int ij  = i + j*jj;
                   b[ijk] = buoyancy(ex, thl[ijk], qt[ijk], ql[ij], thvref[k]);
               }
       }
   }

  template<typename TF>
  void calc_buoyancy_h(TF* restrict bh, TF* restrict thl,  TF* restrict qt,
                                            TF* restrict ph, TF* restrict thvrefh, TF* restrict thlh, TF* restrict qth,
                                            TF* restrict ql, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int jj, const int kk)
    {
        using Finite_difference::O2::interp2;
        TF tl, exnh;

        for (int k=kstart; k<kend; k++)
        {
            exnh = exner(ph[k]);
            if (k>=kstart)
            {

                for (int j=jstart; j<jend; j++)
                    #pragma ivdep
                    for (int i=istart; i<iend; i++)
                    {
                        const int ijk = i + j*jj + k*kk;
                        const int ij  = i + j*jj;
                        thlh[ij] = interp2(thl[ijk-kk], thl[ijk]);
                        qth[ij]  = interp2(qt[ijk-kk], qt[ijk]);
                    }
                    for (int j=jstart; j<jend; j++)
                        #pragma ivdep
                        for (int i=istart; i<iend; i++)
                        {
                            const int ij  = i + j*jj;
                            ql[ij] = sat_adjust(thlh[ij], qth[ij], ph[k], exnh).ql;
                        }
            }
            else
            {
                for (int j=jstart; j<jend; j++)
                   #pragma ivdep
                   for (int i=istart; i<iend; i++)
                   {
                       const int ij  = i + j*jj;
                       ql[ij] = 0.;
                   }

            }
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ij  = i + j*jj;
                    bh[ijk] = buoyancy(exnh, thlh[ijk], qth[ijk], ql[ij], thvrefh[k]);
                }
        }
   }

   template<typename TF>
   void calc_liquid_water(TF* restrict ql, TF* restrict thl, TF* restrict qt, TF* restrict p, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int jj, const int kk)
   {
       TF ex;
       // Calculate the ql field
       for (int k=kstart; k<kend; k++)
       {
           ex = exner(p[k]);
           for (int j=jstart; j<jend; j++)
               #pragma ivdep
               for (int i=istart; i<iend; i++)
               {
                   const int ijk = i + j*jj + k*kk;
                   ql[ijk] = sat_adjust(thl[ijk], qt[ijk], p[k], ex).ql;
               }
       }
   }

  template<typename TF>
  void calc_liquid_water_h(TF* restrict qlh, TF* restrict thl,  TF* restrict qt,
                                            TF* restrict ph, TF* restrict thlh, TF* restrict qth, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int jj, const int kk)
    {
          using Finite_difference::O2::interp2;
          TF tl, exnh;

          for (int k=kstart+1; k<kend; k++)
          {
              exnh = exner(ph[k]);
              for (int j=jstart; j<jend; j++)
                  #pragma ivdep
                  for (int i=istart; i<iend; i++)
                  {
                      const int ijk = i + j*jj + k*kk;
                      const int ij  = i + j*jj;
                      thlh[ij] = interp2(thl[ijk-kk], thl[ijk]);
                      qth[ij]  = interp2(qt[ijk-kk], qt[ijk]);
                  }
              for (int j=jstart; j<jend; j++)
                  #pragma ivdep
                  for (int i=istart; i<iend; i++)
                  {
                      const int ij  = i + j*jj;
                      const int ijk  = i + j*jj+k*kk;
                      qlh[ijk] = sat_adjust(thlh[ij], qth[ij], ph[k], exnh).ql;
                   }
          }
          for (int j=jstart; j<jend; j++)
              #pragma ivdep
              for (int i=istart; i<iend; i++)
              {
                  const int ijk  = i + j*jj+kstart*kk;
                  qlh[ijk] = 0.;
               }

   }

   template<typename TF>
   void calc_N2(TF* restrict N2, const TF* const restrict thl, const TF* const restrict dzi,
                              TF* restrict thvref, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int jj, const int kk)
   {
       for (int k=kstart; k<kend; ++k)
           for (int j=jstart; j<jend; ++j)
               #pragma ivdep
               for (int i=istart; i<iend; ++i)
               {
                   const int ijk = i + j*jj + k*kk;
                   N2[ijk] = grav/thvref[k]*0.5*(thl[ijk+kk] - thl[ijk-kk])*dzi[k];
               }
   }

   template<typename TF>
   void calc_T(TF* const restrict T, const TF* const restrict thl, const TF* const restrict qt,
               const TF* const restrict pref, const TF* const restrict exnref,
               const int istart, const int iend, const int jstart, const int jend,
               const int jj, const int kk, const int kcells)
   {
       for (int k=0; k<kcells; ++k)
           for (int j=jstart; j<jend; ++j)
               #pragma ivdep
               for (int i=istart; i<iend; ++i)
               {
                   const int ijk = i + j*jj+ k*kk;
                   T[ijk] = sat_adjust(thl[ijk],qt[ijk],pref[k],exnref[k]).t;
                }
   }

   template<typename TF>
   void calc_T_h(TF* restrict Th, TF* restrict thl,  TF* restrict qt,
                                             TF* restrict ph, TF* restrict thlh, TF* restrict qth,
                                             TF* restrict ql, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int jj, const int kk)
     {
           using Finite_difference::O2::interp2;
           TF tl, exnh;

           for (int k=kstart+1; k<kend; k++)
           {
               exnh = exner(ph[k]);
               for (int j=jstart; j<jend; j++)
                   #pragma ivdep
                   for (int i=istart; i<iend; i++)
                   {
                       const int ijk = i + j*jj + k*kk;
                       const int ij  = i + j*jj;
                       thlh[ij] = interp2(thl[ijk-kk], thl[ijk]);
                       qth[ij]  = interp2(qt[ijk-kk], qt[ijk]);
                   }
               for (int j=jstart; j<jend; j++)
                   #pragma ivdep
                   for (int i=istart; i<iend; i++)
                   {
                       const int ij  = i + j*jj;
                       const int ijk  = i + j*jj+k*kk;
                       Th[ijk] = sat_adjust(thlh[ij], qth[ij], ph[k], exnh).t;
                    }
           }
    }

   template<typename TF>
   void calc_T_bot(TF* const restrict T_bot, const TF* const restrict th,
                   const TF* const restrict exnrefh, const TF* const restrict threfh,
                   const int istart, const int iend, const int jstart, const int jend, const int kstart,
                   const int jj, const int kk)
   {
       using Finite_difference::O2::interp2;

       for (int j=jstart; j<jend; ++j)
           #pragma ivdep
           for (int i=istart; i<iend; ++i)
           {
               const int ij = i + j*jj;
               const int ijk = i + j*jj + kstart*kk;
               T_bot[ij] = exnrefh[kstart]*threfh[kstart] + (interp2(th[ijk-kk], th[ijk]) - threfh[kstart]);
           }
   }

   template<typename TF>
   void calc_buoyancy_bot(TF* restrict b,      TF* restrict bbot,
                                        TF* restrict thl,    TF* restrict thlbot,
                                        TF* restrict qt,     TF* restrict qtbot,
                                        TF* restrict thvref, TF* restrict thvrefh, const int icells, const int jcells, const int ijcells, const int kstart)
   {
       // assume no liquid water at the lowest model level
       for (int j=0; j<jcells; j++)
           #pragma ivdep
           for (int i=0; i<icells; i++)
           {
               const int ij  = i + j*icells;
               const int ijk = i + j*icells + kstart*ijcells;
               bbot[ij ] = buoyancy_no_ql(thlbot[ij], qtbot[ij], thvrefh[kstart]);
               b   [ijk] = buoyancy_no_ql(thl[ijk], qt[ijk], thvref[kstart]);
           }
   }

   template<typename TF>
   void calc_buoyancy_fluxbot(TF* restrict bfluxbot, TF* restrict thlbot, TF* restrict thlfluxbot, TF* restrict qtbot, TF* restrict qtfluxbot, TF* restrict thvrefh, const int icells, const int jcells, const int kstart)
   {

       // assume no liquid water at the lowest model level
       for (int j=0; j<jcells; j++)
           #pragma ivdep
           for (int i=0; i<icells; i++)
           {
               const int ij = i + j*icells;
               bfluxbot[ij] = buoyancy_flux_no_ql(thlbot[ij], thlfluxbot[ij], qtbot[ij], qtfluxbot[ij], thvrefh[kstart]);
           }
   }

   template<typename TF>
   void calc_buoyancy_tend_4th(TF* restrict wt, TF* restrict thl,  TF* restrict qt,
                                             TF* restrict ph, TF* restrict thlh, TF* restrict qth,
                                             TF* restrict ql, TF* restrict thvrefh, const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
   {
       const int jj  = icells;
       const int kk1 = 1*ijcells;
       const int kk2 = 2*ijcells;

       TF tl, exnh;

       for (int k=kstart+1; k<kend; k++)
       {
           exnh = exner(ph[k]);
           for (int j=jstart; j<jend; j++)
               #pragma ivdep
               for (int i=istart; i<iend; i++)
               {
                   const int ijk = i + j*jj + k*kk1;
                   const int ij  = i + j*jj;
                   thlh[ij] = interp4(thl[ijk-kk2], thl[ijk-kk1], thl[ijk], thl[ijk+kk1]);
                   qth[ij]  = interp4(qt[ijk-kk2],  qt[ijk-kk1],  qt[ijk],  qt[ijk+kk1]);
                   tl       = thlh[ij] * exnh;
                   // Calculate first estimate of ql using Tl
                   // if ql(Tl)>0, saturation adjustment routine needed
                   ql[ij]  = qth[ij]-qsat(ph[k],tl);
               }
           for (int j=jstart; j<jend; j++)
               #pragma ivdep
               for (int i=istart; i<iend; i++)
               {
                   const int ij = i + j*jj;
                   if (ql[ij]>0)   // already doesn't vectorize because of iteration in sat_adjust()
                       ql[ij] = sat_adjust(thlh[ij], qth[ij], ph[k], exnh).ql;
                   else
                       ql[ij] = 0.;
               }
           for (int j=jstart; j<jend; j++)
               #pragma ivdep
               for (int i=istart; i<iend; i++)
               {
                   const int ijk = i + j*jj + k*kk1;
                   const int ij  = i + j*jj;
                   wt[ijk] += buoyancy(exnh, thlh[ij], qth[ij], ql[ij], thvrefh[k]);
               }
       }
   }

    // Rational tanh approximation
    template<typename TF>
    inline TF tanh2(const TF x)
    {
        return x * (27 + x * x) / (27 + 9 * x * x);
    }

    // Given rain water content (qr), number density (nr) and density (rho)
    // calculate mean mass of rain drop
    template<typename TF>
    inline TF calc_rain_mass(const TF qr, const TF nr, const TF rho)
    {
        //TF mr = rho * qr / (nr + dsmall);
        TF mr = rho * qr / std::max(nr, TF(1.));
        mr        = std::min(std::max(mr, mr_min<TF>), mr_max<TF>);
        return mr;
    }

    // Given mean mass rain drop, calculate mean diameter
    template<typename TF>
    inline TF calc_rain_diameter(const TF mr)
    {
        return pow(mr/pirhow<TF>, TF(1.)/TF(3.));
    }

    // Shape parameter mu_r
    template<typename TF>
    inline TF calc_mu_r(const TF dr)
    {
        //return 1./3.; // SB06
        //return 10. * (1. + tanh(1200 * (dr - 0.0015))); // SS08 (Milbrandt&Yau, 2005) -> similar as UCLA
        return 10. * (TF(1.) + tanh2(TF(1200.) * (dr - TF(0.0015)))); // SS08 (Milbrandt&Yau, 2005) -> similar as UCLA
    }

    // Slope parameter lambda_r
    template<typename TF>
    inline TF calc_lambda_r(const TF mur, const TF dr)
    {
        return pow((mur+3)*(mur+2)*(mur+1), TF(1.)/TF(3.)) / dr;
    }

    template<typename TF>
    inline TF minmod(const TF a, const TF b)
    {
        return copysign(TF(1.), a) * std::max(TF(0.), std::min(std::abs(a), TF(copysign(TF(1.), a))*b));
    }

    template<typename TF>
    inline TF min3(const TF a, const TF b, const TF c)
    {
        return std::min(a,std::min(b,c));
    }

    template<typename TF>
    inline TF max3(const TF a, const TF b, const TF c)
    {
        return std::max(a,std::max(b,c));
    }
}

// Microphysics per 2D (xz) slices
namespace mp2d
{
    // Calculate microphysics properties which are used in multiple routines
    template<typename TF>
    void prepare_microphysics_slice(TF* const restrict rain_mass, TF* const restrict rain_diameter,
                                    TF* const restrict mu_r, TF* const restrict lambda_r,
                                    const TF* const restrict qr, const TF* const restrict nr,
                                    const TF* const restrict rho,
                                    const int istart, const int iend,
                                    const int kstart, const int kend,
                                    const int icells, const int ijcells, const int j)
    {

        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                if(qr[ijk] > qr_min<TF>)
                {
                    rain_mass[ik]     = calc_rain_mass(qr[ijk], nr[ijk], rho[k]);
                    rain_diameter[ik] = calc_rain_diameter(rain_mass[ik]);
                    mu_r[ik]          = calc_mu_r(rain_diameter[ik]);
                    lambda_r[ik]      = calc_lambda_r(mu_r[ik], rain_diameter[ik]);
                }
                else
                {
                    rain_mass[ik]     = 0;
                    rain_diameter[ik] = 0;
                    mu_r[ik]          = 0;
                    lambda_r[ik]      = 0;
                }
            }
    }

    // Evaporation: evaporation of rain drops in unsaturated environment
    template<typename TF>
    void evaporation(TF* const restrict qrt, TF* const restrict nrt,
                     TF* const restrict qtt, TF* const restrict thlt,
                     const TF* const restrict qr, const TF* const restrict nr,
                     const TF* const restrict ql, const TF* const restrict qt, const TF* const restrict thl,
                     const TF* const restrict rho, const TF* const restrict exner, const TF* const restrict p,
                     const TF* const restrict rain_mass, const TF* const restrict rain_diameter,
                     const int istart, const int jstart, const int kstart,
                     const int iend,   const int jend,   const int kend,
                     const int jj, const int kk, const int j)
    {
        const TF lambda_evap = 1.; // 1.0 in UCLA, 0.7 in DALES

        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ik  = i + k*jj;
                const int ijk = i + j*jj + k*kk;

                if(qr[ijk] > qr_min<TF>)
                {
                    const TF mr  = rain_mass[ik];
                    const TF dr  = rain_diameter[ik];

                    const TF T   = thl[ijk] * exner[k] + (Lv * ql[ijk]) / (cp * exner[k]); // Absolute temperature [K]
                    const TF Glv = pow(Rv * T / (esat(T) * D_v<TF>) + (Lv / (K_t<TF> * T)) * (Lv / (Rv * T) - 1), -1); // Cond/evap rate (kg m-1 s-1)?

                    const TF S   = (qt[ijk] - ql[ijk]) / qsat(p[k], T) - 1; // Saturation
                    const TF F   = 1.; // Evaporation excludes ventilation term from SB06 (like UCLA, unimportant term? TODO: test)

                    const TF ev_tend = TF(2.) * pi<TF> * dr * Glv * S * F * nr[ijk] / rho[k];

                    qrt[ijk]  += ev_tend;
                    nrt[ijk]  += lambda_evap * ev_tend * rho[k] / mr;
                    qtt[ijk]  -= ev_tend;
                    thlt[ijk] += Lv / (cp * exner[k]) * ev_tend;
                }
            }
    }

    // Selfcollection & breakup: growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
    template<typename TF>
    void selfcollection_breakup(TF* const restrict nrt, const TF* const restrict qr, const TF* const restrict nr, const TF* const restrict rho,
                                const TF* const restrict rain_mass, const TF* const restrict rain_diameter,
                                const TF* const restrict lambda_r,
                                const int istart, const int jstart, const int kstart,
                                const int iend,   const int jend,   const int kend,
                                const int jj, const int kk, const int j)
    {
        const TF k_rr     = 7.12;   // SB06, p49
        const TF kappa_rr = 60.7;   // SB06, p49
        const TF D_eq     = 0.9e-3; // SB06, list of symbols
        const TF k_br1    = 1.0e3;  // SB06, p50, for 0.35e-3 <= Dr <= D_eq
        const TF k_br2    = 2.3e3;  // SB06, p50, for Dr > D_eq

        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ik  = i + k*jj;
                const int ijk = i + j*jj + k*kk;

                if(qr[ijk] > qr_min<TF>)
                {
                    // Calculate mean rain drop mass and diameter
                    const TF mr      = rain_mass[ik];
                    const TF dr      = rain_diameter[ik];
                    const TF lambdar = lambda_r[ik];

                    // Selfcollection
                    const TF sc_tend = -k_rr * nr[ijk] * qr[ijk]*rho[k] * pow(TF(1.) + kappa_rr / lambdar * pow(pirhow<TF>, TF(1.)/TF(3.)), -9) * pow(rho_0<TF> / rho[k], TF(0.5));
                    nrt[ijk] += sc_tend;

                    // Breakup
                    const TF dDr = dr - D_eq;
                    if(dr > 0.35e-3)
                    {
                        TF phi_br;
                        if(dr <= D_eq)
                            phi_br = k_br1 * dDr;
                        else
                            phi_br = TF(2.) * exp(k_br2 * dDr) - TF(1.);
                        const TF br_tend = -(phi_br + TF(1.)) * sc_tend;
                        nrt[ijk] += br_tend;
                    }
                }
            }
    }

    // Sedimentation from Stevens and Seifert (2008)
    template<typename TF>
    void sedimentation_ss08(TF* const restrict qrt, TF* const restrict nrt,
                            TF* const restrict w_qr,TF* const restrict w_nr,
                            TF* const restrict c_qr,TF* const restrict c_nr,
                            TF* const restrict slope_qr,TF* const restrict slope_nr,
                            TF* const restrict flux_qr,TF* const restrict flux_nr,
                            const TF* const restrict mu_r, const TF* const restrict lambda_r,
                            const TF* const restrict qr, const TF* const restrict nr,
                            const TF* const restrict rho, const TF* const restrict dzi,
                            const TF* const restrict dz, const double dt,
                            const int istart, const int jstart, const int kstart,
                            const int iend,   const int jend,   const int kend,
                            const int icells, const int kcells, const int ijcells, const int j)
    {
        const TF w_max = 9.65; // 9.65=UCLA, 20=SS08, appendix A
        const TF a_R = 9.65;   // SB06, p51
        const TF c_R = 600;    // SB06, p51
        const TF Dv  = 25.0e-6;
        const TF b_R = a_R * exp(c_R*Dv); // UCLA-LES

        const int ikcells = icells * kcells;

        const int kk3d = ijcells;
        const int kk2d = icells;

        // 1. Calculate sedimentation velocity at cell center

        for (int k=kstart; k<kend; k++)
        {
            const TF rho_n = pow(TF(1.2) / rho[k], TF(0.5));
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                if(qr[ijk] > qr_min<TF>)
                {
                    // SS08:
                    w_qr[ik] = std::min(w_max, std::max(TF(0.1), rho_n * a_R - b_R * TF(pow(TF(1.) + c_R/lambda_r[ik], TF(-1.)*(mu_r[ik]+TF(4.))))));
                    w_nr[ik] = std::min(w_max, std::max(TF(0.1), rho_n * a_R - b_R * TF(pow(TF(1.) + c_R/lambda_r[ik], TF(-1.)*(mu_r[ik]+TF(1.))))));
                }
                else
                {
                    w_qr[ik] = 0.;
                    w_nr[ik] = 0.;
                }
            }
        }

        // 1.1 Set one ghost cell to zero
        for (int i=istart; i<iend; i++)
        {
            const int ik1  = i + (kstart-1)*icells;
            const int ik2  = i + (kend    )*icells;
            w_qr[ik1] = w_qr[ik1+kk2d];
            w_nr[ik1] = w_nr[ik1+kk2d];
            w_qr[ik2] = 0;
            w_nr[ik2] = 0;
        }

         // 2. Calculate CFL number using interpolated sedimentation velocity

        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ik  = i + k*icells;
                c_qr[ik] = TF(0.25) * (w_qr[ik-kk2d] + TF(2.)*w_qr[ik] + w_qr[ik+kk2d]) * dzi[k] * dt;
                c_nr[ik] = TF(0.25) * (w_nr[ik-kk2d] + TF(2.)*w_nr[ik] + w_nr[ik+kk2d]) * dzi[k] * dt;
            }

        // 3. Calculate slopes

        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                slope_qr[ik] = minmod(qr[ijk]-qr[ijk-kk3d], qr[ijk+kk3d]-qr[ijk]);
                slope_nr[ik] = minmod(nr[ijk]-nr[ijk-kk3d], nr[ijk+kk3d]-nr[ijk]);
            }

        // Calculate flux
        // Set the fluxes at the top of the domain (kend) to zero
        for (int i=istart; i<iend; i++)
        {
            const int ik  = i + kend*icells;
            flux_qr[ik] = 0;
            flux_nr[ik] = 0;
        }

        for (int k=kend-1; k>kstart-1; k--)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                int kk;
                TF ftot, dzz, cc;

                // q_rain
                kk    = k;  // current grid level
                ftot  = 0;  // cumulative 'flux' (kg m-2)
                dzz   = 0;  // distance from zh[k]
                cc    = std::min(TF(1.), c_qr[ik]);
                while (cc > 0 && kk < kend)
                {
                    const int ikk  = i + kk*icells;
                    const int ijkk = i + j*icells + kk*ijcells;

                    ftot  += rho[kk] * (qr[ijkk] + TF(0.5) * slope_qr[ikk] * (TF(1.)-cc)) * cc * dz[kk];

                    dzz   += dz[kk];
                    kk    += 1;
                    cc     = std::min(TF(1.), c_qr[ikk] - dzz*dzi[kk]);
                }

                // Given flux at top, limit bottom flux such that the total rain content stays >= 0.
                ftot = std::min(ftot, rho[k] * dz[k] * qr[ijk] - flux_qr[ik+icells] * TF(dt));
                flux_qr[ik] = -ftot / dt;

                // number density
                kk    = k;  // current grid level
                ftot  = 0;  // cumulative 'flux'
                dzz   = 0;  // distance from zh[k]
                cc    = std::min(TF(1.), c_nr[ik]);
                while (cc > 0 && kk < kend)
                {
                    const int ikk  = i + kk*icells;
                    const int ijkk = i + j*icells + kk*ijcells;

                    ftot += rho[kk] * (nr[ijkk] + TF(0.5) * slope_nr[ikk] * (TF(1.)-cc)) * cc * dz[kk];

                    dzz   += dz[kk];
                    kk    += 1;
                    cc     = std::min(TF(1.), c_nr[ikk] - dzz*dzi[k]);
                }

                // Given flux at top, limit bottom flux such that the number density stays >= 0.
                ftot = std::min(ftot, rho[k] * dz[k] * nr[ijk] - flux_nr[ik+icells] * TF(dt));
                flux_nr[ik] = -ftot / dt;
            }

        // Calculate tendency
        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                qrt[ijk] += -(flux_qr[ik+kk2d] - flux_qr[ik]) / rho[k] * dzi[k];
                nrt[ijk] += -(flux_nr[ik+kk2d] - flux_nr[ik]) / rho[k] * dzi[k];
            }
    }
}

namespace mp
{
    // Remove negative values from 3D field
    template<typename TF>
    void remove_neg_values(TF* const restrict field,
                           const int istart, const int jstart, const int kstart,
                           const int iend,   const int jend,   const int kend,
                           const int jj, const int kk)
    {
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    field[ijk] = std::max(TF(0.), field[ijk]);
                }
    }

    template<typename TF>
    void zero(TF* const restrict field, const int ncells)
    {
        for (int n=0; n<ncells; ++n)
            field[n] = 0.;
    }


    // Autoconversion: formation of rain drop by coagulating cloud droplets
    template<typename TF>
    void autoconversion(TF* const restrict qrt, TF* const restrict nrt,
                        TF* const restrict qtt, TF* const restrict thlt,
                        const TF* const restrict qr,  const TF* const restrict ql,
                        const TF* const restrict rho, const TF* const restrict exner,
                        const int istart, const int jstart, const int kstart,
                        const int iend,   const int jend,   const int kend,
                        const int jj, const int kk)
    {
        const TF x_star = 2.6e-10;       // SB06, list of symbols, same as UCLA-LES
        const TF k_cc   = 9.44e9;        // UCLA-LES (Long, 1974), 4.44e9 in SB06, p48
        const TF nu_c   = 1;             // SB06, Table 1., same as UCLA-LES
        const TF kccxs  = k_cc / (TF(20.) * x_star) * (nu_c+2)*(nu_c+4) / pow(nu_c+1, 2);

        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    if(ql[ijk] > ql_min<TF>)
                    {
                        const TF xc      = rho[k] * ql[ijk] / Nc0<TF>; // Mean mass of cloud drops [kg]
                        const TF tau     = TF(1.) - ql[ijk] / (ql[ijk] + qr[ijk] + dsmall); // SB06, Eq 5
                        const TF phi_au  = TF(600.) * pow(tau, TF(0.68)) * pow(TF(1.) - pow(tau, TF(0.68)), 3); // UCLA-LES
                        //const TF phi_au  = 400. * pow(tau, 0.7) * pow(1. - pow(tau, 0.7), 3); // SB06, Eq 6
                        const TF au_tend = rho[k] * kccxs * pow(ql[ijk], 2) * pow(xc, 2) *
                                               (TF(1.) + phi_au / pow(TF(1.)-tau, 2)); // SB06, eq 4

                        qrt[ijk]  += au_tend;
                        nrt[ijk]  += au_tend * rho[k] / x_star;
                        qtt[ijk]  -= au_tend;
                        thlt[ijk] += Lv / (cp * exner[k]) * au_tend;
                    }
                }
    }

    // Evaporation: evaporation of rain drops in unsaturated environment
    template<typename TF>
    void evaporation(TF* const restrict qrt, TF* const restrict nrt,
                     TF* const restrict qtt, TF* const restrict thlt,
                     const TF* const restrict qr, const TF* const restrict nr,
                     const TF* const restrict ql, const TF* const restrict qt, const TF* const restrict thl,
                     const TF* const restrict rho, const TF* const restrict exner, const TF* const restrict p,
                     const int istart, const int jstart, const int kstart,
                     const int iend,   const int jend,   const int kend,
                     const int jj, const int kk)
    {
        const TF lambda_evap = 1.; // 1.0 in UCLA, 0.7 in DALES

        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    if(qr[ijk] > qr_min<TF>)
                    {
                        // Calculate mean rain drop mass and diameter
                        const TF mr  = calc_rain_mass(qr[ijk], nr[ijk], rho[k]);
                        const TF dr  = calc_rain_diameter(mr);

                        const TF T   = thl[ijk] * exner[k] + (Lv * ql[ijk]) / (cp * exner[k]); // Absolute temperature [K]
                        const TF Glv = pow(Rv * T / (esat(T) * D_v<TF>) + (Lv / (K_t<TF> * T)) * (Lv / (Rv * T) - TF(1.)), -1); // Cond/evap rate (kg m-1 s-1)?

                        const TF S   = (qt[ijk] - ql[ijk]) / qsat(p[k], T) - TF(1.); // Saturation
                        const TF F   = 1.; // Evaporation excludes ventilation term from SB06 (like UCLA, unimportant term? TODO: test)

                        const TF ev_tend = TF(2.) * pi<TF> * dr * Glv * S * F * nr[ijk] / rho[k];

                        qrt[ijk]  += ev_tend;
                        nrt[ijk]  += lambda_evap * ev_tend * rho[k] / mr;
                        qtt[ijk]  -= ev_tend;
                        thlt[ijk] += Lv / (cp * exner[k]) * ev_tend;
                    }
                }
    }

    // Accreation: growth of raindrops collecting cloud droplets
    template<typename TF>
    void accretion(TF* const restrict qrt, TF* const restrict qtt, TF* const restrict thlt,
                   const TF* const restrict qr,  const TF* const restrict ql,
                   const TF* const restrict rho, const TF* const restrict exner,
                   const int istart, const int jstart, const int kstart,
                   const int iend,   const int jend,   const int kend,
                   const int jj, const int kk)
    {
        const TF k_cr  = 5.25; // SB06, p49

        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    if(ql[ijk] > ql_min<TF> && qr[ijk] > qr_min<TF>)
                    {
                        const TF tau     = TF(1.) - ql[ijk] / (ql[ijk] + qr[ijk]); // SB06, Eq 5
                        const TF phi_ac  = pow(tau / (tau + TF(5e-5)), 4); // SB06, Eq 8
                        const TF ac_tend = k_cr * ql[ijk] *  qr[ijk] * phi_ac * pow(rho_0<TF> / rho[k], TF(0.5)); // SB06, Eq 7

                        qrt[ijk]  += ac_tend;
                        qtt[ijk]  -= ac_tend;
                        thlt[ijk] += Lv / (cp * exner[k]) * ac_tend;
                    }
                }
    }


    // Selfcollection & breakup: growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
    template<typename TF>
    void selfcollection_breakup(TF* const restrict nrt, const TF* const restrict qr, const TF* const restrict nr, const TF* const restrict rho,
                                const int istart, const int jstart, const int kstart,
                                const int iend,   const int jend,   const int kend,
                                const int jj, const int kk)
    {
        const TF k_rr     = 7.12;   // SB06, p49
        const TF kappa_rr = 60.7;   // SB06, p49
        const TF D_eq     = 0.9e-3; // SB06, list of symbols
        const TF k_br1    = 1.0e3;  // SB06, p50, for 0.35e-3 <= Dr <= D_eq
        const TF k_br2    = 2.3e3;  // SB06, p50, for Dr > D_eq

        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    if(qr[ijk] > qr_min<TF>)
                    {
                        // Calculate mean rain drop mass and diameter
                        const TF mr      = calc_rain_mass(qr[ijk], nr[ijk], rho[k]);
                        const TF dr      = calc_rain_diameter(mr);
                        const TF mur     = calc_mu_r(dr);
                        const TF lambdar = calc_lambda_r(mur, dr);

                        // Selfcollection
                        const TF sc_tend = -k_rr * nr[ijk] * qr[ijk]*rho[k] * pow(TF(1.) + kappa_rr / lambdar * pow(pirhow<TF>, TF(1.)/TF(3.)), -9) * pow(rho_0<TF> / rho[k], 0.5);
                        nrt[ijk] += sc_tend;

                        // Breakup
                        const TF dDr = dr - D_eq;
                        if(dr > 0.35e-3)
                        {
                            TF phi_br;
                            if(dr <= D_eq)
                                phi_br = k_br1 * dDr;
                            else
                                phi_br = TF(2.) * exp(k_br2 * dDr) - TF(1.);
                            const TF br_tend = -(phi_br + TF(1.)) * sc_tend;
                            nrt[ijk] += br_tend;
                        }
                    }
                }
    }

    // Sedimentation from Stevens and Seifert (2008)
    template<typename TF>
    void sedimentation_ss08(TF* const restrict qrt, TF* const restrict nrt,
                            TF* const restrict w_qr, TF* const restrict w_nr,
                            TF* const restrict c_qr, TF* const restrict c_nr,
                            TF* const restrict slope_qr, TF* const restrict slope_nr,
                            TF* const restrict flux_qr, TF* const restrict flux_nr,
                            const TF* const restrict qr, const TF* const restrict nr,
                            const TF* const restrict rho, const TF* const restrict dzi,
                            const TF* const restrict dz, const double dt,
                            const int istart, const int jstart, const int kstart,
                            const int iend,   const int jend,   const int kend,
                            const int icells, const int kcells, const int ijcells)
    {
        const TF w_max = 9.65; // 9.65=UCLA, 20=SS08, appendix A
        const TF a_R = 9.65;   // SB06, p51
        const TF c_R = 600;    // SB06, p51
        const TF Dv  = 25.0e-6;
        const TF b_R = a_R * exp(c_R*Dv); // UCLA-LES

        const int ikcells = icells * kcells;

        for (int j=jstart; j<jend; ++j)
        {
            // 1. Calculate sedimentation velocity at cell center
            for (int k=kstart; k<kend; k++)
            {
                const TF rho_n = pow(TF(1.2) / rho[k], TF(0.5));
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int ik  = i + k*icells;

                    if(qr[ijk] > qr_min<TF>)
                    {
                        // Calculate mean rain drop mass and diameter
                        const TF mr      = calc_rain_mass(qr[ijk], nr[ijk], rho[k]);
                        const TF dr      = calc_rain_diameter(mr);
                        const TF mur     = calc_mu_r(dr);
                        const TF lambdar = calc_lambda_r(mur, dr);

                        // SS08:
                        w_qr[ik] = std::min(w_max, std::max(TF(0.1), rho_n * a_R - b_R * TF(pow(TF(1.) + c_R/lambdar, TF(-1.)*(mur+TF(4.))))));
                        w_nr[ik] = std::min(w_max, std::max(TF(0.1), rho_n * a_R - b_R * TF(pow(TF(1.) + c_R/lambdar, TF(-1.)*(mur+TF(1.))))));
                    }
                    else
                    {
                        w_qr[ik] = 0.;
                        w_nr[ik] = 0.;
                    }
                }
            }

            // 1.1 Set one ghost cell to zero
            for (int i=istart; i<iend; i++)
            {
                const int ik1  = i + (kstart-1)*icells;
                const int ik2  = i + (kend    )*icells;
                w_qr[ik1] = w_qr[ik1+icells];
                w_nr[ik1] = w_nr[ik1+icells];
                w_qr[ik2] = 0;
                w_nr[ik2] = 0;
            }

             // 2. Calculate CFL number using interpolated sedimentation velocity

            for (int k=kstart; k<kend; k++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ik  = i + k*icells;
                    c_qr[ik] = TF(0.25) * (w_qr[ik-icells] + TF(2.)*w_qr[ik] + w_qr[ik+icells]) * dzi[k] * dt;
                    c_nr[ik] = TF(0.25) * (w_nr[ik-icells] + TF(2.)*w_nr[ik] + w_nr[ik+icells]) * dzi[k] * dt;
                }

            // 3. Calculate slopes with slope limiter to prevent forming new maxima/minima using slopes

            // Dissable slope limiter near surface, i.e. assume that e.g. qr[kstart]-qr[kstart-1] == qr[kstart+1]-qr[kstart]
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*icells + kstart*ijcells;
                const int ik  = i + kstart*icells;

                slope_qr[ik] = qr[ijk+ijcells] - qr[ijk];
                slope_nr[ik] = nr[ijk+ijcells] - nr[ijk];
            }

            for (int k=kstart+1; k<kend; k++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int ik  = i + k*icells;

                    slope_qr[ik] = minmod(qr[ijk]-qr[ijk-ijcells], qr[ijk+ijcells]-qr[ijk]);
                    slope_nr[ik] = minmod(nr[ijk]-nr[ijk-ijcells], nr[ijk+ijcells]-nr[ijk]);
                }

            // Calculate flux
            // Set the fluxes at the top of the domain (kend) to zero
            for (int i=istart; i<iend; i++)
            {
                const int ik  = i + kend*icells;
                flux_qr[ik] = 0;
                flux_nr[ik] = 0;
            }

            for (int k=kend-1; k>kstart-1; k--)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int ik  = i + k*icells;

                    int kk;
                    TF ftot, dzz, cc;

                    // q_rain
                    kk    = k;  // current grid level
                    ftot  = 0;  // cumulative 'flux' (kg m-2)
                    dzz   = 0;  // distance from zh[k]
                    cc    = std::min(TF(1.), c_qr[ik]);
                    while (cc > 0 && kk < kend)
                    {
                        const int ikk  = i + kk*icells;
                        const int ijkk = i + j*icells + kk*ijcells;

                        ftot += rho[kk] * (qr[ijkk] + TF(0.5) * slope_qr[ikk] * (TF(1.)-cc)) * cc * dz[kk];

                        dzz   += dz[kk];
                        kk    += 1;
                        cc     = std::min(TF(1.), c_qr[ikk] - dzz*dzi[kk]);
                    }

                    // Given flux at top, limit bottom flux such that the total rain content stays >= 0.
                    ftot = std::min(ftot, rho[k] * dz[k] * qr[ijk] - flux_qr[ik+icells] * TF(dt) - TF(dsmall));
                    flux_qr[ik] = -ftot / dt;

                    // number density
                    kk    = k;  // current grid level
                    ftot  = 0;  // cumulative 'flux'
                    dzz   = 0;  // distance from zh[k]
                    cc    = std::min(TF(1.), c_nr[ik]);
                    while (cc > 0 && kk < kend)
                    {
                        const int ikk  = i + kk*icells;
                        const int ijkk = i + j*icells + kk*ijcells;

                        ftot += rho[kk] * (nr[ijkk] + TF(0.5) * slope_nr[ikk] * (TF(1.)-cc)) * cc * dz[kk];

                        dzz   += dz[kk];
                        kk    += 1;
                        cc     = std::min(TF(1.), c_nr[ikk] - dzz*dzi[kk]);
                    }

                    // Given flux at top, limit bottom flux such that the number density stays >= 0.
                    ftot = std::min(ftot, rho[k] * dz[k] * nr[ijk] - flux_nr[ik+icells] * TF(dt) - TF(dsmall));
                    flux_nr[ik] = -ftot / TF(dt);
                }

            // Calculate tendency
            for (int k=kstart; k<kend; k++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int ik  = i + k*icells;

                    qrt[ijk] += -(flux_qr[ik+icells] - flux_qr[ik]) / rho[k] * dzi[k];
                    nrt[ijk] += -(flux_nr[ik+icells] - flux_nr[ik]) / rho[k] * dzi[k];
                }
        }
    }

    // Calculate maximum sedimentation velocity
    template<typename TF>
    TF calc_max_sedimentation_cfl(TF* const restrict w_qr,
                                      const TF* const restrict qr, const TF* const restrict nr,
                                      const TF* const restrict rho, const TF* const restrict dzi,
                                      const double dt,
                                      const int istart, const int jstart, const int kstart,
                                      const int iend,   const int jend,   const int kend,
                                      const int icells, const int ijcells)
    {
        const TF w_max = 9.65; // 9.65=UCLA, 20=SS08, appendix A
        const TF a_R = 9.65;   // SB06, p51
        const TF c_R = 600;    // SB06, p51
        const TF Dv  = 25.0e-6;
        const TF b_R = a_R * exp(c_R*Dv); // UCLA-LES

        // Calculate sedimentation velocity at cell centre
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;

                    if(qr[ijk] > qr_min<TF>)
                    {
                        // Calculate mean rain drop mass and diameter
                        const TF mr      = calc_rain_mass(qr[ijk], nr[ijk], rho[k]);
                        const TF dr      = calc_rain_diameter(mr);
                        const TF mur     = calc_mu_r(dr);
                        const TF lambdar = calc_lambda_r(mur, dr);

                        w_qr[ijk] = std::min(w_max, std::max(TF(0.1), a_R - b_R * TF(pow(TF(1.) + c_R/lambdar, TF(-1.)*(mur+TF(4.))))));
                    }
                    else
                    {
                        w_qr[ijk] = 0.;
                    }
                }

        // Calculate maximum CFL based on interpolated velocity
        TF cfl_max = 1e-5;
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;

                    const TF cfl_qr = TF(0.25) * (w_qr[ijk-ijcells] + TF(2.)*w_qr[ijk] + w_qr[ijk+ijcells]) * dzi[k] * TF(dt);
                    cfl_max = std::max(cfl_max, cfl_qr);
                }

        return cfl_max;
    }

} // End namespace


template<typename TF>
Thermo_moist<TF>::Thermo_moist(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
Thermo<TF>(masterin, gridin, fieldsin, inputin),
boundary_cyclic(masterin, gridin),
field3d_operators(master, grid, fieldsin)
{
    swthermo = "moist";

    // BvS:micro Get microphysics switch, and init rain and number density
    std::string swmicro_in = inputin.get_item<std::string>("thermo", "swmicro", "", "0");
    if(swmicro_in == "0")
    {
        swmicro = Micro_type::disabled;
    }
    else if(swmicro_in == "2mom_warm")
    {
        swmicro = Micro_type::two_mom_warm;

        #ifdef USECUDA
        throw std::runtime_error("swmicro = \"2mom_warm\" not (yet) implemented in CUDA\n");
        #endif

        swmicrobudget = inputin.get_item<bool>("thermo", "swmicrobudget", "", false);
        cflmax_micro  = inputin.get_item<TF>("thermo", "cflmax_micro", "", 2.);

        fields.init_prognostic_field("qr", "Rain water mixing ratio", "kg kg-1");
        fields.init_prognostic_field("nr", "Number density rain", "m-3");
    }
    else
    {
        throw std::runtime_error("Invalid option for \"swmicro\"");
    }

    // Initialize the prognostic fields
    fields.init_prognostic_field("thl", "Liquid water potential temperature", "K");
    fields.init_prognostic_field("qt", "Total water mixing ratio", "kg kg-1");

    // Get the diffusivities of temperature and moisture
    fields.sp.at("thl")->visc = inputin.get_item<TF>("fields", "svisc", "thl");
    fields.sp.at("qt")->visc = inputin.get_item<TF>("fields", "svisc", "qt");


    // Test if the diffusivities of theta and qt are equal, else throw error
    if (fields.sp.at("thl")->visc != fields.sp.at("qt")->visc)
    {
        throw std::runtime_error("The diffusivities of temperature and moisture must be equal\n");
    }

    bs.pbot = inputin.get_item<TF>("thermo", "pbot", "");

    // Get base state option (boussinesq or anelastic)
    std::string swbasestate_in = inputin.get_item<std::string>("thermo", "swbasestate", "", "");
    if (swbasestate_in == "boussinesq")
        bs.swbasestate = Basestate_type::boussinesq;
    else if (swbasestate_in == "anelastic")
        bs.swbasestate = Basestate_type::anelastic;
    else
        throw std::runtime_error("Invalid option for \"swbasestate\"");

    if (grid.swspatialorder == "4" && bs.swbasestate == Basestate_type::anelastic)
    {
        master.print_error("Anelastic mode is not supported for swspatialorder=4\n");
        throw std::runtime_error("Illegal options swbasestate");
    }


    // BvS test for updating hydrostatic prssure during run
    // swupdate..=0 -> initial base state pressure used in saturation calculation
    // swupdate..=1 -> base state pressure updated before saturation calculation
    bs.swupdatebasestate = inputin.get_item<bool>("thermo", "swupdatebasestate", "", false);

    // Time variable surface pressure
    swtimedep_pbot = inputin.get_item<bool>("thermo", "swtimedep_pbot", "", false);
    available_masks.insert(available_masks.end(), {"ql", "qlcore", "qr"});

}

template<typename TF>
Thermo_moist<TF>::~Thermo_moist()
{

}

template<typename TF>
void Thermo_moist<TF>::init()
{
    auto& gd = grid.get_grid_data();

    bs.thl0.resize(gd.kcells);
    bs.qt0.resize(gd.kcells);
    bs.thvref.resize(gd.kcells);
    bs.thvrefh.resize(gd.kcells);
    bs.exnref.resize(gd.kcells);
    bs.exnrefh.resize(gd.kcells);
    bs.pref.resize(gd.kcells);
    bs.prefh.resize(gd.kcells);
}

template<typename TF>
void Thermo_moist<TF>::create(Input& inputin, Data_block& data_block, Stats<TF>& stats, Column<TF>& column, Cross<TF>& cross, Dump<TF>& dump)
{
    auto& gd = grid.get_grid_data();

    // Enable automated calculation of horizontally averaged fields
    if (bs.swupdatebasestate)
        fields.set_calc_mean_profs(true);

    // Calculate the base state profiles. With swupdatebasestate=1, these profiles are updated on every iteration.
    // 1. Take the initial profile as the reference
    data_block.get_vector(bs.thl0, "thl", gd.ktot, 0, 0);
    data_block.get_vector(bs.qt0, "qt", gd.ktot, 0, 0);

    calc_top_and_bot(bs.thl0.data(), bs.qt0.data(), gd.z.data(), gd.zh.data(), gd.dzhi.data(), gd.kstart, gd.kend);

    // 4. Calculate the initial/reference base state
    calc_base_state(bs.pref.data(), bs.prefh.data(), fields.rhoref.data(), fields.rhorefh.data(), bs.thvref.data(), bs.thvrefh.data(), bs.exnref.data(), bs.exnrefh.data(), bs.thl0.data(), bs.qt0.data(), bs.pbot, gd.kstart, gd.kend, gd.z.data(), gd.dz.data(), gd.dzh.data());

    // 5. In Boussinesq mode, overwrite reference temperature and density
    if (bs.swbasestate == Basestate_type::boussinesq)
    {
        bs.thvref0 = inputin.get_item<TF>("thermo", "thvref0", "");

        for (int k=0; k<gd.kcells; ++k)
        {
            fields.rhoref[k]  = 1.;
            fields.rhorefh[k] = 1.;
            bs.thvref[k]      = bs.thvref0;
            bs.thvrefh[k]     = bs.thvref0;
        }
    }

    // 6. Process the time dependent surface pressure
/*    if (swtimedep_pbot == 1)
    {
        const int nerror = inputin->get_time(&timedeppbot, &timedeptime, "pbot");
        if (nerror > 0)
            throw 1;
    }
*/
    // Init the toolbox classes.
    boundary_cyclic.init();

    // Set up output classes
    create_stats(stats);
    create_column(column);
    create_dump(dump);
    create_cross(cross);
}

#ifndef USECUDA
template<typename TF>
void Thermo_moist<TF>::exec(const double dt)
{
    auto& gd = grid.get_grid_data();

    // Re-calculate hydrostatic pressure and exner, pass dummy as rhoref,thvref to prevent overwriting base state
    auto tmp = fields.get_tmp();
    if (bs.swupdatebasestate)
        calc_base_state(bs.pref.data(), bs.prefh.data(),
                        &tmp->fld[0*gd.kcells], &tmp->fld[1*gd.kcells], &tmp->fld[2*gd.kcells], &tmp->fld[3*gd.kcells],
                        bs.exnref.data(), bs.exnrefh.data(), fields.sp.at("thl")->fld_mean.data(), fields.sp.at("qt")->fld_mean.data(), bs.pbot, gd.kstart, gd.kend, gd.z.data(), gd.dz.data(), gd.dzh.data());

    // extend later for gravity vector not normal to surface
    calc_buoyancy_tend_2nd(fields.mt.at("w")->fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), bs.prefh.data(),
                           &tmp->fld[0*gd.ijcells], &tmp->fld[1*gd.ijcells],
                           &tmp->fld[2*gd.ijcells], bs.thvrefh.data(), gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);

   fields.release_tmp(tmp);
    // 2-moment warm microphysics
    if(swmicro == Micro_type::two_mom_warm)
        exec_microphysics(dt);
}
#endif

template<typename TF>
unsigned long Thermo_moist<TF>::get_time_limit(unsigned long idt, const double dt)
{
    auto& gd = grid.get_grid_data();

    if(swmicro == Micro_type::two_mom_warm)
    {
        auto w_qr = fields.get_tmp();
        TF cfl = mp::calc_max_sedimentation_cfl(w_qr->fld.data(), fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),
                                                    fields.rhoref.data(), gd.dzi.data(), dt,
                                                    gd.istart, gd.jstart, gd.kstart,
                                                    gd.iend,   gd.jend,   gd.kend,
                                                    gd.icells, gd.ijcells);
        fields.release_tmp(w_qr);
        master.max(&cfl, 1);
        return idt * cflmax_micro / cfl;
    }
    else
    {
        return Constants::ulhuge;
    }
}

// BvS:micro
template<typename TF>
void Thermo_moist<TF>::exec_microphysics(const double dt)
{
    auto& gd = grid.get_grid_data();

    // Switch to solve certain routines over xz-slices, to reduce calculations recurring in several microphysics routines
    bool per_slice = true;

    // Remove the negative values from the precipitation fields
    mp::remove_neg_values(fields.sp.at("qr")->fld.data(), gd.istart, gd.jstart, gd.kstart, gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);
    mp::remove_neg_values(fields.sp.at("nr")->fld.data(), gd.istart, gd.jstart, gd.kstart, gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);

    auto ql = fields.get_tmp();
    // Calculate the cloud liquid water concent using the saturation adjustment method
    calc_liquid_water(ql->fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), bs.pref.data(), gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);

    // xz tmp slices for quantities which are used by multiple microphysics routines
    const int ikslice = gd.icells * gd.kcells;

    int tot_tmp_slice;
    int n_tmp_fld;
    int n_slice_in_tmp;
    if (per_slice)
        tot_tmp_slice = 12;
    else
        tot_tmp_slice = 8;
    int tot_tmp_fld = std::ceil(tot_tmp_slice/gd.jcells);
    std::shared_ptr<Field3d<TF>> atmp[tot_tmp_fld];
    for (int n=0; tot_tmp_fld; ++n)
        atmp[n] = fields.get_tmp();

    int n2d = 0;

    n2d++;
    n_tmp_fld = std::floor(n2d/gd.jcells);
    n_slice_in_tmp = n2d%gd.jcells * ikslice;
    TF* w_qr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
    n2d++;
    n_tmp_fld = std::floor(n2d/gd.jcells);
    n_slice_in_tmp = n2d%gd.jcells * ikslice;
    TF* w_nr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
    n2d++;
    n_tmp_fld = std::floor(n2d/gd.jcells);
    n_slice_in_tmp = n2d%gd.jcells * ikslice;
    TF* c_qr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
    n2d++;
    n_tmp_fld = std::floor(n2d/gd.jcells);
    n_slice_in_tmp = n2d%gd.jcells * ikslice;
    TF* c_nr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
    n2d++;
    n_tmp_fld = std::floor(n2d/gd.jcells);
    n_slice_in_tmp = n2d%gd.jcells * ikslice;
    TF* slope_qr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
    n2d++;
    n_tmp_fld = std::floor(n2d/gd.jcells);
    n_slice_in_tmp = n2d%gd.jcells * ikslice;
    TF* slope_nr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
    n2d++;
    n_tmp_fld = std::floor(n2d/gd.jcells);
    n_slice_in_tmp = n2d%gd.jcells * ikslice;
    TF* flux_qr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
    n2d++;
    n_tmp_fld = std::floor(n2d/gd.jcells);
    n_slice_in_tmp = n2d%gd.jcells * ikslice;
    TF* flux_nr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);

    // Autoconversion; formation of rain drop by coagulating cloud droplets
    mp::autoconversion(fields.st.at("qr")->fld.data(), fields.st.at("nr")->fld.data(), fields.st.at("qt")->fld.data(), fields.st.at("thl")->fld.data(),
                       fields.sp.at("qr")->fld.data(), ql->fld.data(), fields.rhoref.data(), bs.exnref.data(),
                       gd.istart, gd.jstart, gd.kstart,
                       gd.iend,   gd.jend,   gd.kend,
                       gd.icells, gd.ijcells);

    // Accretion; growth of raindrops collecting cloud droplets
    mp::accretion(fields.st.at("qr")->fld.data(), fields.st.at("qt")->fld.data(), fields.st.at("thl")->fld.data(),
                  fields.sp.at("qr")->fld.data(), ql->fld.data(), fields.rhoref.data(), bs.exnref.data(),
                  gd.istart, gd.jstart, gd.kstart,
                  gd.iend,   gd.jend,   gd.kend,
                  gd.icells, gd.ijcells);

    if(per_slice)
    {
        n2d++;
        n_tmp_fld = std::floor(n2d/gd.jcells);
        n_slice_in_tmp = n2d%gd.jcells * ikslice;
        TF* rain_mass = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
        n2d++;
        n_tmp_fld = std::floor(n2d/gd.jcells);
        n_slice_in_tmp = n2d%gd.jcells * ikslice;
        TF* rain_diam = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
        n2d++;
        n_tmp_fld = std::floor(n2d/gd.jcells);
        n_slice_in_tmp = n2d%gd.jcells * ikslice;
        TF* lambda_r  = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
        n2d++;
        n_tmp_fld = std::floor(n2d/gd.jcells);
        n_slice_in_tmp = n2d%gd.jcells * ikslice;
        TF* mu_r      = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);

        for (int j=gd.jstart; j<gd.jend; ++j)
        {
            mp2d::prepare_microphysics_slice(rain_mass, rain_diam, mu_r, lambda_r, fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(), fields.rhoref.data(),
                                             gd.istart, gd.iend, gd.kstart, gd.kend, gd.icells, gd.ijcells, j);

            // Evaporation; evaporation of rain drops in unsaturated environment
            mp2d::evaporation(fields.st.at("qr")->fld.data(), fields.st.at("nr")->fld.data(),  fields.st.at("qt")->fld.data(), fields.st.at("thl")->fld.data(),
                              fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),  ql->fld.data(),
                              fields.sp.at("qt")->fld.data(), fields.sp.at("thl")->fld.data(), fields.rhoref.data(), bs.exnref.data(), bs.pref.data(),
                              rain_mass, rain_diam,
                              gd.istart, gd.jstart, gd.kstart,
                              gd.iend,   gd.jend,   gd.kend,
                              gd.icells, gd.ijcells, j);

            // Self collection and breakup; growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
            mp2d::selfcollection_breakup(fields.st.at("nr")->fld.data(), fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(), fields.rhoref.data(),
                                         rain_mass, rain_diam, lambda_r,
                                         gd.istart, gd.jstart, gd.kstart,
                                         gd.iend,   gd.jend,   gd.kend,
                                         gd.icells, gd.ijcells, j);

            // Sedimentation; sub-grid sedimentation of rain
            mp2d::sedimentation_ss08(fields.st.at("qr")->fld.data(), fields.st.at("nr")->fld.data(),
                                     w_qr, w_nr, c_qr, c_nr, slope_qr, slope_nr, flux_qr, flux_nr, mu_r, lambda_r,
                                     fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),
                                     fields.rhoref.data(), gd.dzi.data(), gd.dz.data(), dt,
                                     gd.istart, gd.jstart, gd.kstart,
                                     gd.iend,   gd.jend,   gd.kend,
                                     gd.icells, gd.kcells, gd.ijcells, j);
        }
    }
    else
    {
        // Evaporation; evaporation of rain drops in unsaturated environment
        mp::evaporation(fields.st.at("qr")->fld.data(), fields.st.at("nr")->fld.data(),  fields.st.at("qt")->fld.data(), fields.st.at("thl")->fld.data(),
                        fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),  ql->fld.data(),
                        fields.sp.at("qt")->fld.data(), fields.sp.at("thl")->fld.data(), fields.rhoref.data(), bs.exnref.data(), bs.pref.data(),
                        gd.istart, gd.jstart, gd.kstart,
                        gd.iend,   gd.jend,   gd.kend,
                        gd.icells, gd.ijcells);

        // Self collection and breakup; growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
        mp::selfcollection_breakup(fields.st.at("nr")->fld.data(), fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(), fields.rhoref.data(),
                                   gd.istart, gd.jstart, gd.kstart,
                                   gd.iend,   gd.jend,   gd.kend,
                                   gd.icells, gd.ijcells);

        // Sedimentation; sub-grid sedimentation of rain
        mp::sedimentation_ss08(fields.st.at("qr")->fld.data(), fields.st.at("nr")->fld.data(), w_qr, w_nr, c_qr, c_nr, slope_qr, slope_nr, flux_qr, flux_nr,
                               fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),
                               fields.rhoref.data(), gd.dzi.data(), gd.dz.data(), dt,
                               gd.istart, gd.jstart, gd.kstart,
                               gd.iend,   gd.jend,   gd.kend,
                               gd.icells, gd.kcells, gd.ijcells);
    }

    for (int n=0;n<tot_tmp_fld;++n)
        fields.release_tmp(atmp[n]);
    fields.release_tmp(ql);
}

template<typename TF>
void Thermo_moist<TF>::get_mask(Field3d<TF>& mfield, Field3d<TF>& mfieldh, Stats<TF>& stats, std::string mask_name)
{
    auto& gd = grid.get_grid_data();
    #ifndef USECUDA
        bs_stats = bs;
    #endif

    if (mask_name == "ql")
    {
        auto ql = fields.get_tmp();
        auto qlh = fields.get_tmp();
        get_thermo_field(*ql,"ql", true, true);
        get_thermo_field(*qlh,"ql_h", true, true);
        stats.set_mask_true(mfield, mfieldh);
        stats.set_mask_thres(mfield, mfieldh, *ql,*qlh,0.,stats.Mask_type::Plus);
        stats.get_nmask(mfield, mfieldh);
        fields.release_tmp(ql);
        fields.release_tmp(qlh);
    }
    else if (mask_name == "qr")
    {
        TF threshold = 1e-8;
        const int wloc[] = {0,0,1};
        const int sloc[] = {0,0,0};

        auto qrh = fields.get_tmp();
        grid.interpolate_2nd(qrh->fld.data(),  fields.sp.at("qr")->fld.data(), sloc, wloc);
        stats.set_mask_true(mfield, mfieldh);
        stats.set_mask_thres(mfield, mfieldh, *fields.sp.at("qr"),*qrh,threshold,stats.Mask_type::Plus);
        stats.get_nmask(mfield, mfieldh);
        fields.release_tmp(qrh);
    }
    else if (mask_name == "qlcore")
    {
        stats.set_mask_true(mfield, mfieldh);

        auto ql = fields.get_tmp();
        auto qlh = fields.get_tmp();
        get_thermo_field(*ql,"ql", true, true);
        get_thermo_field(*qlh,"ql_h", true, true);
        stats.set_mask_thres(mfield, mfieldh, *ql,*qlh,0.,stats.Mask_type::Plus);
        fields.release_tmp(ql);
        fields.release_tmp(qlh);

        auto b = fields.get_tmp();
        auto bh = fields.get_tmp();
        get_thermo_field(*b,"b", true, true);
        get_thermo_field(*bh,"b_h", true, true);
        field3d_operators.calc_mean_profile(b->fld_mean.data(),b->fld.data());
        field3d_operators.calc_mean_profile(bh->fld_mean.data(),bh->fld.data());
        stats.set_mask_thres_pert(mfield, mfieldh, *b,*bh,0.,stats.Mask_type::Plus);
        stats.get_nmask(mfield, mfieldh);
        fields.release_tmp(b);
        fields.release_tmp(bh);
    }
    boundary_cyclic.exec(mfield.fld.data());
    boundary_cyclic.exec(mfieldh.fld.data());
    boundary_cyclic.exec_2d(mfieldh.fld_bot.data());

}


template<typename TF>
bool Thermo_moist<TF>::has_mask(std::string mask_name)
{
    if (std::find(available_masks.begin(), available_masks.end(), mask_name) != available_masks.end())
        return true;
    else
        return false;
}

template<typename TF>
bool Thermo_moist<TF>::check_field_exists(const std::string name)
{
    if (name == "b" || name == "ql" || name == "T")
        return true;
    else
        return false;
}

template<typename TF>
void Thermo_moist<TF>::update_time_dependent()
{
/*    if (swtimedep_pbot == 0)
        return;

    // Get/calculate the interpolation indexes/factors. Assing to zero to avoid compiler warnings.
    int index0 = 0, index1 = 0;
    TF fac0 = 0., fac1 = 0.;

    timeloop.get_interpolation_factors(index0, index1, fac0, fac1, timedep.time[it.first]);

    bs.pbot = fac0 * timedeppbot[index0] + fac1 * timedeppbot[index1];
*/}

template<typename TF>
void Thermo_moist<TF>::get_thermo_field(Field3d<TF>& fld, std::string name, bool cyclic, bool is_stat)
{
    auto& gd = grid.get_grid_data();
    background_state base;
    if (is_stat)
        base = bs_stats;
    else
        base = bs;

    // BvS: getThermoField() is called from subgrid-model, before thermo(), so re-calculate the hydrostatic pressure
    // Pass dummy as rhoref,bs.thvref to prevent overwriting base state
    if (bs.swupdatebasestate)
    {
        auto tmp = fields.get_tmp();
        calc_base_state(base.pref.data(), base.prefh.data(), &tmp->fld[0*gd.kcells], &tmp->fld[1*gd.kcells], &tmp->fld[2*gd.kcells], &tmp->fld[3*gd.kcells], base.exnref.data(), base.exnrefh.data(),
                fields.sp.at("thl")->fld_mean.data(), fields.sp.at("qt")->fld_mean.data(), base.pbot, gd.kstart, gd.kend, gd.z.data(), gd.dz.data(), gd.dzh.data());
        fields.release_tmp(tmp);
    }
    if (name == "b")
    {
        auto tmp = fields.get_tmp();
        calc_buoyancy(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.pref.data(), tmp->fld.data(), base.thvref.data(), gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kcells, gd.icells, gd.ijcells);
        fields.release_tmp(tmp);
    }
    if (name == "b_h")
    {
        auto tmp = fields.get_tmp();
        calc_buoyancy_h(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.prefh.data(), base.thvrefh.data(), &tmp->fld[0*gd.ijcells], &tmp->fld[1*gd.ijcells], &tmp->fld[2*gd.ijcells], gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        fields.release_tmp(tmp);
    }
    else if (name == "ql")
    {
        calc_liquid_water(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.pref.data(), gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
    }
    else if (name == "ql_h")
    {
        auto tmp = fields.get_tmp();
        calc_liquid_water_h(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.prefh.data(), &tmp->fld[0*gd.ijcells], &tmp->fld[1*gd.ijcells], gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        fields.release_tmp(tmp);
    }
    else if (name == "N2")
    {
        calc_N2(fld.fld.data(), fields.sp.at("thl")->fld.data(), gd.dzi.data(), base.thvref.data(), gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kcells, gd.icells, gd.ijcells);
    }
    else if (name == "T")
    {
        calc_T(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.pref.data(), base.exnref.data(), gd.istart, gd.iend, gd.jstart, gd.jend, gd.icells, gd.ijcells, gd.kcells);
    }
    else if (name == "T_h")
    {
        auto tmp = fields.get_tmp();
        calc_T_h(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.prefh.data(),&tmp->fld[0*gd.ijcells], &tmp->fld[1*gd.ijcells],
        &tmp->fld[2*gd.ijcells], gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        fields.release_tmp(tmp);
    }
    else
    {
        throw std::runtime_error("Illegal thermo field");
    }

    if (cyclic)
        boundary_cyclic.exec(fld.fld.data());
}

template<typename TF>
void Thermo_moist<TF>::get_buoyancy_surf(Field3d<TF>& b, bool is_stat)
{
    auto& gd = grid.get_grid_data();
    background_state base;
    if (is_stat)
        base = bs_stats;
    else
        base = bs;

    calc_buoyancy_bot(b.fld.data(), b.fld_bot.data(),
                      fields.sp.at("thl")->fld.data(), fields.sp.at("thl")->fld_bot.data(),
                      fields.sp.at("qt")->fld.data(), fields.sp.at("qt")->fld_bot.data(),
                      base.thvref.data(), base.thvrefh.data(), gd.icells, gd.jcells, gd.ijcells, gd.kstart);
    calc_buoyancy_fluxbot(b.flux_bot.data(), fields.sp.at("thl")->fld_bot.data(), fields.sp.at("thl")->flux_bot.data(), fields.sp.at("qt")->fld_bot.data(), fields.sp.at("qt")->flux_bot.data(), base.thvrefh.data(), gd.icells, gd.jcells, gd.kstart);
}

template<typename TF>
void Thermo_moist<TF>::get_buoyancy_fluxbot(Field3d<TF>& b, bool is_stat)
{
    auto& gd = grid.get_grid_data();
    background_state base;
    if (is_stat)
        base = bs_stats;
    else
        base = bs;

    calc_buoyancy_fluxbot(b.flux_bot.data(), fields.sp.at("thl")->fld_bot.data(), fields.sp.at("thl")->flux_bot.data(), fields.sp.at("qt")->fld_bot.data(), fields.sp.at("qt")->flux_bot.data(), base.thvrefh.data(), gd.icells, gd.jcells, gd.kstart);
}

template<typename TF>
void Thermo_moist<TF>::get_T_bot(Field3d<TF>& T_bot, bool is_stat)
{
    auto& gd = grid.get_grid_data();
    background_state base;
    if (is_stat)
        base = bs_stats;
    else
        base = bs;

    calc_T_bot(T_bot.fld_bot.data(), fields.sp.at("thl")->fld.data(), base.exnrefh.data(), base.thl0.data(),
               gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.icells, gd.ijcells);
}

template<typename TF>
void Thermo_moist<TF>::get_prog_vars(std::vector<std::string>& list)
{
    list.push_back("thl");
    list.push_back("qt");
}

template<typename TF>
TF Thermo_moist<TF>::get_buoyancy_diffusivity()
{
    // Use the diffusivity from the liquid water potential temperature
    return fields.sp.at("thl")->visc;
}

template<typename TF>
void Thermo_moist<TF>::create_stats(Stats<TF>& stats)
{
    // Add variables to the statistics
    if (stats.get_switch())
    {
        /* Add fixed base-state density and temperature profiles. Density should probably be in fields (?), but
           there the statistics are initialized before thermo->create() is called */
        stats.add_fixed_prof("rhoref",  "Full level basic state density", "kg m-3", "z",  fields.rhoref.data() );
        stats.add_fixed_prof("rhorefh", "Half level basic state density", "kg m-3", "zh", fields.rhorefh.data());
        stats.add_fixed_prof("thvref",  "Full level basic state virtual potential temperature", "K", "z", bs.thvref.data() );
        stats.add_fixed_prof("thvrefh", "Half level basic state virtual potential temperature", "K", "zh",bs.thvrefh.data());

        if (bs_stats.swupdatebasestate)
        {
            stats.add_prof("ph",   "Full level hydrostatic pressure", "Pa",     "z" );
            stats.add_prof("phh",  "Half level hydrostatic pressure", "Pa",     "zh");
            stats.add_prof("rho",  "Full level density",  "kg m-3", "z" );
            stats.add_prof("rhoh", "Half level density",  "kg m-3", "zh");
        }
        else
        {
            stats.add_fixed_prof("ph",  "Full level hydrostatic pressure", "Pa", "z",  bs.pref.data() );
            stats.add_fixed_prof("phh", "Half level hydrostatic pressure", "Pa", "zh", bs.prefh.data());
        }

        stats.add_prof("b", "Buoyancy", "m s-2", "z");
        for (int n=2; n<5; ++n)
        {
            std::stringstream ss;
            ss << n;
            std::string sn = ss.str();
            stats.add_prof("b"+sn, "Moment " +sn+" of the buoyancy", "(m s-2)"+sn,"z");
        }

        stats.add_prof("bgrad", "Gradient of the buoyancy", "m s-3", "zh");
        stats.add_prof("bw"   , "Turbulent flux of the buoyancy", "m2 s-3", "zh");
        stats.add_prof("bdiff", "Diffusive flux of the buoyancy", "m2 s-3", "zh");
        stats.add_prof("bflux", "Total flux of the buoyancy", "m2 s-3", "zh");

        stats.add_prof("ql", "Liquid water mixing ratio", "kg kg-1", "z");
        stats.add_prof("cfrac", "Cloud fraction", "-","z");

        stats.add_time_series("lwp", "Liquid water path", "kg m-2");
        stats.add_time_series("ccover", "Projected cloud cover", "-");

        // BvS:micro
        if(swmicro == Micro_type::two_mom_warm)
        {
            stats.add_time_series("rwp", "Rain water path", "kg m-2");

            if(swmicrobudget)
            {
                stats.add_prof("auto_qrt" , "Autoconversion tendency qr", "kg kg-1 s-1", "z");
                stats.add_prof("auto_nrt" , "Autoconversion tendency nr", "m-3 s-1", "z");
                stats.add_prof("auto_thlt", "Autoconversion tendency thl", "K s-1", "z");
                stats.add_prof("auto_qtt" , "Autoconversion tendency qt", "kg kg-1 s-1", "z");

                stats.add_prof("evap_qrt" , "Evaporation tendency qr", "kg kg-1 s-1", "z");
                stats.add_prof("evap_nrt" , "Evaporation tendency nr", "m-3 s-1", "z");
                stats.add_prof("evap_thlt", "Evaporation tendency thl", "K s-1", "z");
                stats.add_prof("evap_qtt" , "Evaporation tendency qt", "kg kg-1 s-1", "z");

                stats.add_prof("accr_qrt" , "Accretion tendency qr", "kg kg-1 s-1", "z");
                stats.add_prof("accr_thlt", "Accretion tendency thl", "K s-1", "z");
                stats.add_prof("accr_qtt" , "Accretion tendency qt", "kg kg-1 s-1", "z");

                stats.add_prof("scbr_nrt" , "Selfcollection and breakup tendency nr", "m-3 s-1", "z");

                stats.add_prof("sed_qrt"  , "Sedimentation tendency qr", "kg kg-1 s-1", "z");
                stats.add_prof("sed_nrt"  , "Sedimentation tendency nr", "m-3 s-1", "z");
            }
        }
    }
}



template<typename TF>
void Thermo_moist<TF>::create_column(Column<TF>& column)
{
    // add the profiles to the columns
    if (column.get_switch())
        column.add_prof("b", "Buoyancy", "m s-2", "z");
        column.add_prof("ql", "Liquid water mixing ratio", "kg kg-1", "z");
}
template<typename TF>
void Thermo_moist<TF>::create_cross(Cross<TF>& cross)
{
    if (cross.get_switch())
    {
        swcross_b = false;
        swcross_ql = false;
        std::vector<std::string> allowed_crossvars_b; ///< List with allowed cross variables
        std::vector<std::string> allowed_crossvars_ql; ///< List with allowed cross variables

        allowed_crossvars_b.push_back("b");
        allowed_crossvars_b.push_back("bbot");
        allowed_crossvars_b.push_back("bfluxbot");
        if (grid.swspatialorder == "4")
            allowed_crossvars_b.push_back("blngrad");
        allowed_crossvars_ql.push_back("ql");
        allowed_crossvars_ql.push_back("qlpath");
        allowed_crossvars_ql.push_back("qlbase");
        allowed_crossvars_ql.push_back("qltop");

        // Get global cross-list from cross.cxx
        std::vector<std::string> *crosslist_global = cross.get_crosslist();

        // Check input list of cross variables (crosslist)
        std::vector<std::string>::iterator it=crosslist_global->begin();
        while (it != crosslist_global->end())
        {
            if (std::count(allowed_crossvars_b.begin(),allowed_crossvars_b.end(),*it))
            {
                // Remove variable from global list, put in local list
                crosslist.push_back(*it);
                crosslist_global->erase(it); // erase() returns iterator of next element..
                swcross_b = true;
            }
            else if (std::count(allowed_crossvars_ql.begin(),allowed_crossvars_ql.end(),*it))
            {
                // Remove variable from global list, put in local list
                crosslist.push_back(*it);
                crosslist_global->erase(it); // erase() returns iterator of next element..
                swcross_ql = true;
            }
            else
                ++it;
        }

        // Sort crosslist to group ql and b variables
        std::sort(crosslist.begin(),crosslist.end());
    }
}

template<typename TF>
void Thermo_moist<TF>::create_dump(Dump<TF>& dump)
{
    if (dump.get_switch())
    {
        // Get global cross-list from cross.cxx
        std::vector<std::string> *dumplist_global = dump.get_dumplist();

        // Check if fields in dumplist are retrievable thermo fields
        std::vector<std::string>::iterator dumpvar=dumplist_global->begin();
        while (dumpvar != dumplist_global->end())
        {
            if (check_field_exists(*dumpvar))
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
void Thermo_moist<TF>::exec_stats(Stats<TF>& stats, std::string mask_name, Field3d<TF>& mask_field, Field3d<TF>& mask_fieldh,
        const Diff<TF>& diff, const double dt)
{
    auto& gd = grid.get_grid_data();
    #ifndef USECUDA
        bs_stats = bs;
    #endif

    Mask<TF>& m = stats.masks[mask_name];

    const TF no_offset = 0.;

    // calculate the buoyancy and its surface flux for the profiles
    auto b = fields.get_tmp();
    get_thermo_field(*b, "b",true, true);
    get_buoyancy_surf(*b, true);
    get_buoyancy_fluxbot(*b, true);


    // define the location
    const int sloc[] = {0,0,0};

    // calculate the mean
    stats.calc_mean(m.profs["b"].data.data(), b->fld.data(), no_offset, mask_field.fld.data(), stats.nmask.data());

    // calculate the moments
    for (int n=2; n<5; ++n)
    {
        std::string sn = std::to_string(n);
        stats.calc_moment(b->fld.data(), m.profs["b"].data.data(), m.profs["b"+sn].data.data(), n, mask_field.fld.data(), stats.nmask.data());
    }

    if (grid.swspatialorder == "2")
    {
        auto tmp = fields.get_tmp();
        stats.calc_grad_2nd(b->fld.data(), m.profs["bgrad"].data.data(), gd.dzhi.data(), mask_fieldh.fld.data(), stats.nmaskh.data());
        stats.calc_flux_2nd(b->fld.data(), m.profs["b"].data.data(), fields.mp["w"]->fld.data(), m.profs["w"].data.data(),
                            m.profs["bw"].data.data(), tmp->fld.data(), sloc, mask_fieldh.fld.data(), stats.nmaskh.data());
        if (diff.get_switch() == Diffusion_type::Diff_smag2)
        {
            stats.calc_diff_2nd(b->fld.data(), fields.mp["w"]->fld.data(), fields.sd["evisc"]->fld.data(),
                                m.profs["bdiff"].data.data(), gd.dzhi.data(),
                                b->flux_bot.data(), b->flux_top.data(), diff.tPr, sloc,
                                mask_fieldh.fld.data(), stats.nmaskh.data());
        }
        else
            stats.calc_diff_2nd(b->fld.data(), m.profs["bdiff"].data.data(), gd.dzhi.data(),
                                get_buoyancy_diffusivity(), sloc, mask_fieldh.fld.data(), stats.nmaskh.data());
        fields.release_tmp(tmp);
    }
    // calculate the total fluxes
    stats.add_fluxes(m.profs["bflux"].data.data(), m.profs["bw"].data.data(), m.profs["bdiff"].data.data());

    // calculate the sorted buoyancy profile
    //stats->calc_sorted_prof(fields.sd["tmp1"]->data, fields.sd["tmp2"]->data, m->profs["bsort"].data);
    fields.release_tmp(b);


    // calculate the liquid water stats
    auto ql = fields.get_tmp();
    get_thermo_field(*ql,"ql", true, true);
    stats.calc_mean(m.profs["ql"].data.data(), ql->fld.data(), no_offset, mask_field.fld.data(), stats.nmask.data());
    //stats.calc_count(m.profs["ccover"].data.data(), ql->fld.data(), no_offset, mask_field.fld.data(), stats.nmask.data());
    //stats.calc_cover(m.profs["cfrac"].data.data(), ql->fld.data(), no_offset, mask_field.fld.data(), stats.nmask.data());
    //stats.calc_path(m.profs["lwp"].data.data(), ql->fld.data(), no_offset, mask_field.fld.data(), stats.nmask.data());

    // BvS:micro
    if(swmicro == Micro_type::two_mom_warm)
    {
        //stats.calc_path(m.profs["rwp"].data.data(), fields.sp.at("qr")->fld.data(), no_offset, mask_field.fld.data(), stats.nmask.data());

        if(swmicrobudget)
        {
            auto qrt  = fields.get_tmp();
            auto nrt  = fields.get_tmp();
            auto qtt  = fields.get_tmp();
            auto thlt = fields.get_tmp();

            // Autoconversion
            mp::zero(qrt->fld.data(), gd.ncells);
            mp::zero(nrt->fld.data(), gd.ncells);
            mp::zero(qtt->fld.data(), gd.ncells);
            mp::zero(thlt->fld.data(), gd.ncells);

            mp::autoconversion(qrt->fld.data(), nrt->fld.data(), qtt->fld.data(), thlt->fld.data(),
                               fields.sp.at("qr")->fld.data(), ql->fld.data(), fields.rhoref.data(), bs_stats.exnref.data(),
                               gd.istart, gd.jstart, gd.kstart,
                               gd.iend,   gd.jend,   gd.kend,
                               gd.icells, gd.ijcells);

            stats.calc_mean(m.profs["auto_qrt" ].data.data(), qrt->fld.data(), no_offset, mask_field.fld.data(),  stats.nmask.data());
            stats.calc_mean(m.profs["auto_nrt" ].data.data(), nrt->fld.data(), no_offset, mask_field.fld.data(),  stats.nmask.data());
            stats.calc_mean(m.profs["auto_qtt" ].data.data(), qtt->fld.data(), no_offset, mask_field.fld.data(),  stats.nmask.data());
            stats.calc_mean(m.profs["auto_thlt"].data.data(), thlt->fld.data(), no_offset, mask_field.fld.data(),  stats.nmask.data());

            // Evaporation
            mp::zero(qrt->fld.data(), gd.ncells);
            mp::zero(nrt->fld.data(), gd.ncells);
            mp::zero(qtt->fld.data(), gd.ncells);
            mp::zero(thlt->fld.data(), gd.ncells);

            mp::evaporation(qrt->fld.data(), nrt->fld.data(),  qtt->fld.data(), thlt->fld.data(),
                            fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),  ql->fld.data(),
                            fields.sp.at("qt")->fld.data(), fields.sp.at("thl")->fld.data(), fields.rhoref.data(), bs_stats.exnref.data(), bs_stats.pref.data(),
                            gd.istart, gd.jstart, gd.kstart,
                            gd.iend,   gd.jend,   gd.kend,
                            gd.icells, gd.ijcells);

            stats.calc_mean(m.profs["evap_qrt" ].data.data(), qrt->fld.data(), no_offset, mask_field.fld.data(),  stats.nmask.data());
            stats.calc_mean(m.profs["evap_nrt" ].data.data(), nrt->fld.data(), no_offset, mask_field.fld.data(),  stats.nmask.data());
            stats.calc_mean(m.profs["evap_qtt" ].data.data(), qtt->fld.data(), no_offset, mask_field.fld.data(),  stats.nmask.data());
            stats.calc_mean(m.profs["evap_thlt"].data.data(), thlt->fld.data(), no_offset, mask_field.fld.data(),  stats.nmask.data());

            // Accretion
            mp::zero(qrt->fld.data(), gd.ncells);
            mp::zero(nrt->fld.data(), gd.ncells);
            mp::zero(qtt->fld.data(), gd.ncells);

            mp::accretion(qrt->fld.data(), nrt->fld.data(), qtt->fld.data(),
                          fields.sp.at("qr")->fld.data(), ql->fld.data(), fields.rhoref.data(), bs_stats.exnref.data(),
                          gd.istart, gd.jstart, gd.kstart,
                          gd.iend,   gd.jend,   gd.kend,
                          gd.icells, gd.ijcells);

            stats.calc_mean(m.profs["accr_qrt" ].data.data(), qrt->fld.data(), no_offset, mask_field.fld.data(),  stats.nmask.data());
            stats.calc_mean(m.profs["accr_qtt" ].data.data(), nrt->fld.data(), no_offset, mask_field.fld.data(),  stats.nmask.data());
            stats.calc_mean(m.profs["accr_thlt"].data.data(), qtt->fld.data(), no_offset, mask_field.fld.data(),  stats.nmask.data());

            // Selfcollection and breakup
            mp::zero(qrt->fld.data(), gd.ncells);

            mp::selfcollection_breakup(qrt->fld.data(), fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(), fields.rhoref.data(),
                                       gd.istart, gd.jstart, gd.kstart,
                                       gd.iend,   gd.jend,   gd.kend,
                                       gd.icells, gd.ijcells);

            stats.calc_mean(m.profs["scbr_nrt" ].data.data(), qrt->fld.data(), no_offset, mask_field.fld.data(),  stats.nmask.data());

            // Sedimentation
            mp::zero(qrt->fld.data(), gd.ncells);
            mp::zero(nrt->fld.data(), gd.ncells);
            mp::zero(qtt->fld.data(), gd.ncells);
            mp::zero(thlt->fld.data(), gd.ncells);

            // 1. Get number of substeps based on sedimentation with CFL=1
            //const TF dt = model->timeloop->get_sub_time_step();
            //int nsubstep = mp::get_sedimentation_steps(fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),
            //                                           fields.rhoref, gd.dzi, dt,
            //                                           gd.istart, gd.jstart, gd.kstart,
            //                                           gd.iend,   gd.jend,   gd.kend,
            //                                           gd.icells, gd.ijcells);

            //// 2. Synchronize over all MPI processes
            //gd.get_max(&nsubstep);

            //// Stay a bit further from CFL=1
            //nsubstep *= 2;

            //// Sedimentation in nsubstep steps:
            //mp::sedimentation_sub(qrt->fld.data(), nrt->fld.data(),
            //                      qtt->fld.data(), thlt->fld.data(),
            //                      fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),
            //                      fields.rhoref, gd.dzi, gd.dzhi, dt,
            //                      gd.istart, gd.jstart, gd.kstart,
            //                      gd.iend,   gd.jend,   gd.kend,
            //                      gd.icells, gd.kcells, gd.ijcells,
            //                      nsubstep);
            const int ikslice = gd.icells * gd.kcells;

            int tot_tmp_slice;
            int n_tmp_fld;
            int n_slice_in_tmp;
            tot_tmp_slice = 8;
            int tot_tmp_fld = std::ceil(tot_tmp_slice/gd.jcells);
            std::shared_ptr<Field3d<TF>> atmp[tot_tmp_fld];
            for (int n=0; tot_tmp_fld; ++n)
                atmp[n] = fields.get_tmp();

            int n2d = 0;

            n2d++;
            n_tmp_fld = std::floor(n2d/gd.jcells);
            n_slice_in_tmp = n2d%gd.jcells * ikslice;
            TF* w_qr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
            n2d++;
            n_tmp_fld = std::floor(n2d/gd.jcells);
            n_slice_in_tmp = n2d%gd.jcells * ikslice;
            TF* w_nr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
            n2d++;
            n_tmp_fld = std::floor(n2d/gd.jcells);
            n_slice_in_tmp = n2d%gd.jcells * ikslice;
            TF* c_qr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
            n2d++;
            n_tmp_fld = std::floor(n2d/gd.jcells);
            n_slice_in_tmp = n2d%gd.jcells * ikslice;
            TF* c_nr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
            n2d++;
            n_tmp_fld = std::floor(n2d/gd.jcells);
            n_slice_in_tmp = n2d%gd.jcells * ikslice;
            TF* slope_qr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
            n2d++;
            n_tmp_fld = std::floor(n2d/gd.jcells);
            n_slice_in_tmp = n2d%gd.jcells * ikslice;
            TF* slope_nr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
            n2d++;
            n_tmp_fld = std::floor(n2d/gd.jcells);
            n_slice_in_tmp = n2d%gd.jcells * ikslice;
            TF* flux_qr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);
            n2d++;
            n_tmp_fld = std::floor(n2d/gd.jcells);
            n_slice_in_tmp = n2d%gd.jcells * ikslice;
            TF* flux_nr = &(atmp[n_tmp_fld]->fld[n_slice_in_tmp]);

            // Sedimentation; sub-grid sedimentation of rain
            mp::sedimentation_ss08(qrt->fld.data(), nrt->fld.data(), w_qr, w_nr, c_qr, c_nr, slope_qr, slope_nr, flux_qr, flux_nr,
                                  fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),
                                  fields.rhoref.data(), gd.dzi.data(), gd.dz.data(), dt,
                                  gd.istart, gd.jstart, gd.kstart,
                                  gd.iend,   gd.jend,   gd.kend,
                                  gd.icells, gd.kcells, gd.ijcells);

            for (int n=0;n<tot_tmp_fld;++n)
                fields.release_tmp(atmp[n]);

            stats.calc_mean(m.profs["sed_qrt"].data.data(), qrt->fld.data(), no_offset, mask_field.fld.data(),  stats.nmask.data());
            stats.calc_mean(m.profs["sed_nrt"].data.data(), nrt->fld.data(), no_offset, mask_field.fld.data(),  stats.nmask.data());
        }
    }

    fields.release_tmp(ql);

    // Calculate base state in tmp array
    if (bs_stats.swupdatebasestate)
    {
            m.profs["ph"  ].data = bs_stats.pref;
            m.profs["phh" ].data = bs_stats.prefh;
            m.profs["rho" ].data = fields.rhoref;
            m.profs["rhoh"].data = fields.rhorefh;
    }
}


template<typename TF>
void Thermo_moist<TF>::exec_column(Column<TF>& column)
{
    auto& gd = grid.get_grid_data();
    #ifndef USECUDA
        bs_stats = bs;
    #endif
    const TF no_offset = 0.;
    auto output = fields.get_tmp();

    for (auto& it : dumplist)
    {
        if (it == "b")
            get_thermo_field(*output, "b",false, true);
        else if (it == "ql")
            get_thermo_field(*output, "ql",false, true);
        else if (it == "T")
            get_thermo_field(*output, "T",false, true);
        else
        {
            master.print_error("Thermo dump of field \"%s\" not supported\n",it.c_str());
            throw std::runtime_error("Error in Thermo Dump");
        }
        column.calc_column(it, output->fld.data(), no_offset);
    }
    fields.release_tmp(output);
}


template<typename TF>
void Thermo_moist<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    auto& gd = grid.get_grid_data();
    #ifndef USECUDA
        bs_stats = bs;
    #endif
    auto output = fields.get_tmp();

    if(swcross_b)
    {
        get_thermo_field(*output, "b",false, true);
        get_buoyancy_fluxbot(*output, true);
    }
    for (auto& it : crosslist)
    {
        if (it == "b")
            cross.cross_simple(output->fld.data(), "b", iotime);
        else if (it == "blngrad")
            cross.cross_lngrad(output->fld.data(), "blngrad", iotime);
        else if (it == "bbot")
            cross.cross_plane(output->fld_bot.data(), "bbot", iotime);
        else if (it == "bfluxbot")
            cross.cross_plane(output->flux_bot.data(), "bfluxbot", iotime);
    }

    if(swcross_ql)
    {
        get_thermo_field(*output, "ql",false, true);
    }
    for (auto& it : crosslist)
    {
        if (it == "ql")
            cross.cross_simple(output->fld.data(), "ql", iotime);
        if (it == "qlpath")
            cross.cross_path(output->fld.data(), "ql", iotime);
        if (it == "qlbase")
            cross.cross_height_threshold(output->fld.data(), 0., cross.Direction::Bottom_to_top, "ql", iotime);
        if (it == "qltop")
            cross.cross_height_threshold(output->fld.data(), 0., cross.Direction::Top_to_bottom, "ql", iotime);
    }

    fields.release_tmp(output);
}

template<typename TF>
void Thermo_moist<TF>::exec_dump(Dump<TF>& dump, unsigned long iotime)
{
    #ifndef USECUDA
        bs_stats = bs;
    #endif
    auto output = fields.get_tmp();

    for (auto& it : dumplist)
    {
        if (it == "b")
            get_thermo_field(*output, "b",false, true);
        else if (it == "ql")
            get_thermo_field(*output, "ql",false, true);
        else if (it == "T")
            get_thermo_field(*output, "T",false, true);
        else
        {
            master.print_error("Thermo dump of field \"%s\" not supported\n",it.c_str());
            throw std::runtime_error("Error in Thermo Dump");
        }
        dump.save_dump(output->fld.data(), it, iotime);
    }
    fields.release_tmp(output);
}


template class Thermo_moist<double>;
template class Thermo_moist<float>;
