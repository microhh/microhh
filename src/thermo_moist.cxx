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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <netcdf>
#include "grid.h"
#include "fields.h"
#include "thermo_moist.h"
#include "diff.h"
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
using namespace Constants;
using namespace Thermo_moist_functions;

namespace
{

   template<typename TF>
   void calc_top_and_bot(TF* restrict thl0, TF* restrict qt0,
                         const TF* const z, const TF* const zh,
                         const TF* const dzhi,
                         const int kstart, const int kend)
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
                               TF* restrict ql, TF* restrict thvrefh,
                               const int istart, const int iend,
                               const int jstart, const int jend,
                               const int kstart, const int kend,
                               const int jj, const int kk)
   {
       #pragma omp parallel for
       for (int k=kstart+1; k<kend; k++)
       {
           const TF exnh = exner(ph[k]);
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
                      TF* restrict p, TF* restrict ql, TF* restrict thvref,
                      const int istart, const int iend,
                      const int jstart, const int jend,
                      const int kstart, const int kend,
                      const int kcells, const int jj, const int kk)
   {
       #pragma omp parallel for
       for (int k=0; k<kcells; k++)
       {
           const TF ex = exner(p[k]);
           if (k >= kstart && k < kend)
           {
               for (int j=jstart; j<jend; j++)
                   #pragma ivdep
                   for (int i=istart; i<iend; i++)
                   {
                       const int ijk = i + j*jj + k*kk;
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
                   b[ijk] = buoyancy(ex, thl[ijk], qt[ijk], ql[ijk], thvref[k]);
               }
       }
   }

  template<typename TF>
  void calc_buoyancy_h(TF* restrict bh, TF* restrict thl, TF* restrict qt,
                       TF* restrict ph, TF* restrict thvrefh, TF* restrict thlh, TF* restrict qth, TF* restrict ql,
                       const int istart, const int iend,
                       const int jstart, const int jend,
                       const int kstart, const int kend,
                       const int jj, const int kk)
    {
        using Finite_difference::O2::interp2;

        for (int k=kstart; k<kend+1; k++)
        {
            const TF exnh = exner(ph[k]);

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

                    bh[ijk] = buoyancy(exnh, thlh[ij], qth[ij], ql[ij], thvrefh[k]);
                }
        }
   }

   template<typename TF>
   void calc_liquid_water(TF* restrict ql, TF* restrict thl, TF* restrict qt, TF* restrict p,
                          const int istart, const int iend,
                          const int jstart, const int jend,
                          const int kstart, const int kend,
                          const int jj, const int kk)
   {
       // Calculate the ql field
       #pragma omp parallel for
       for (int k=kstart; k<kend; k++)
       {
           const TF ex = exner(p[k]);
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
                           TF* restrict ph, TF* restrict thlh, TF* restrict qth,
                           const int istart, const int iend,
                           const int jstart, const int jend,
                           const int kstart, const int kend,
                           const int jj, const int kk)
   {
       using Finite_difference::O2::interp2;

       for (int k=kstart+1; k<kend+1; k++)
       {
           const TF exnh = exner(ph[k]);

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
   void calc_N2(TF* restrict N2, const TF* const restrict thl, const TF* const restrict dzi, TF* restrict thvref,
                const int istart, const int iend,
                const int jstart, const int jend,
                const int kstart, const int kend,
                const int jj, const int kk)
   {
       #pragma omp parallel for
       for (int k=kstart; k<kend; ++k)
           for (int j=jstart; j<jend; ++j)
               #pragma ivdep
               for (int i=istart; i<iend; ++i)
               {
                   const int ijk = i + j*jj + k*kk;
                   N2[ijk] = grav<TF>/thvref[k]*TF(0.5)*(thl[ijk+kk] - thl[ijk-kk])*dzi[k];
               }
   }

   template<typename TF>
   void calc_T(TF* const restrict T, const TF* const restrict thl, const TF* const restrict qt,
               const TF* const restrict pref, const TF* const restrict exnref,
               const int istart, const int iend,
               const int jstart, const int jend,
               const int kstart, const int kend,
               const int jj, const int kk)
   {
       #pragma omp parallel for
       for (int k=kstart; k<kend; ++k)
       {
           for (int j=jstart; j<jend; ++j)
               #pragma ivdep
               for (int i=istart; i<iend; ++i)
               {
                   const int ijk = i + j*jj+ k*kk;
                   T[ijk] = sat_adjust(thl[ijk], qt[ijk], pref[k], exnref[k]).t;
               }
       }
   }

   template<typename TF>
   void calc_T_h(TF* restrict Th, TF* restrict thl,  TF* restrict qt,
                 TF* restrict ph, TF* restrict thlh, TF* restrict qth, TF* restrict ql,
                 const int istart, const int iend,
                 const int jstart, const int jend,
                 const int kstart, const int kend,
                 const int jj, const int kk)
     {
           using Finite_difference::O2::interp2;

           for (int k=kstart+1; k<kend+1; k++)
           {
               const TF exnh = exner(ph[k]);
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
                          TF* restrict thvref, TF* restrict thvrefh,
                          const int icells, const int jcells,
                          const int ijcells, const int kstart)
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
   void calc_buoyancy_fluxbot(TF* restrict bfluxbot, TF* restrict thlbot, TF* restrict thlfluxbot,
                              TF* restrict qtbot, TF* restrict qtfluxbot, TF* restrict thvrefh,
                              const int icells, const int jcells, const int kstart)
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
                               TF* restrict ql, TF* restrict thvrefh,
                               const int istart, const int iend,
                               const int jstart, const int jend,
                               const int kstart, const int kend,
                               const int icells, const int ijcells)
   {
       const int jj  = icells;
       const int kk1 = 1*ijcells;
       const int kk2 = 2*ijcells;

       for (int k=kstart+1; k<kend; k++)
       {
           const TF exnh = exner(ph[k]);
           for (int j=jstart; j<jend; j++)
               #pragma ivdep
               for (int i=istart; i<iend; i++)
               {
                   const int ijk = i + j*jj + k*kk1;
                   const int ij  = i + j*jj;

                   thlh[ij]    = interp4c(thl[ijk-kk2], thl[ijk-kk1], thl[ijk], thl[ijk+kk1]);
                   qth[ij]     = interp4c(qt[ijk-kk2],  qt[ijk-kk1],  qt[ijk],  qt[ijk+kk1]);
                   const TF tl = thlh[ij] * exnh;

                   // Calculate first estimate of ql using Tl
                   // if ql(Tl)>0, saturation adjustment routine needed
                   ql[ij]  = qth[ij]-qsat(ph[k], tl);
               }

           for (int j=jstart; j<jend; j++)
               #pragma ivdep
               for (int i=istart; i<iend; i++)
               {
                   const int ij = i + j*jj;

                   if (ql[ij] > 0)   // already doesn't vectorize because of iteration in sat_adjust()
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
}


template<typename TF>
Thermo_moist<TF>::Thermo_moist(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Thermo<TF>(masterin, gridin, fieldsin, inputin),
    boundary_cyclic(masterin, gridin),
    field3d_operators(master, grid, fieldsin)
{
    auto& gd = grid.get_grid_data();

    swthermo = "moist";

    // 4th order code is not implemented in thermo_moist
    if (grid.get_spatial_order() == Grid_order::Fourth)
        throw std::runtime_error("swthermo=moist is not supported for swspatialorder=4\n");

    // Initialize the prognostic fields
    fields.init_prognostic_field("thl", "Liquid water potential temperature", "K", gd.sloc);
    fields.init_prognostic_field("qt", "Total water mixing ratio", "kg kg-1", gd.sloc);

    // Get the diffusivities of temperature and moisture
    fields.sp.at("thl")->visc = inputin.get_item<TF>("fields", "svisc", "thl");
    fields.sp.at("qt")->visc = inputin.get_item<TF>("fields", "svisc", "qt");

    // Test if the diffusivities of theta and qt are equal, else throw error
    if (fields.sp.at("thl")->visc != fields.sp.at("qt")->visc)
        throw std::runtime_error("The diffusivities of temperature and moisture must be equal\n");

    bs.pbot = inputin.get_item<TF>("thermo", "pbot", "");

    // Get base state option (boussinesq or anelastic)
    std::string swbasestate_in = inputin.get_item<std::string>("thermo", "swbasestate", "", "");
    if (swbasestate_in == "boussinesq")
        bs.swbasestate = Basestate_type::boussinesq;
    else if (swbasestate_in == "anelastic")
        bs.swbasestate = Basestate_type::anelastic;
    else
        throw std::runtime_error("Invalid option for \"swbasestate\"");


    // BvS test for updating hydrostatic prssure during run
    // swupdate..=0 -> initial base state pressure used in saturation calculation
    // swupdate..=1 -> base state pressure updated before saturation calculation
    bs.swupdatebasestate = inputin.get_item<bool>("thermo", "swupdatebasestate", "", false);

    // Time variable surface pressure
    tdep_pbot = std::make_unique<Timedep<TF>>(master, grid, "p_sbot", inputin.get_item<bool>("thermo", "swtimedep_pbot", "", false));

    available_masks.insert(available_masks.end(), {"ql", "qlcore"});
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
    bs.rhoref.resize(gd.kcells);
    bs.rhorefh.resize(gd.kcells);
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
    data_block.get_vector(bs.thl0, "thl", gd.ktot, 0, gd.kstart);
    data_block.get_vector(bs.qt0, "qt", gd.ktot, 0, gd.kstart);

    calc_top_and_bot(bs.thl0.data(), bs.qt0.data(), gd.z.data(), gd.zh.data(), gd.dzhi.data(), gd.kstart, gd.kend);

    // 4. Calculate the initial/reference base state
    calc_base_state(bs.pref.data(), bs.prefh.data(), bs.rhoref.data(), bs.rhorefh.data(), bs.thvref.data(),
                    bs.thvrefh.data(), bs.exnref.data(), bs.exnrefh.data(), bs.thl0.data(), bs.qt0.data(), bs.pbot,
                    gd.kstart, gd.kend, gd.z.data(), gd.dz.data(), gd.dzh.data());

    // 5. In Boussinesq mode, overwrite reference temperature and density
    if (bs.swbasestate == Basestate_type::boussinesq)
    {
        bs.thvref0 = inputin.get_item<TF>("thermo", "thvref0", "");

        for (int k=0; k<gd.kcells; ++k)
        {
            bs.rhoref[k]  = 1.;
            bs.rhorefh[k] = 1.;
            bs.thvref[k]  = bs.thvref0;
            bs.thvrefh[k] = bs.thvref0;
        }
    }

    // 6. Copy the initial reference to the fields. This is the reference used in the dynamics.
    //    This one is not updated throughout the simulation to be consistent with the anelastic approximation.
    fields.rhoref = bs.rhoref;
    fields.rhorefh = bs.rhorefh;

    // 7. Process the time dependent surface pressure
    tdep_pbot->create_timedep();


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

    // Re-calculate hydrostatic pressure and exner, pass dummy as thvref to prevent overwriting base state
    auto tmp = fields.get_tmp();
    if (bs.swupdatebasestate)
    {
        calc_base_state(bs.pref.data(), bs.prefh.data(),
                        bs.rhoref.data(), bs.rhorefh.data(), &tmp->fld[0*gd.kcells], &tmp->fld[1*gd.kcells],
                        bs.exnref.data(), bs.exnrefh.data(), fields.sp.at("thl")->fld_mean.data(), fields.sp.at("qt")->fld_mean.data(),
                        bs.pbot, gd.kstart, gd.kend, gd.z.data(), gd.dz.data(), gd.dzh.data());
    }

    // extend later for gravity vector not normal to surface
    calc_buoyancy_tend_2nd(fields.mt.at("w")->fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), bs.prefh.data(),
                           &tmp->fld[0*gd.ijcells], &tmp->fld[1*gd.ijcells],
                           &tmp->fld[2*gd.ijcells], bs.thvrefh.data(),
                           gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                           gd.icells, gd.ijcells);

    fields.release_tmp(tmp);
}
#endif

template<typename TF>
unsigned long Thermo_moist<TF>::get_time_limit(unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}

template<typename TF>
void Thermo_moist<TF>::get_mask(Stats<TF>& stats, std::string mask_name)
{
    #ifndef USECUDA
    bs_stats = bs;
    #endif

    if (mask_name == "ql")
    {
        auto ql = fields.get_tmp();
        auto qlh = fields.get_tmp();

        get_thermo_field(*ql, "ql", true, true);
        get_thermo_field(*qlh, "ql_h", true, true);

        stats.set_mask_thres(mask_name, *ql, *qlh, 0., Stats_mask_type::Plus);

        fields.release_tmp(ql);
        fields.release_tmp(qlh);
    }
    else if (mask_name == "qlcore")
    {
        auto ql = fields.get_tmp();
        auto qlh = fields.get_tmp();

        get_thermo_field(*ql, "ql", true, true);
        get_thermo_field(*qlh, "ql_h", true, true);

        stats.set_mask_thres(mask_name, *ql, *qlh, 0., Stats_mask_type::Plus);

        fields.release_tmp(ql);
        fields.release_tmp(qlh);

        auto b = fields.get_tmp();
        auto bh = fields.get_tmp();

        get_thermo_field(*b, "b", true, true);
        get_thermo_field(*bh, "b_h", true, true);

        field3d_operators.calc_mean_profile(b->fld_mean.data(), b->fld.data());
        field3d_operators.calc_mean_profile(bh->fld_mean.data(), bh->fld.data());
        field3d_operators.subtract_mean_profile(b->fld.data(), b->fld_mean.data());
        field3d_operators.subtract_mean_profile(bh->fld.data(), bh->fld_mean.data());

        stats.set_mask_thres(mask_name, *b, *bh, 0., Stats_mask_type::Plus);

        fields.release_tmp(b);
        fields.release_tmp(bh);
    }
    else
    {
        std::string message = "Moist thermodynamics can not provide mask: \"" + mask_name +"\"";
        throw std::runtime_error(message);
    }
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
void Thermo_moist<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    tdep_pbot->update_time_dependent(bs.pbot, timeloop);
}

template<typename TF>
void Thermo_moist<TF>::get_thermo_field(Field3d<TF>& fld, std::string name, bool cyclic, bool is_stat)
{
    auto& gd = grid.get_grid_data();

    Background_state base;
    if (is_stat)
        base = bs_stats;
    else
        base = bs;

    // BvS: getThermoField() is called from subgrid-model, before thermo(), so re-calculate the hydrostatic pressure
    // Pass dummy as rhoref,bs.thvref to prevent overwriting base state
    if (bs.swupdatebasestate)
    {
        auto tmp = fields.get_tmp();
        calc_base_state(base.pref.data(), base.prefh.data(), &tmp->fld[0*gd.kcells], &tmp->fld[1*gd.kcells], &tmp->fld[2*gd.kcells],
                        &tmp->fld[3*gd.kcells], base.exnref.data(), base.exnrefh.data(), fields.sp.at("thl")->fld_mean.data(),
                        fields.sp.at("qt")->fld_mean.data(), base.pbot, gd.kstart, gd.kend, gd.z.data(), gd.dz.data(), gd.dzh.data());
        fields.release_tmp(tmp);
    }

    if (name == "b")
    {
        auto tmp = fields.get_tmp();
        calc_buoyancy(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.pref.data(), tmp->fld.data(), base.thvref.data(),
                      gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.kcells, gd.icells, gd.ijcells);
        fields.release_tmp(tmp);
    }
    else if (name == "b_h")
    {
        auto tmp = fields.get_tmp();
        calc_buoyancy_h(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.prefh.data(), base.thvrefh.data(),
                        &tmp->fld[0*gd.ijcells], &tmp->fld[1*gd.ijcells], &tmp->fld[2*gd.ijcells],
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        fields.release_tmp(tmp);
    }
    else if (name == "ql")
    {
        calc_liquid_water(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.pref.data(),
                          gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
    }
    else if (name == "ql_h")
    {
        auto tmp = fields.get_tmp();
        calc_liquid_water_h(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.prefh.data(), &tmp->fld[0*gd.ijcells], &tmp->fld[1*gd.ijcells],
                            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        fields.release_tmp(tmp);
    }
    else if (name == "N2")
    {
        calc_N2(fld.fld.data(), fields.sp.at("thl")->fld.data(), gd.dzi.data(), base.thvref.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
    }
    else if (name == "T")
    {
        calc_T(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.pref.data(), base.exnref.data(),
               gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
    }
    else if (name == "T_h")
    {
        auto tmp = fields.get_tmp();
        calc_T_h(fld.fld.data(), fields.sp.at("thl")->fld.data(), fields.sp.at("qt")->fld.data(), base.prefh.data(), &tmp->fld[0*gd.ijcells], &tmp->fld[1*gd.ijcells],
                 &tmp->fld[2*gd.ijcells], gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        fields.release_tmp(tmp);
    }
    else
    {
        std::string error_message = "Can not get thermo field: \"" + name + "\"";
        throw std::runtime_error(error_message);
    }

    if (cyclic)
        boundary_cyclic.exec(fld.fld.data());
}

template<typename TF>
void Thermo_moist<TF>::get_buoyancy_surf(Field3d<TF>& b, bool is_stat)
{
    auto& gd = grid.get_grid_data();
    Background_state base;
    if (is_stat)
        base = bs_stats;
    else
        base = bs;

    calc_buoyancy_bot(b.fld.data(), b.fld_bot.data(),
                      fields.sp.at("thl")->fld.data(), fields.sp.at("thl")->fld_bot.data(),
                      fields.sp.at("qt")->fld.data(), fields.sp.at("qt")->fld_bot.data(),
                      base.thvref.data(), base.thvrefh.data(), gd.icells, gd.jcells, gd.ijcells, gd.kstart);

    calc_buoyancy_fluxbot(b.flux_bot.data(), fields.sp.at("thl")->fld_bot.data(), fields.sp.at("thl")->flux_bot.data(),
                          fields.sp.at("qt")->fld_bot.data(), fields.sp.at("qt")->flux_bot.data(), base.thvrefh.data(),
                          gd.icells, gd.jcells, gd.kstart);
}

template<typename TF>
void Thermo_moist<TF>::get_buoyancy_fluxbot(Field3d<TF>& b, bool is_stat)
{
    auto& gd = grid.get_grid_data();
    Background_state base;
    if (is_stat)
        base = bs_stats;
    else
        base = bs;

    calc_buoyancy_fluxbot(b.flux_bot.data(), fields.sp.at("thl")->fld_bot.data(), fields.sp.at("thl")->flux_bot.data(),
                          fields.sp.at("qt")->fld_bot.data(), fields.sp.at("qt")->flux_bot.data(), base.thvrefh.data(),
                          gd.icells, gd.jcells, gd.kstart);
}

template<typename TF>
void Thermo_moist<TF>::get_T_bot(Field3d<TF>& T_bot, bool is_stat)
{
    auto& gd = grid.get_grid_data();
    Background_state base;
    if (is_stat)
        base = bs_stats;
    else
        base = bs;

    calc_T_bot(T_bot.fld_bot.data(), fields.sp.at("thl")->fld.data(), base.exnrefh.data(), base.thl0.data(),
               gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.icells, gd.ijcells);
}

template<typename TF>
const std::vector<TF>& Thermo_moist<TF>::get_p_vector() const
{
    return bs.pref;
}

template<typename TF>
const std::vector<TF>& Thermo_moist<TF>::get_ph_vector() const
{
    return bs.prefh;
}

template<typename TF>
const std::vector<TF>& Thermo_moist<TF>::get_exner_vector() const
{
    return bs.exnref;
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
    bs_stats = bs;

    // Add variables to the statistics
    if (stats.get_switch())
    {
        /* Add fixed base-state density and temperature profiles. Density should probably be in fields (?), but
           there the statistics are initialized before thermo->create() is called */
        stats.add_fixed_prof("rhoref",  "Full level basic state density", "kg m-3", "z",  bs.rhoref.data() );
        stats.add_fixed_prof("rhorefh", "Half level basic state density", "kg m-3", "zh", bs.rhorefh.data());
        stats.add_fixed_prof("thvref",  "Full level basic state virtual potential temperature", "K", "z", bs.thvref.data() );
        stats.add_fixed_prof("thvrefh", "Half level basic state virtual potential temperature", "K", "zh", bs.thvrefh.data());

        if (bs_stats.swupdatebasestate)
        {
            stats.add_prof("phydro", "Full level hydrostatic pressure", "Pa", "z" );
            stats.add_prof("phydroh","Half level hydrostatic pressure", "Pa", "zh");
            stats.add_prof("rho",  "Full level density", "kg m-3", "z" );
            stats.add_prof("rhoh", "Half level density", "kg m-3", "zh");
        }
        else
        {
            stats.add_fixed_prof("pydroh",  "Full level hydrostatic pressure", "Pa", "z",  bs.pref.data() );
            stats.add_fixed_prof("phydroh", "Half level hydrostatic pressure", "Pa", "zh", bs.prefh.data());
        }

        stats.add_prof("b", "Buoyancy", "m s-2", "z", Stats_whitelist_type::White);
        for (int n=2; n<5; ++n)
        {
            std::stringstream ss;
            ss << n;
            std::string sn = ss.str();
            stats.add_prof("b"+sn, "Moment " +sn+" of the buoyancy", "(m s-2)"+sn, "z");
        }

        stats.add_prof("bgrad", "Gradient of the buoyancy", "m s-3", "zh");
        stats.add_prof("bw"   , "Turbulent flux of the buoyancy", "m2 s-3", "zh");
        stats.add_prof("bdiff", "Diffusive flux of the buoyancy", "m2 s-3", "zh");
        stats.add_prof("bflux", "Total flux of the buoyancy", "m2 s-3", "zh");

        stats.add_prof("ql", "Liquid water mixing ratio", "kg kg-1", "z", Stats_whitelist_type::White);
        stats.add_prof("qlfrac", "Cloud fraction", "-", "z");

        stats.add_time_series("qlpath", "Liquid water path", "kg m-2");
        stats.add_time_series("qlcover", "Projected cloud cover", "-");
    }
}

template<typename TF>
void Thermo_moist<TF>::create_column(Column<TF>& column)
{
    // add the profiles to the columns
    if (column.get_switch())
    {
        column.add_prof("b", "Buoyancy", "m s-2", "z");
        column.add_prof("ql", "Liquid water mixing ratio", "kg kg-1", "z");
    }
}

template<typename TF>
void Thermo_moist<TF>::create_cross(Cross<TF>& cross)
{
    if (cross.get_switch())
    {
        swcross_b = false;
        swcross_ql = false;

        // Vectors with allowed cross variables for buoyancy and liquid water
        std::vector<std::string> allowed_crossvars_b = {"b", "bbot", "bfluxbot"};
        std::vector<std::string> allowed_crossvars_ql = {"ql", "qlpath", "qlbase", "qltop"};

        std::vector<std::string> bvars  = cross.get_enabled_variables(allowed_crossvars_b);
        std::vector<std::string> qlvars = cross.get_enabled_variables(allowed_crossvars_ql);

        if (bvars.size() > 0)
            swcross_b  = true;
        if (qlvars.size() > 0)
            swcross_ql = true;

        // Merge into one vector
        crosslist = bvars;
        crosslist.insert(crosslist.end(), qlvars.begin(), qlvars.end());
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
void Thermo_moist<TF>::exec_stats(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    #ifndef USECUDA
    bs_stats = bs;
    #endif

    const TF no_offset = 0.;
    const TF no_threshold = 0.;

    // calculate the buoyancy and its surface flux for the profiles
    auto b = fields.get_tmp();
    b->loc = gd.sloc;
    get_thermo_field(*b, "b", true, true);
    get_buoyancy_surf(*b, true);
    get_buoyancy_fluxbot(*b, true);

    // calculate the mean
    std::vector<std::string> operators = {"mean","2","3","4","w","grad","diff","flux"};

    stats.calc_stats("b", *b, no_offset, no_threshold, operators);

    fields.release_tmp(b);

    // calculate the liquid water stats
    auto ql = fields.get_tmp();
    ql->loc = gd.sloc;

    get_thermo_field(*ql, "ql", true, true);
    stats.calc_stats("ql", *ql, no_offset, no_threshold, {"mean","cover","frac","path"});


    fields.release_tmp(ql);

    if (bs_stats.swupdatebasestate)
    {
        stats.set_prof("phydro" , bs_stats.pref);
        stats.set_prof("phydroh", bs_stats.prefh);
        stats.set_prof("rho"    , bs_stats.rhoref);
        stats.set_prof("rhoh"   , bs_stats.rhorefh);
    }

}


#ifndef USECUDA
template<typename TF>
void Thermo_moist<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;
    auto output = fields.get_tmp();

    get_thermo_field(*output, "b",false, true);
    column.calc_column("b", output->fld.data(), no_offset);

    get_thermo_field(*output, "ql",false, true);
    column.calc_column("ql", output->fld.data(), no_offset);

    fields.release_tmp(output);
}
#endif


template<typename TF>
void Thermo_moist<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    #ifndef USECUDA
    bs_stats = bs;
    #endif

    auto output = fields.get_tmp();

    if(swcross_b)
    {
        get_thermo_field(*output, "b", false, true);
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
        get_thermo_field(*output, "ql", false, true);
    }
    for (auto& it : crosslist)
    {
        if (it == "ql")
            cross.cross_simple(output->fld.data(), "ql", iotime);
        if (it == "qlpath")
            cross.cross_path(output->fld.data(), "qlpath", iotime);
        if (it == "qlbase")
            cross.cross_height_threshold(output->fld.data(), 0., Cross_direction::Bottom_to_top, "qlbase", iotime);
        if (it == "qltop")
            cross.cross_height_threshold(output->fld.data(), 0., Cross_direction::Top_to_bottom, "qltop", iotime);
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
            get_thermo_field(*output, "b", false, true);
        else if (it == "ql")
            get_thermo_field(*output, "ql", false, true);
        else if (it == "T")
            get_thermo_field(*output, "T", false, true);
        else
        {
            master.print_error("Thermo dump of field \"%s\" not supported\n", it.c_str());
            throw std::runtime_error("Error in Thermo Dump");
        }
        dump.save_dump(output->fld.data(), it, iotime);
    }
    fields.release_tmp(output);
}

template class Thermo_moist<double>;
template class Thermo_moist<float>;
