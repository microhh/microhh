/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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
#include "budget.h"
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "fd.h"
#include "model.h"
#include "thermo.h"
#include "stats.h"
#include <netcdfcpp.h>

using namespace fd::o4; // budget is only 4th order

// #define NO_OFFSET 0.

Budget::Budget(Model *modelin, Input *inputin)
{
  model = modelin;

  umodel = 0;
  vmodel = 0;

  int nerror = 0;
  nerror += inputin->getItem(&swbudget, "budget", "swbudget", "", "0");

  if(nerror)
    throw 1;
}

Budget::~Budget()
{
  delete[] umodel;
  delete[] vmodel;
}

void Budget::init()
{
  // copy the pointers
  grid   = model->grid;
  fields = model->fields;
  stats  = model->stats;
  master = model->master;

  // if the stats is disabled, also disable the budget stats
  if(stats->getSwitch() == "0")
    swbudget = "0";

  umodel = new double[grid->kcells];
  vmodel = new double[grid->kcells];

  for(int k=0; k<grid->kcells; ++k)
  {
    umodel[k] = 0.;
    vmodel[k] = 0.;
  }
}

void Budget::create()
{
  if(swbudget == "0")
    return;

  // add the profiles for the kinetic energy to the statistics
  stats->addprof("ke" , "Kinetic energy" , "m2 s-2", "z");
  stats->addprof("tke", "Turbulent kinetic energy" , "m2 s-2", "z");

  // add the profiles for the kinetic energy budget to the statistics
  stats->addprof("u2_shear" , "Shear production term in U2 budget" , "m2 s-3", "z");
  stats->addprof("v2_shear" , "Shear production term in V2 budget" , "m2 s-3", "z");
  stats->addprof("tke_shear", "Shear production term in TKE budget", "m2 s-3", "z");

  stats->addprof("u2_turb" , "Turbulent transport term in U2 budget" , "m2 s-3", "z" );
  stats->addprof("v2_turb" , "Turbulent transport term in V2 budget" , "m2 s-3", "z" );
  stats->addprof("w2_turb" , "Turbulent transport term in W2 budget" , "m2 s-3", "zh");
  stats->addprof("tke_turb", "Turbulent transport term in TKE budget", "m2 s-3", "z" );

  stats->addprof("u2_visc" , "Viscous transport term in U2 budget" , "m2 s-3", "z" );
  stats->addprof("v2_visc" , "Viscous transport term in V2 budget" , "m2 s-3", "z" );
  stats->addprof("w2_visc" , "Viscous transport term in W2 budget" , "m2 s-3", "zh");
  stats->addprof("tke_visc", "Viscous transport term in TKE budget", "m2 s-3", "z" );

  stats->addprof("u2_diss" , "Dissipation term in U2 budget" , "m2 s-3", "z" );
  stats->addprof("v2_diss" , "Dissipation term in V2 budget" , "m2 s-3", "z" );
  stats->addprof("w2_diss" , "Dissipation term in W2 budget" , "m2 s-3", "zh");
  stats->addprof("tke_diss", "Dissipation term in TKE budget", "m2 s-3", "z" );

  stats->addprof("w2_pres" , "Pressure transport term in W2 budget" , "m2 s-3", "zh");
  stats->addprof("tke_pres", "Pressure transport term in TKE budget", "m2 s-3", "z" );

  stats->addprof("u2_rdstr", "Pressure redistribution term in U2 budget", "m2 s-3", "z" );
  stats->addprof("v2_rdstr", "Pressure redistribution term in V2 budget", "m2 s-3", "z" );
  stats->addprof("w2_rdstr", "Pressure redistribution term in W2 budget", "m2 s-3", "zh");

  if(model->thermo->getSwitch() != "0")
  {
    stats->addprof("w2_buoy" , "Buoyancy production/destruction term in W2 budget" , "m2 s-3", "zh");
    stats->addprof("tke_buoy", "Buoyancy production/destruction term in TKE budget", "m2 s-3", "z" );
  }

  if(model->thermo->getSwitch() != "0")
  {
    // add the profiles for the potential energy budget to the statistics
    stats->addprof("bsort", "Sorted buoyancy", "m s-2", "z");
    stats->addprof("zsort", "Height diff buoyancy and sorted buoyancy", "m", "z");
    stats->addprof("pe"   , "Total potential energy", "m2 s-2", "z");
    stats->addprof("ape"  , "Available potential energy", "m2 s-2", "z");
    stats->addprof("bpe"  , "Background potential energy", "m2 s-2", "z");

    // add the budget terms for the potential energy
    stats->addprof("pe_turb", "Turbulent transport term in potential energy budget", "m2 s-3", "z");
    stats->addprof("pe_visc", "Viscous transport term in potential energy budget", "m2 s-3", "z");
    stats->addprof("pe_bous", "Boussinesq term in potential energy budget", "m2 s-3", "z");

    // add the budget terms for the background potential energy
    stats->addprof("bpe_turb", "Turbulent transport term in background potential energy budget", "m2 s-3", "z");
    stats->addprof("bpe_visc", "Viscous transport term in background potential energy budget", "m2 s-3", "z");
    stats->addprof("bpe_diss", "Dissipation term in background potential energy budget", "m2 s-3", "z");
  }
}

int Budget::execStats(mask *m)
{
  if(swbudget == "0")
    return 0;

  // calculate the mean of the fields
  grid->calcMean(umodel, fields->u->data, grid->kcells);
  grid->calcMean(vmodel, fields->v->data, grid->kcells);

  if(grid->swspatialorder == "4")
  {
    // calculate the TKE budget
    calcke(fields->u->data, fields->v->data, fields->w->data,
           umodel, vmodel,
           grid->utrans, grid->vtrans,
           m->profs["ke"].data, m->profs["tke"].data);

    calctkebudget(fields->u->data, fields->v->data, fields->w->data, fields->sd["p"]->data,
                  fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data,
                  umodel, vmodel,
                  m->profs["u2_shear"].data, m->profs["v2_shear"].data, m->profs["tke_shear"].data,
                  m->profs["u2_turb"].data, m->profs["v2_turb"].data, m->profs["w2_turb"].data, m->profs["tke_turb"].data,
                  m->profs["u2_visc"].data, m->profs["v2_visc"].data, m->profs["w2_visc"].data, m->profs["tke_visc"].data,
                  m->profs["u2_diss"].data, m->profs["v2_diss"].data, m->profs["w2_diss"].data, m->profs["tke_diss"].data,
                  m->profs["w2_pres"].data, m->profs["tke_pres"].data,
                  m->profs["u2_rdstr"].data, m->profs["v2_rdstr"].data, m->profs["w2_rdstr"].data,
                  grid->dzi4, grid->dzhi4, fields->visc);

    // calculate the buoyancy term of the TKE budget
    if(model->thermo->getSwitch() != "0")
    {
      // store the buoyancy in the tmp1 field
      model->thermo->getThermoField(fields->atmp["tmp1"], fields->atmp["tmp2"], "b");
      calctkebudget_buoy(fields->w->data, fields->atmp["tmp1"]->data,
                    m->profs["w2_buoy"].data, m->profs["tke_buoy"].data);
    }

    // calculate the potential energy budget
    if(model->thermo->getSwitch() != "0")
    {
      // calculate the sorted buoyancy profile, tmp1 still contains the buoyancy
      stats->calcsortprof(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, m->profs["bsort"].data);

      // calculate the potential energy back, tmp1 contains the buoyancy, tmp2 will contain height that the local buoyancy
      // will reach in the sorted profile
      calcpe(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp2"]->databot, fields->atmp["tmp2"]->datatop,
             grid->z,
             m->profs["bsort"].data,
             m->profs["pe"].data, m->profs["ape"].data, m->profs["bpe"].data,
             m->profs["zsort"].data);


      // calculate the budget of background potential energy, start with this one, because tmp2 contains the needed height
      // which will be overwritten inside of the routine
      calcbpebudget(fields->w->data, fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp2"]->databot, fields->atmp["tmp2"]->datatop,
                    m->profs["bpe_turb"].data, m->profs["bpe_visc"].data, m->profs["bpe_diss"].data,
                    // TODO put the correct value for visc here!!!!!
                    m->profs["bsort"].data,
                    grid->z, grid->dzi4, grid->dzhi4,
                    fields->visc);

      // calculate the budget of potential energy
      calcpebudget(fields->w->data, fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp2"]->datatop,
                   m->profs["pe_turb"].data, m->profs["pe_visc"].data, m->profs["pe_bous"].data,
                   // TODO put the correct value for visc here!!!!!
                   grid->z, grid->zh, grid->dzi4, grid->dzhi4,
                   fields->visc);
    }
  }

  return 0;
}

int Budget::calcke(double * restrict u, double * restrict v, double * restrict w, 
                    double * restrict umodel, double * restrict vmodel,
                    double utrans, double vtrans,
                    double * restrict ke, double * restrict tke)
{
  int ijk,ii1,ii2,jj1,jj2,kk1,kk2;
  double u2,v2,w2;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    ke [k] = 0;
    tke[k] = 0;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj1 + k*kk1;
        u2 = ci0*std::pow(u[ijk-ii1] + utrans, 2) + ci1*std::pow(u[ijk    ] + utrans, 2) 
           + ci2*std::pow(u[ijk+ii1] + utrans, 2) + ci3*std::pow(u[ijk+ii2] + utrans, 2);
        v2 = ci0*std::pow(v[ijk-jj1] + vtrans, 2) + ci1*std::pow(v[ijk    ] + vtrans, 2)
           + ci2*std::pow(v[ijk+jj1] + vtrans, 2) + ci3*std::pow(v[ijk+jj2] + vtrans, 2);
        w2 = ci0*std::pow(w[ijk-kk1], 2) + ci1*std::pow(w[ijk], 2) + ci2*std::pow(w[ijk+kk1], 2) + ci3*std::pow(w[ijk+kk2], 2);
        ke[k] += 0.5*(u2 + v2 + w2);
      }

    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj1 + k*kk1;
        u2 = ci0*std::pow(u[ijk-ii1] - umodel[k], 2) + ci1*std::pow(u[ijk    ] - umodel[k], 2) 
           + ci2*std::pow(u[ijk+ii1] - umodel[k], 2) + ci3*std::pow(u[ijk+ii2] - umodel[k], 2);
        v2 = ci0*std::pow(v[ijk-jj1] - vmodel[k], 2) + ci1*std::pow(v[ijk    ] - vmodel[k], 2)
           + ci2*std::pow(v[ijk+jj1] - vmodel[k], 2) + ci3*std::pow(v[ijk+jj2] - vmodel[k], 2);
        w2 = ci0*std::pow(w[ijk-kk1], 2) + ci1*std::pow(w[ijk], 2) + ci2*std::pow(w[ijk+kk1], 2) + ci3*std::pow(w[ijk+kk2], 2);
        tke[k] += 0.5*(u2 + v2 + w2);
      }
  }

  master->sum(ke , grid->kcells);
  master->sum(tke, grid->kcells);

  int n = grid->itot*grid->jtot;
  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    ke [k] /= n;
    tke[k] /= n;
  }

  return 0;
}

int Budget::calctkebudget(double * restrict u, double * restrict v, double * restrict w, double * restrict p,
                           double * restrict wx, double * restrict wy,
                           double * restrict umean, double * restrict vmean,
                           double * restrict u2_shear, double * restrict v2_shear, double * restrict tke_shear,
                           double * restrict u2_turb, double * restrict v2_turb, double * restrict w2_turb, double * restrict tke_turb,
                           double * restrict u2_visc, double * restrict v2_visc, double * restrict w2_visc, double * restrict tke_visc,
                           double * restrict u2_diss, double * restrict v2_diss, double * restrict w2_diss, double * restrict tke_diss,
                           double * restrict w2_pres, double * restrict tke_pres,
                           double * restrict u2_rdstr, double * restrict v2_rdstr, double * restrict w2_rdstr,
                           double * restrict dzi4, double * restrict dzhi4, double visc)
{
  // 1. INTERPOLATE THE VERTICAL VELOCITY TO U AND V LOCATION
  const int wloc [3] = {0,0,1};
  const int wxloc[3] = {1,0,1};
  const int wyloc[3] = {0,1,1};

  grid->interpolate_4th(wx, w, wloc, wxloc);
  grid->interpolate_4th(wy, w, wloc, wyloc);

  int ijk,ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  int kstart,kend;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;
  kk3 = 3*grid->ijcells;

  kstart = grid->kstart;
  kend   = grid->kend;

  double n = grid->itot*grid->jtot;

  // 2. CALCULATE THE SHEAR TERM u'w*dumean/dz
  // bottom boundary
  int k = grid->kstart;
  u2_shear [k] = 0.;
  v2_shear [k] = 0.;
  tke_shear[k] = 0.;

  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      u2_shear[k] -= 2.*(u[ijk]-umean[k])*(ci0*wx[ijk-kk1] + ci1*wx[ijk] + ci2*wx[ijk+kk1] + ci3*wx[ijk+kk2])
                   * ( cg0*(bi0*umean[k-2] + bi1*umean[k-1] + bi2*umean[k  ] + bi3*umean[k+1])
                     + cg1*(ci0*umean[k-2] + ci1*umean[k-1] + ci2*umean[k  ] + ci3*umean[k+1])
                     + cg2*(ci0*umean[k-1] + ci1*umean[k  ] + ci2*umean[k+1] + ci3*umean[k+2])
                     + cg3*(ci0*umean[k  ] + ci1*umean[k+1] + ci2*umean[k+2] + ci3*umean[k+3])) * dzi4[k];

      v2_shear[k] -= 2.*(v[ijk]-vmean[k])*(ci0*wy[ijk-kk1] + ci1*wy[ijk] + ci2*wy[ijk+kk1] + ci3*wy[ijk+kk2])
                   * ( cg0*(bi0*vmean[k-2] + bi1*vmean[k-1] + bi2*vmean[k  ] + bi3*vmean[k+1])
                     + cg1*(ci0*vmean[k-2] + ci1*vmean[k-1] + ci2*vmean[k  ] + ci3*vmean[k+1])
                     + cg2*(ci0*vmean[k-1] + ci1*vmean[k  ] + ci2*vmean[k+1] + ci3*vmean[k+2])
                     + cg3*(ci0*vmean[k  ] + ci1*vmean[k+1] + ci2*vmean[k+2] + ci3*vmean[k+3])) * dzi4[k];
    }
  tke_shear[k] += 0.5*(u2_shear[k] + v2_shear[k]);

  // interior
  for(int k=grid->kstart+1; k<grid->kend-1; ++k)
  {
    u2_shear [k] = 0.;
    v2_shear [k] = 0.;
    tke_shear[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        u2_shear[k] -= 2.*(u[ijk]-umean[k])*(ci0*wx[ijk-kk1] + ci1*wx[ijk] + ci2*wx[ijk+kk1] + ci3*wx[ijk+kk2])
                     * ( cg0*(ci0*umean[k-3] + ci1*umean[k-2] + ci2*umean[k-1] + ci3*umean[k  ])
                       + cg1*(ci0*umean[k-2] + ci1*umean[k-1] + ci2*umean[k  ] + ci3*umean[k+1])
                       + cg2*(ci0*umean[k-1] + ci1*umean[k  ] + ci2*umean[k+1] + ci3*umean[k+2])
                       + cg3*(ci0*umean[k  ] + ci1*umean[k+1] + ci2*umean[k+2] + ci3*umean[k+3])) * dzi4[k];

        v2_shear[k] -= 2.*(v[ijk]-vmean[k])*(ci0*wy[ijk-kk1] + ci1*wy[ijk] + ci2*wy[ijk+kk1] + ci3*wy[ijk+kk2])
                     * ( cg0*(ci0*vmean[k-3] + ci1*vmean[k-2] + ci2*vmean[k-1] + ci3*vmean[k  ])
                       + cg1*(ci0*vmean[k-2] + ci1*vmean[k-1] + ci2*vmean[k  ] + ci3*vmean[k+1])
                       + cg2*(ci0*vmean[k-1] + ci1*vmean[k  ] + ci2*vmean[k+1] + ci3*vmean[k+2])
                       + cg3*(ci0*vmean[k  ] + ci1*vmean[k+1] + ci2*vmean[k+2] + ci3*vmean[k+3])) * dzi4[k];
      }
    tke_shear[k] += 0.5*(u2_shear[k] + v2_shear[k]);
  }

  // top boundary
  k = grid->kend-1;

  u2_shear [k] = 0.;
  v2_shear [k] = 0.;
  tke_shear[k] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      u2_shear[k] -= 2.*(u[ijk]-umean[k])*(ci0*wx[ijk-kk1] + ci1*wx[ijk] + ci2*wx[ijk+kk1] + ci3*wx[ijk+kk2])
                   * ( cg0*(ci0*umean[k-3] + ci1*umean[k-2] + ci2*umean[k-1] + ci3*umean[k  ])
                     + cg1*(ci0*umean[k-2] + ci1*umean[k-1] + ci2*umean[k  ] + ci3*umean[k+1])
                     + cg2*(ci0*umean[k-1] + ci1*umean[k  ] + ci2*umean[k+1] + ci3*umean[k+2])
                     + cg3*(ti0*umean[k  ] + ti1*umean[k+1] + ti2*umean[k+2] + ti3*umean[k+3])) * dzi4[k];

      v2_shear[k] -= 2.*(v[ijk]-vmean[k])*(ci0*wy[ijk-kk1] + ci1*wy[ijk] + ci2*wy[ijk+kk1] + ci3*wy[ijk+kk2])
                   * ( cg0*(ci0*vmean[k-3] + ci1*vmean[k-2] + ci2*vmean[k-1] + ci3*vmean[k  ])
                     + cg1*(ci0*vmean[k-2] + ci1*vmean[k-1] + ci2*vmean[k  ] + ci3*vmean[k+1])
                     + cg2*(ci0*vmean[k-1] + ci1*vmean[k  ] + ci2*vmean[k+1] + ci3*vmean[k+2])
                     + cg3*(ti0*vmean[k-1] + ti1*vmean[k  ] + ti2*vmean[k+1] + ti3*vmean[k+2])) * dzi4[k];
    }
  tke_shear[k] += 0.5*(u2_shear[k] + v2_shear[k]);

  // create the profiles
  master->sum(u2_shear, grid->kcells);
  master->sum(v2_shear, grid->kcells);
  master->sum(tke_shear, grid->kcells);

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    u2_shear [k] /= n;
    v2_shear [k] /= n;
    tke_shear[k] /= n;
  }

  // 3. CALCULATE TURBULENT FLUXES
  // bottom boundary
  k = grid->kstart;

  u2_turb [k] = 0.;
  v2_turb [k] = 0.;
  tke_turb[k] = 0.;

  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      u2_turb[k]  -= ( cg0*((bi0*std::pow(u[ijk-kk2]-umean[k-2],2) + bi1*std::pow(u[ijk-kk1]-umean[k-1],2) + bi2*std::pow(u[ijk    ]-umean[k  ],2) + bi3*std::pow(u[ijk+kk1]-umean[k+1],2))*wx[ijk-kk1])
                     + cg1*((ci0*std::pow(u[ijk-kk2]-umean[k-2],2) + ci1*std::pow(u[ijk-kk1]-umean[k-1],2) + ci2*std::pow(u[ijk    ]-umean[k  ],2) + ci3*std::pow(u[ijk+kk1]-umean[k+1],2))*wx[ijk    ])
                     + cg2*((ci0*std::pow(u[ijk-kk1]-umean[k-1],2) + ci1*std::pow(u[ijk    ]-umean[k  ],2) + ci2*std::pow(u[ijk+kk1]-umean[k+1],2) + ci3*std::pow(u[ijk+kk2]-umean[k+2],2))*wx[ijk+kk1])
                     + cg3*((ci0*std::pow(u[ijk    ]-umean[k  ],2) + ci1*std::pow(u[ijk+kk1]-umean[k+1],2) + ci2*std::pow(u[ijk+kk2]-umean[k+2],2) + ci3*std::pow(u[ijk+kk3]-umean[k+3],2))*wx[ijk+kk2]) ) * dzi4[k];

      v2_turb[k]  -= ( cg0*((bi0*std::pow(v[ijk-kk2]-vmean[k-2],2) + bi1*std::pow(v[ijk-kk1]-vmean[k-1],2) + bi2*std::pow(v[ijk    ]-vmean[k  ],2) + bi3*std::pow(v[ijk+kk1]-vmean[k+1],2))*wy[ijk-kk1])
                     + cg1*((ci0*std::pow(v[ijk-kk2]-vmean[k-2],2) + ci1*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci2*std::pow(v[ijk    ]-vmean[k  ],2) + ci3*std::pow(v[ijk+kk1]-vmean[k+1],2))*wy[ijk    ]) 
                     + cg2*((ci0*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci1*std::pow(v[ijk    ]-vmean[k  ],2) + ci2*std::pow(v[ijk+kk1]-vmean[k+1],2) + ci3*std::pow(v[ijk+kk2]-vmean[k+2],2))*wy[ijk+kk1]) 
                     + cg3*((ci0*std::pow(v[ijk    ]-vmean[k  ],2) + ci1*std::pow(v[ijk+kk1]-vmean[k+1],2) + ci2*std::pow(v[ijk+kk2]-vmean[k+2],2) + ci3*std::pow(v[ijk+kk3]-vmean[k+3],2))*wy[ijk+kk2]) ) * dzi4[k];

      tke_turb[k] -= 0.5*( cg0*std::pow(w[ijk-kk1], 3) + cg1*std::pow(w[ijk], 3) + cg2*std::pow(w[ijk+kk1], 3) + cg3*std::pow(w[ijk+kk2], 3)) * dzi4[k];
    }
  tke_turb[k] += 0.5*(u2_turb[k] + v2_turb[k]);

  // interior
  for(int k=grid->kstart+1; k<grid->kend-1; ++k)
  {
    u2_turb [k] = 0.;
    v2_turb [k] = 0.;
    tke_turb[k] = 0.;

    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        u2_turb[k]  -= ( cg0*((ci0*std::pow(u[ijk-kk3]-umean[k-3],2) + ci1*std::pow(u[ijk-kk2]-umean[k-2],2) + ci2*std::pow(u[ijk-kk1]-umean[k-1],2) + ci3*std::pow(u[ijk    ]-umean[k  ],2))*wx[ijk-kk1])
                       + cg1*((ci0*std::pow(u[ijk-kk2]-umean[k-2],2) + ci1*std::pow(u[ijk-kk1]-umean[k-1],2) + ci2*std::pow(u[ijk    ]-umean[k  ],2) + ci3*std::pow(u[ijk+kk1]-umean[k+1],2))*wx[ijk    ])
                       + cg2*((ci0*std::pow(u[ijk-kk1]-umean[k-1],2) + ci1*std::pow(u[ijk    ]-umean[k  ],2) + ci2*std::pow(u[ijk+kk1]-umean[k+1],2) + ci3*std::pow(u[ijk+kk2]-umean[k+2],2))*wx[ijk+kk1])
                       + cg3*((ci0*std::pow(u[ijk    ]-umean[k  ],2) + ci1*std::pow(u[ijk+kk1]-umean[k+1],2) + ci2*std::pow(u[ijk+kk2]-umean[k+2],2) + ci3*std::pow(u[ijk+kk3]-umean[k+3],2))*wx[ijk+kk2]) ) * dzi4[k];

        v2_turb[k]  -= ( cg0*((ci0*std::pow(v[ijk-kk3]-vmean[k-3],2) + ci1*std::pow(v[ijk-kk2]-vmean[k-2],2) + ci2*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci3*std::pow(v[ijk    ]-vmean[k  ],2))*wy[ijk-kk1]) 
                       + cg1*((ci0*std::pow(v[ijk-kk2]-vmean[k-2],2) + ci1*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci2*std::pow(v[ijk    ]-vmean[k  ],2) + ci3*std::pow(v[ijk+kk1]-vmean[k+1],2))*wy[ijk    ]) 
                       + cg2*((ci0*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci1*std::pow(v[ijk    ]-vmean[k  ],2) + ci2*std::pow(v[ijk+kk1]-vmean[k+1],2) + ci3*std::pow(v[ijk+kk2]-vmean[k+2],2))*wy[ijk+kk1]) 
                       + cg3*((ci0*std::pow(v[ijk    ]-vmean[k  ],2) + ci1*std::pow(v[ijk+kk1]-vmean[k+1],2) + ci2*std::pow(v[ijk+kk2]-vmean[k+2],2) + ci3*std::pow(v[ijk+kk3]-vmean[k+3],2))*wy[ijk+kk2]) ) * dzi4[k];

        tke_turb[k] -= 0.5*( cg0*std::pow(w[ijk-kk1], 3) + cg1*std::pow(w[ijk], 3) + cg2*std::pow(w[ijk+kk1], 3) + cg3*std::pow(w[ijk+kk2], 3)) * dzi4[k];
      }
    tke_turb[k] += 0.5*(u2_turb[k] + v2_turb[k]);
  }

  // top boundary
  k = grid->kend-1;

  u2_turb [k] = 0.;
  v2_turb [k] = 0.;
  tke_turb[k] = 0.;

  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      u2_turb[k]  -= ( cg0*((ci0*std::pow(u[ijk-kk3]-umean[k-3],2) + ci1*std::pow(u[ijk-kk2]-umean[k-2],2) + ci2*std::pow(u[ijk-kk1]-umean[k-1],2) + ci3*std::pow(u[ijk    ]-umean[k  ],2))*wx[ijk-kk1])
                     + cg1*((ci0*std::pow(u[ijk-kk2]-umean[k-2],2) + ci1*std::pow(u[ijk-kk1]-umean[k-1],2) + ci2*std::pow(u[ijk    ]-umean[k  ],2) + ci3*std::pow(u[ijk+kk1]-umean[k+1],2))*wx[ijk    ])
                     + cg2*((ci0*std::pow(u[ijk-kk1]-umean[k-1],2) + ci1*std::pow(u[ijk    ]-umean[k  ],2) + ci2*std::pow(u[ijk+kk1]-umean[k+1],2) + ci3*std::pow(u[ijk+kk2]-umean[k+2],2))*wx[ijk+kk1])
                     + cg3*((ti0*std::pow(u[ijk-kk1]-umean[k-1],2) + ti1*std::pow(u[ijk    ]-umean[k  ],2) + ti2*std::pow(u[ijk+kk1]-umean[k+1],2) + ti3*std::pow(u[ijk+kk2]-umean[k+2],2))*wx[ijk+kk1]) ) * dzi4[k];

      v2_turb[k]  -= ( cg0*((ci0*std::pow(v[ijk-kk3]-vmean[k-3],2) + ci1*std::pow(v[ijk-kk2]-vmean[k-2],2) + ci2*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci3*std::pow(v[ijk    ]-vmean[k  ],2))*wy[ijk-kk1]) 
                     + cg1*((ci0*std::pow(v[ijk-kk2]-vmean[k-2],2) + ci1*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci2*std::pow(v[ijk    ]-vmean[k  ],2) + ci3*std::pow(v[ijk+kk1]-vmean[k+1],2))*wy[ijk    ]) 
                     + cg2*((ci0*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci1*std::pow(v[ijk    ]-vmean[k  ],2) + ci2*std::pow(v[ijk+kk1]-vmean[k+1],2) + ci3*std::pow(v[ijk+kk2]-vmean[k+2],2))*wy[ijk+kk1]) 
                     + cg3*((ti0*std::pow(v[ijk-kk1]-vmean[k-1],2) + ti1*std::pow(v[ijk    ]-vmean[k  ],2) + ti2*std::pow(v[ijk+kk1]-vmean[k+1],2) + ti3*std::pow(v[ijk+kk2]-vmean[k+2],2))*wy[ijk+kk1]) ) * dzi4[k];

      tke_turb[k] -= 0.5*( cg0*std::pow(w[ijk-kk1], 3) + cg1*std::pow(w[ijk], 3) + cg2*std::pow(w[ijk+kk1], 3) + cg3*std::pow(w[ijk+kk2], 3)) * dzi4[k];
    }
  tke_turb[k] += 0.5*(u2_turb[k] + v2_turb[k]);

  // calculate the vertical velocity term
  k = grid->kstart;
  w2_turb[k] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      w2_turb[k] -= ( cg0*(bi0*std::pow(w[ijk-kk2],3) + bi1*std::pow(w[ijk-kk1],3) + bi2*std::pow(w[ijk    ],3) + bi3*std::pow(w[ijk+kk1],3))
                    + cg1*(ci0*std::pow(w[ijk-kk2],3) + ci1*std::pow(w[ijk-kk1],3) + ci2*std::pow(w[ijk    ],3) + ci3*std::pow(w[ijk+kk1],3))
                    + cg2*(ci0*std::pow(w[ijk-kk1],3) + ci1*std::pow(w[ijk    ],3) + ci2*std::pow(w[ijk+kk1],3) + ci3*std::pow(w[ijk+kk2],3))
                    + cg3*(ci0*std::pow(w[ijk    ],3) + ci1*std::pow(w[ijk+kk1],3) + ci2*std::pow(w[ijk+kk2],3) + ci3*std::pow(w[ijk+kk3],3)) ) * dzhi4[k];
    }

  for(int k=grid->kstart+1; k<grid->kend; ++k)
  {
    w2_turb[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        w2_turb[k] -= ( cg0*(ci0*std::pow(w[ijk-kk3],3) + ci1*std::pow(w[ijk-kk2],3) + ci2*std::pow(w[ijk-kk1],3) + ci3*std::pow(w[ijk    ],3))
                      + cg1*(ci0*std::pow(w[ijk-kk2],3) + ci1*std::pow(w[ijk-kk1],3) + ci2*std::pow(w[ijk    ],3) + ci3*std::pow(w[ijk+kk1],3))
                      + cg2*(ci0*std::pow(w[ijk-kk1],3) + ci1*std::pow(w[ijk    ],3) + ci2*std::pow(w[ijk+kk1],3) + ci3*std::pow(w[ijk+kk2],3))
                      + cg3*(ci0*std::pow(w[ijk    ],3) + ci1*std::pow(w[ijk+kk1],3) + ci2*std::pow(w[ijk+kk2],3) + ci3*std::pow(w[ijk+kk3],3)) ) * dzhi4[k];
      }
  }

  k = grid->kend;
  w2_turb[k] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      w2_turb[k] -= ( cg0*(ci0*std::pow(w[ijk-kk3],3) + ci1*std::pow(w[ijk-kk2],3) + ci2*std::pow(w[ijk-kk1],3) + ci3*std::pow(w[ijk    ],3))
                    + cg1*(ci0*std::pow(w[ijk-kk2],3) + ci1*std::pow(w[ijk-kk1],3) + ci2*std::pow(w[ijk    ],3) + ci3*std::pow(w[ijk+kk1],3))
                    + cg2*(ci0*std::pow(w[ijk-kk1],3) + ci1*std::pow(w[ijk    ],3) + ci2*std::pow(w[ijk+kk1],3) + ci3*std::pow(w[ijk+kk2],3))
                    + cg3*(ti0*std::pow(w[ijk-kk1],3) + ti1*std::pow(w[ijk    ],3) + ti2*std::pow(w[ijk+kk1],3) + ti3*std::pow(w[ijk+kk2],3)) ) * dzhi4[k];
    }

  // calculate the profiles
  master->sum(u2_turb , grid->kcells);
  master->sum(v2_turb , grid->kcells);
  master->sum(w2_turb , grid->kcells);
  master->sum(tke_turb, grid->kcells);

  for(k=grid->kstart; k<grid->kend; ++k)
  {
    u2_turb [k] /= n;
    v2_turb [k] /= n;
    tke_turb[k] /= n;
  }

  for(k=grid->kstart; k<grid->kend+1; ++k)
    w2_turb [k] /= n;

  // 4. CALCULATE THE PRESSURE TRANSPORT TERM
  // bottom boundary
  k = grid->kstart;
  tke_pres[k] = 0.;

  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      tke_pres[k] -= ( cg0*((bi0*p[ijk-kk2] + bi1*p[ijk-kk1] + bi2*p[ijk    ] + bi3*p[ijk+kk1])*w[ijk-kk1])
                     + cg1*((ci0*p[ijk-kk2] + ci1*p[ijk-kk1] + ci2*p[ijk    ] + ci3*p[ijk+kk1])*w[ijk    ])
                     + cg2*((ci0*p[ijk-kk1] + ci1*p[ijk    ] + ci2*p[ijk+kk1] + ci3*p[ijk+kk2])*w[ijk+kk1])
                     + cg3*((ci0*p[ijk    ] + ci1*p[ijk+kk1] + ci2*p[ijk+kk2] + ci3*p[ijk+kk3])*w[ijk+kk2]) ) * dzi4[k];
    }

  // interior
  for(int k=grid->kstart+1; k<grid->kend-1; ++k)
  {
    tke_pres[k] = 0.;

    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        tke_pres[k] -= ( cg0*((ci0*p[ijk-kk3] + ci1*p[ijk-kk2] + ci2*p[ijk-kk1] + ci3*p[ijk    ])*w[ijk-kk1])
                       + cg1*((ci0*p[ijk-kk2] + ci1*p[ijk-kk1] + ci2*p[ijk    ] + ci3*p[ijk+kk1])*w[ijk    ])
                       + cg2*((ci0*p[ijk-kk1] + ci1*p[ijk    ] + ci2*p[ijk+kk1] + ci3*p[ijk+kk2])*w[ijk+kk1])
                       + cg3*((ci0*p[ijk    ] + ci1*p[ijk+kk1] + ci2*p[ijk+kk2] + ci3*p[ijk+kk3])*w[ijk+kk2]) ) * dzi4[k];
      }
  }

  // top boundary
  k = grid->kend-1;
  tke_pres[k] = 0.;

  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      tke_pres[k] -= ( cg0*((ci0*p[ijk-kk3] + ci1*p[ijk-kk2] + ci2*p[ijk-kk1] + ci3*p[ijk    ])*w[ijk-kk1])
                     + cg1*((ci0*p[ijk-kk2] + ci1*p[ijk-kk1] + ci2*p[ijk    ] + ci3*p[ijk+kk1])*w[ijk    ])
                     + cg2*((ci0*p[ijk-kk1] + ci1*p[ijk    ] + ci2*p[ijk+kk1] + ci3*p[ijk+kk2])*w[ijk+kk1])
                     + cg3*((ti0*p[ijk-kk1] + ti1*p[ijk    ] + ti2*p[ijk+kk1] + ti3*p[ijk+kk2])*w[ijk+kk2]) ) * dzi4[k];
    }
 
  // calculate the vertical velocity pressure transport term
  // \TODO implement the proper BC as soon as the full BC's for pressure are added
  // bottom boundary
  k = grid->kstart;
  w2_pres[k] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      w2_pres[k] -= 0.*( cg0*((bi0*w[ijk-kk2] + bi1*w[ijk-kk1] + bi2*w[ijk    ] + bi3*w[ijk+kk1])*p[ijk-kk2])
                       + cg1*((ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1])*p[ijk-kk1])
                       + cg2*((ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])*p[ijk    ])
                       + cg3*((ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3])*p[ijk+kk1]) ) * dzhi4[k];
    }

  // interior
  for(int k=grid->kstart+1; k<grid->kend; ++k)
  {
    w2_pres[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        w2_pres[k] -= 2.*( cg0*((ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ])*p[ijk-kk2])
                         + cg1*((ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1])*p[ijk-kk1])
                         + cg2*((ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])*p[ijk    ])
                         + cg3*((ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3])*p[ijk+kk1]) ) * dzhi4[k];
      }
  }

  // top boundary
  k = grid->kend;
  w2_pres[k] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      w2_pres[k] -= 0.*( cg0*((ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ])*p[ijk-kk2])
                       + cg1*((ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1])*p[ijk-kk1])
                       + cg2*((ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])*p[ijk    ])
                       + cg3*((ti0*w[ijk-kk1] + ti1*w[ijk    ] + ti2*w[ijk+kk1] + ti3*w[ijk+kk2])*p[ijk+kk1]) ) * dzhi4[k];
    }

  master->sum(w2_pres , grid->kcells);
  master->sum(tke_pres, grid->kcells);

  for(int k=grid->kstart; k<grid->kend; ++k)
    tke_pres[k] /= n;

  for(int k=grid->kstart; k<grid->kend+1; ++k)
    w2_pres [k] /= n;
  
  // 5. CALCULATE THE VISCOUS TRANSPORT TERM
  // first, interpolate the vertical velocity to the scalar levels using temporary array wx
  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj1 + k*kk1;
        wx[ijk] = ci0*w[ijk-kk1] + ci1*w[ijk] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2];
      }

  // calculate the ghost cells at the bottom
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk = i + j*jj1 + kstart*kk1;
      wx[ijk-kk1] = - 2.*wx[ijk] + (1./3.)*wx[ijk+kk1];
      wx[ijk-kk2] = - 9.*wx[ijk] + 2.*wx[ijk+kk1];
    }

  // calculate the ghost cells at the top
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      wx[ijk+kk1] = - 2.*wx[ijk] + (1./3.)*wx[ijk-kk1];
      wx[ijk+kk2] = - 9.*wx[ijk] + 2.*wx[ijk-kk1];
    }

  // bottom boundary
  k = grid->kstart;

  u2_visc [k] = 0.;
  v2_visc [k] = 0.;
  tke_visc[k] = 0.;

  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      u2_visc[k]  += visc * ( cg0*((bg0*std::pow(u[ijk-kk2]-umean[k-2],2) + bg1*std::pow(u[ijk-kk1]-umean[k-1],2) + bg2*std::pow(u[ijk    ]-umean[k  ],2) + bg3*std::pow(u[ijk+kk1]-umean[k+1],2)) * dzhi4[k-1])
                            + cg1*((cg0*std::pow(u[ijk-kk2]-umean[k-2],2) + cg1*std::pow(u[ijk-kk1]-umean[k-1],2) + cg2*std::pow(u[ijk    ]-umean[k  ],2) + cg3*std::pow(u[ijk+kk1]-umean[k+1],2)) * dzhi4[k  ])
                            + cg2*((cg0*std::pow(u[ijk-kk1]-umean[k-1],2) + cg1*std::pow(u[ijk    ]-umean[k  ],2) + cg2*std::pow(u[ijk+kk1]-umean[k+1],2) + cg3*std::pow(u[ijk+kk2]-umean[k+2],2)) * dzhi4[k+1])
                            + cg3*((cg0*std::pow(u[ijk    ]-umean[k  ],2) + cg1*std::pow(u[ijk+kk1]-umean[k+1],2) + cg2*std::pow(u[ijk+kk2]-umean[k+2],2) + cg3*std::pow(u[ijk+kk3]-umean[k+3],2)) * dzhi4[k+2]) ) * dzi4[k];

      v2_visc[k]  += visc * ( cg0*((bg0*std::pow(v[ijk-kk2]-vmean[k-2],2) + bg1*std::pow(v[ijk-kk1]-vmean[k-1],2) + bg2*std::pow(v[ijk    ]-vmean[k  ],2) + bg3*std::pow(v[ijk+kk1]-vmean[k+1],2)) * dzhi4[k-1])
                            + cg1*((cg0*std::pow(v[ijk-kk2]-vmean[k-2],2) + cg1*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg2*std::pow(v[ijk    ]-vmean[k  ],2) + cg3*std::pow(v[ijk+kk1]-vmean[k+1],2)) * dzhi4[k  ])
                            + cg2*((cg0*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg1*std::pow(v[ijk    ]-vmean[k  ],2) + cg2*std::pow(v[ijk+kk1]-vmean[k+1],2) + cg3*std::pow(v[ijk+kk2]-vmean[k+2],2)) * dzhi4[k+1])
                            + cg3*((cg0*std::pow(v[ijk    ]-vmean[k  ],2) + cg1*std::pow(v[ijk+kk1]-vmean[k+1],2) + cg2*std::pow(v[ijk+kk2]-vmean[k+2],2) + cg3*std::pow(v[ijk+kk3]-vmean[k+3],2)) * dzhi4[k+2]) ) * dzi4[k];

      tke_visc[k] += 0.5 * visc * ( cg0*((bg0*std::pow(wx[ijk-kk2],2) + bg1*std::pow(wx[ijk-kk1],2) + bg2*std::pow(wx[ijk    ],2) + bg3*std::pow(wx[ijk+kk1],2)) * dzhi4[k-1])
                                  + cg1*((cg0*std::pow(wx[ijk-kk2],2) + cg1*std::pow(wx[ijk-kk1],2) + cg2*std::pow(wx[ijk    ],2) + cg3*std::pow(wx[ijk+kk1],2)) * dzhi4[k  ])
                                  + cg2*((cg0*std::pow(wx[ijk-kk1],2) + cg1*std::pow(wx[ijk    ],2) + cg2*std::pow(wx[ijk+kk1],2) + cg3*std::pow(wx[ijk+kk2],2)) * dzhi4[k+1])
                                  + cg3*((cg0*std::pow(wx[ijk    ],2) + cg1*std::pow(wx[ijk+kk1],2) + cg2*std::pow(wx[ijk+kk2],2) + cg3*std::pow(wx[ijk+kk3],2)) * dzhi4[k+2]) ) * dzi4[k];
      // tke_visc[k] += 0.5*visc * ( cg0*std::pow(w[ijk-kk1],2) + cg1*std::pow(w[ijk],2) + cg2*std::pow(w[ijk+kk1],2) + cg3*std::pow(w[ijk+kk2],2)) * dzi4[k];
    }
    tke_visc[k] += 0.5*(u2_visc[k] + v2_visc[k]);

  // interior
  for(int k=grid->kstart+1; k<grid->kend-1; ++k)
  {
    u2_visc [k] = 0.;
    v2_visc [k] = 0.;
    tke_visc[k] = 0.;

    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        u2_visc[k]  += visc * ( cg0*((cg0*std::pow(u[ijk-kk3]-umean[k-3],2) + cg1*std::pow(u[ijk-kk2]-umean[k-2],2) + cg2*std::pow(u[ijk-kk1]-umean[k-1],2) + cg3*std::pow(u[ijk    ]-umean[k  ],2)) * dzhi4[k-1])
                              + cg1*((cg0*std::pow(u[ijk-kk2]-umean[k-2],2) + cg1*std::pow(u[ijk-kk1]-umean[k-1],2) + cg2*std::pow(u[ijk    ]-umean[k  ],2) + cg3*std::pow(u[ijk+kk1]-umean[k+1],2)) * dzhi4[k  ])
                              + cg2*((cg0*std::pow(u[ijk-kk1]-umean[k-1],2) + cg1*std::pow(u[ijk    ]-umean[k  ],2) + cg2*std::pow(u[ijk+kk1]-umean[k+1],2) + cg3*std::pow(u[ijk+kk2]-umean[k+2],2)) * dzhi4[k+1])
                              + cg3*((cg0*std::pow(u[ijk    ]-umean[k  ],2) + cg1*std::pow(u[ijk+kk1]-umean[k+1],2) + cg2*std::pow(u[ijk+kk2]-umean[k+2],2) + cg3*std::pow(u[ijk+kk3]-umean[k+3],2)) * dzhi4[k+2]) ) * dzi4[k];

        v2_visc[k]  += visc * ( cg0*((cg0*std::pow(v[ijk-kk3]-vmean[k-3],2) + cg1*std::pow(v[ijk-kk2]-vmean[k-2],2) + cg2*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg3*std::pow(v[ijk    ]-vmean[k  ],2)) * dzhi4[k-1])
                              + cg1*((cg0*std::pow(v[ijk-kk2]-vmean[k-2],2) + cg1*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg2*std::pow(v[ijk    ]-vmean[k  ],2) + cg3*std::pow(v[ijk+kk1]-vmean[k+1],2)) * dzhi4[k  ])
                              + cg2*((cg0*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg1*std::pow(v[ijk    ]-vmean[k  ],2) + cg2*std::pow(v[ijk+kk1]-vmean[k+1],2) + cg3*std::pow(v[ijk+kk2]-vmean[k+2],2)) * dzhi4[k+1])
                              + cg3*((cg0*std::pow(v[ijk    ]-vmean[k  ],2) + cg1*std::pow(v[ijk+kk1]-vmean[k+1],2) + cg2*std::pow(v[ijk+kk2]-vmean[k+2],2) + cg3*std::pow(v[ijk+kk3]-vmean[k+3],2)) * dzhi4[k+2]) ) * dzi4[k];

        tke_visc[k] += 0.5 * visc * ( cg0*((cg0*std::pow(wx[ijk-kk3],2) + cg1*std::pow(wx[ijk-kk2],2) + cg2*std::pow(wx[ijk-kk1],2) + cg3*std::pow(wx[ijk    ],2)) * dzhi4[k-1])
                                    + cg1*((cg0*std::pow(wx[ijk-kk2],2) + cg1*std::pow(wx[ijk-kk1],2) + cg2*std::pow(wx[ijk    ],2) + cg3*std::pow(wx[ijk+kk1],2)) * dzhi4[k  ])
                                    + cg2*((cg0*std::pow(wx[ijk-kk1],2) + cg1*std::pow(wx[ijk    ],2) + cg2*std::pow(wx[ijk+kk1],2) + cg3*std::pow(wx[ijk+kk2],2)) * dzhi4[k+1])
                                    + cg3*((cg0*std::pow(wx[ijk    ],2) + cg1*std::pow(wx[ijk+kk1],2) + cg2*std::pow(wx[ijk+kk2],2) + cg3*std::pow(wx[ijk+kk3],2)) * dzhi4[k+2]) ) * dzi4[k];
        // tke_visc[k] += 0.5*visc * ( cg0*std::pow(w[ijk-kk1],2) + cg1*std::pow(w[ijk],2) + cg2*std::pow(w[ijk+kk1],2) + cg3*std::pow(w[ijk+kk2],2)) * dzi4[k];
      }
    tke_visc[k] += 0.5*(u2_visc[k] + v2_visc[k]);
  }

  // top boundary
  k = grid->kend-1;
  u2_visc [k] = 0.;
  v2_visc [k] = 0.;
  tke_visc[k] = 0.;

  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      u2_visc[k]  += visc * ( cg0*((cg0*std::pow(u[ijk-kk3]-umean[k-3],2) + cg1*std::pow(u[ijk-kk2]-umean[k-2],2) + cg2*std::pow(u[ijk-kk1]-umean[k-1],2) + cg3*std::pow(u[ijk    ]-umean[k  ],2)) * dzhi4[k-1])
                            + cg1*((cg0*std::pow(u[ijk-kk2]-umean[k-2],2) + cg1*std::pow(u[ijk-kk1]-umean[k-1],2) + cg2*std::pow(u[ijk    ]-umean[k  ],2) + cg3*std::pow(u[ijk+kk1]-umean[k+1],2)) * dzhi4[k  ])
                            + cg2*((cg0*std::pow(u[ijk-kk1]-umean[k-1],2) + cg1*std::pow(u[ijk    ]-umean[k  ],2) + cg2*std::pow(u[ijk+kk1]-umean[k+1],2) + cg3*std::pow(u[ijk+kk2]-umean[k+2],2)) * dzhi4[k+1])
                            + cg3*((tg0*std::pow(u[ijk-kk1]-umean[k-1],2) + tg1*std::pow(u[ijk    ]-umean[k  ],2) + tg2*std::pow(u[ijk+kk1]-umean[k+1],2) + tg3*std::pow(u[ijk+kk2]-umean[k+2],2)) * dzhi4[k+2]) ) * dzi4[k];

      v2_visc[k]  += visc * ( cg0*((cg0*std::pow(v[ijk-kk3]-vmean[k-3],2) + cg1*std::pow(v[ijk-kk2]-vmean[k-2],2) + cg2*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg3*std::pow(v[ijk    ]-vmean[k  ],2)) * dzhi4[k-1])
                            + cg1*((cg0*std::pow(v[ijk-kk2]-vmean[k-2],2) + cg1*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg2*std::pow(v[ijk    ]-vmean[k  ],2) + cg3*std::pow(v[ijk+kk1]-vmean[k+1],2)) * dzhi4[k  ])
                            + cg2*((cg0*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg1*std::pow(v[ijk    ]-vmean[k  ],2) + cg2*std::pow(v[ijk+kk1]-vmean[k+1],2) + cg3*std::pow(v[ijk+kk2]-vmean[k+2],2)) * dzhi4[k+1])
                            + cg3*((tg0*std::pow(v[ijk-kk1]-vmean[k-1],2) + tg1*std::pow(v[ijk    ]-vmean[k  ],2) + tg2*std::pow(v[ijk+kk1]-vmean[k+1],2) + tg3*std::pow(v[ijk+kk2]-vmean[k+2],2)) * dzhi4[k+2]) ) * dzi4[k];

      tke_visc[k] += 0.5 * visc * ( cg0*((cg0*std::pow(wx[ijk-kk3],2) + cg1*std::pow(wx[ijk-kk2],2) + cg2*std::pow(wx[ijk-kk1],2) + cg3*std::pow(wx[ijk    ],2)) * dzhi4[k-1])
                                  + cg1*((cg0*std::pow(wx[ijk-kk2],2) + cg1*std::pow(wx[ijk-kk1],2) + cg2*std::pow(wx[ijk    ],2) + cg3*std::pow(wx[ijk+kk1],2)) * dzhi4[k  ])
                                  + cg2*((cg0*std::pow(wx[ijk-kk1],2) + cg1*std::pow(wx[ijk    ],2) + cg2*std::pow(wx[ijk+kk1],2) + cg3*std::pow(wx[ijk+kk2],2)) * dzhi4[k+1])
                                  + cg3*((tg0*std::pow(wx[ijk-kk1],2) + tg1*std::pow(wx[ijk    ],2) + tg2*std::pow(wx[ijk+kk1],2) + tg3*std::pow(wx[ijk+kk2],2)) * dzhi4[k+2]) ) * dzi4[k];

      // tke_visc[k] += 0.5*visc * ( cg0*std::pow(w[ijk-kk1],2) + cg1*std::pow(w[ijk],2) + cg2*std::pow(w[ijk+kk1],2) + cg3*std::pow(w[ijk+kk2],2)) * dzi4[k];
    }
  tke_visc[k] += 0.5*(u2_visc[k] + v2_visc[k]);

  // calculate the viscous transport of vertical velocity variance
  // bottom
  k = grid->kstart;
  w2_visc[k] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      w2_visc[k] += visc * ( cg0*((bg0*std::pow(w[ijk-kk2],2) + bg1*std::pow(w[ijk-kk1],2) + bg2*std::pow(w[ijk    ],2) + bg3*std::pow(w[ijk+kk1],2)) * dzi4[k-2])
                           + cg1*((cg0*std::pow(w[ijk-kk2],2) + cg1*std::pow(w[ijk-kk1],2) + cg2*std::pow(w[ijk    ],2) + cg3*std::pow(w[ijk+kk1],2)) * dzi4[k-1])
                           + cg2*((cg0*std::pow(w[ijk-kk1],2) + cg1*std::pow(w[ijk    ],2) + cg2*std::pow(w[ijk+kk1],2) + cg3*std::pow(w[ijk+kk2],2)) * dzi4[k  ])
                           + cg3*((cg0*std::pow(w[ijk    ],2) + cg1*std::pow(w[ijk+kk1],2) + cg2*std::pow(w[ijk+kk2],2) + cg3*std::pow(w[ijk+kk3],2)) * dzi4[k+1]) ) * dzhi4[k];
    }

  // interior
  for(int k=grid->kstart+1; k<grid->kend; ++k)
  {
    w2_visc[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        w2_visc[k] += visc * ( cg0*((cg0*std::pow(w[ijk-kk3],2) + cg1*std::pow(w[ijk-kk2],2) + cg2*std::pow(w[ijk-kk1],2) + cg3*std::pow(w[ijk    ],2)) * dzi4[k-2])
                             + cg1*((cg0*std::pow(w[ijk-kk2],2) + cg1*std::pow(w[ijk-kk1],2) + cg2*std::pow(w[ijk    ],2) + cg3*std::pow(w[ijk+kk1],2)) * dzi4[k-1])
                             + cg2*((cg0*std::pow(w[ijk-kk1],2) + cg1*std::pow(w[ijk    ],2) + cg2*std::pow(w[ijk+kk1],2) + cg3*std::pow(w[ijk+kk2],2)) * dzi4[k  ])
                             + cg3*((cg0*std::pow(w[ijk    ],2) + cg1*std::pow(w[ijk+kk1],2) + cg2*std::pow(w[ijk+kk2],2) + cg3*std::pow(w[ijk+kk3],2)) * dzi4[k+1]) ) * dzhi4[k];
      }
  }

  // top boundary
  k = grid->kend;
  w2_visc[k] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      w2_visc[k] += visc * ( cg0*((cg0*std::pow(w[ijk-kk3],2) + cg1*std::pow(w[ijk-kk2],2) + cg2*std::pow(w[ijk-kk1],2) + cg3*std::pow(w[ijk    ],2)) * dzi4[k-2])
                           + cg1*((cg0*std::pow(w[ijk-kk2],2) + cg1*std::pow(w[ijk-kk1],2) + cg2*std::pow(w[ijk    ],2) + cg3*std::pow(w[ijk+kk1],2)) * dzi4[k-1])
                           + cg2*((cg0*std::pow(w[ijk-kk1],2) + cg1*std::pow(w[ijk    ],2) + cg2*std::pow(w[ijk+kk1],2) + cg3*std::pow(w[ijk+kk2],2)) * dzi4[k  ])
                           + cg3*((tg0*std::pow(w[ijk-kk1],2) + tg1*std::pow(w[ijk    ],2) + tg2*std::pow(w[ijk+kk1],2) + tg3*std::pow(w[ijk+kk2],2)) * dzi4[k+1]) ) * dzhi4[k];
    }

  master->sum(u2_visc , grid->kcells);
  master->sum(v2_visc , grid->kcells);
  master->sum(w2_visc , grid->kcells);
  master->sum(tke_visc, grid->kcells);

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    u2_visc [k] /= n;
    v2_visc [k] /= n;
    tke_visc[k] /= n;
  }
  for(int k=grid->kstart; k<grid->kend+1; ++k)
    w2_visc [k] /= n;

  // 6. CALCULATE THE DISSIPATION TERM
  double dxi,dyi;
  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  // bottom boundary
  k = grid->kstart;

  u2_diss [k] = 0.;
  v2_diss [k] = 0.;
  tke_diss[k] = 0.;

  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      u2_diss[k]  -= 2.*visc * (
                       std::pow( ( cg0*((ci0*(u[ijk-ii3]-umean[k]) + ci1*(u[ijk-ii2]-umean[k]) + ci2*(u[ijk-ii1]-umean[k]) + ci3*(u[ijk    ]-umean[k])))
                                 + cg1*((ci0*(u[ijk-ii2]-umean[k]) + ci1*(u[ijk-ii1]-umean[k]) + ci2*(u[ijk    ]-umean[k]) + ci3*(u[ijk+ii1]-umean[k])))
                                 + cg2*((ci0*(u[ijk-ii1]-umean[k]) + ci1*(u[ijk    ]-umean[k]) + ci2*(u[ijk+ii1]-umean[k]) + ci3*(u[ijk+ii2]-umean[k])))
                                 + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+ii1]-umean[k]) + ci2*(u[ijk+ii2]-umean[k]) + ci3*(u[ijk+ii3]-umean[k]))) ) * cgi*dxi, 2)

                     + std::pow( ( cg0*((ci0*(u[ijk-jj3]-umean[k]) + ci1*(u[ijk-jj2]-umean[k]) + ci2*(u[ijk-jj1]-umean[k]) + ci3*(u[ijk    ]-umean[k])))
                                 + cg1*((ci0*(u[ijk-jj2]-umean[k]) + ci1*(u[ijk-jj1]-umean[k]) + ci2*(u[ijk    ]-umean[k]) + ci3*(u[ijk+jj1]-umean[k])))
                                 + cg2*((ci0*(u[ijk-jj1]-umean[k]) + ci1*(u[ijk    ]-umean[k]) + ci2*(u[ijk+jj1]-umean[k]) + ci3*(u[ijk+jj2]-umean[k])))
                                 + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+jj1]-umean[k]) + ci2*(u[ijk+jj2]-umean[k]) + ci3*(u[ijk+jj3]-umean[k]))) ) * cgi*dyi, 2)

                     + std::pow( ( cg0*((bi0*(u[ijk-kk2]-umean[k-2]) + bi1*(u[ijk-kk1]-umean[k-1]) + bi2*(u[ijk    ]-umean[k  ]) + bi3*(u[ijk+kk1]-umean[k+1])))
                                 + cg1*((ci0*(u[ijk-kk2]-umean[k-2]) + ci1*(u[ijk-kk1]-umean[k-1]) + ci2*(u[ijk    ]-umean[k  ]) + ci3*(u[ijk+kk1]-umean[k+1])))
                                 + cg2*((ci0*(u[ijk-kk1]-umean[k-1]) + ci1*(u[ijk    ]-umean[k  ]) + ci2*(u[ijk+kk1]-umean[k+1]) + ci3*(u[ijk+kk2]-umean[k+2])))
                                 + cg3*((ci0*(u[ijk    ]-umean[k  ]) + ci1*(u[ijk+kk1]-umean[k+1]) + ci2*(u[ijk+kk2]-umean[k+2]) + ci3*(u[ijk+kk3]-umean[k+3]))) ) * dzi4[k], 2) );

      v2_diss[k]  -= 2.*visc * (
                       std::pow( ( cg0*((ci0*(v[ijk-ii3]-vmean[k]) + ci1*(v[ijk-ii2]-vmean[k]) + ci2*(v[ijk-ii1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                                 + cg1*((ci0*(v[ijk-ii2]-vmean[k]) + ci1*(v[ijk-ii1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+ii1]-vmean[k])))
                                 + cg2*((ci0*(v[ijk-ii1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+ii1]-vmean[k]) + ci3*(v[ijk+ii2]-vmean[k])))
                                 + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+ii1]-vmean[k]) + ci2*(v[ijk+ii2]-vmean[k]) + ci3*(v[ijk+ii3]-vmean[k]))) ) * cgi*dxi, 2)

                     + std::pow( ( cg0*((ci0*(v[ijk-jj3]-vmean[k]) + ci1*(v[ijk-jj2]-vmean[k]) + ci2*(v[ijk-jj1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                                 + cg1*((ci0*(v[ijk-jj2]-vmean[k]) + ci1*(v[ijk-jj1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+jj1]-vmean[k])))
                                 + cg2*((ci0*(v[ijk-jj1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+jj1]-vmean[k]) + ci3*(v[ijk+jj2]-vmean[k])))
                                 + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+jj1]-vmean[k]) + ci2*(v[ijk+jj2]-vmean[k]) + ci3*(v[ijk+jj3]-vmean[k]))) ) * cgi*dyi, 2)

                     + std::pow( ( cg0*((bi0*(v[ijk-kk2]-vmean[k-2]) + bi1*(v[ijk-kk1]-vmean[k-1]) + bi2*(v[ijk    ]-vmean[k  ]) + bi3*(v[ijk+kk1]-vmean[k+1])))
                                 + cg1*((ci0*(v[ijk-kk2]-vmean[k-2]) + ci1*(v[ijk-kk1]-vmean[k-1]) + ci2*(v[ijk    ]-vmean[k  ]) + ci3*(v[ijk+kk1]-vmean[k+1])))
                                 + cg2*((ci0*(v[ijk-kk1]-vmean[k-1]) + ci1*(v[ijk    ]-vmean[k  ]) + ci2*(v[ijk+kk1]-vmean[k+1]) + ci3*(v[ijk+kk2]-vmean[k+2])))
                                 + cg3*((ci0*(v[ijk    ]-vmean[k  ]) + ci1*(v[ijk+kk1]-vmean[k+1]) + ci2*(v[ijk+kk2]-vmean[k+2]) + ci3*(v[ijk+kk3]-vmean[k+3]))) ) * dzi4[k], 2) );

      tke_diss[k] -= visc * (
                       std::pow( (cg0*w[ijk-ii1] + cg1*w[ijk] + cg2*w[ijk+ii1] + cg3*w[ijk+ii2]) * cgi*dxi, 2)
                     + std::pow( (cg0*w[ijk-jj1] + cg1*w[ijk] + cg2*w[ijk+jj1] + cg3*w[ijk+jj2]) * cgi*dyi, 2)
                     + std::pow( (cg0*w[ijk-kk1] + cg1*w[ijk] + cg2*w[ijk+kk1] + cg3*w[ijk+kk2]) * dzi4[k], 2) );
    }
  tke_diss[k] += 0.5*(u2_diss[k] + v2_diss[k]);

  // interior
  for(int k=grid->kstart+1; k<grid->kend-1; ++k)
  {
    u2_diss [k] = 0.;
    v2_diss [k] = 0.;
    tke_diss[k] = 0.;

    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        u2_diss[k]  -= 2.*visc * (
                         std::pow( ( cg0*((ci0*(u[ijk-ii3]-umean[k]) + ci1*(u[ijk-ii2]-umean[k]) + ci2*(u[ijk-ii1]-umean[k]) + ci3*(u[ijk    ]-umean[k])))
                                   + cg1*((ci0*(u[ijk-ii2]-umean[k]) + ci1*(u[ijk-ii1]-umean[k]) + ci2*(u[ijk    ]-umean[k]) + ci3*(u[ijk+ii1]-umean[k])))
                                   + cg2*((ci0*(u[ijk-ii1]-umean[k]) + ci1*(u[ijk    ]-umean[k]) + ci2*(u[ijk+ii1]-umean[k]) + ci3*(u[ijk+ii2]-umean[k])))
                                   + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+ii1]-umean[k]) + ci2*(u[ijk+ii2]-umean[k]) + ci3*(u[ijk+ii3]-umean[k]))) ) * cgi*dxi, 2)

                       + std::pow( ( cg0*((ci0*(u[ijk-jj3]-umean[k]) + ci1*(u[ijk-jj2]-umean[k]) + ci2*(u[ijk-jj1]-umean[k]) + ci3*(u[ijk    ]-umean[k])))
                                   + cg1*((ci0*(u[ijk-jj2]-umean[k]) + ci1*(u[ijk-jj1]-umean[k]) + ci2*(u[ijk    ]-umean[k]) + ci3*(u[ijk+jj1]-umean[k])))
                                   + cg2*((ci0*(u[ijk-jj1]-umean[k]) + ci1*(u[ijk    ]-umean[k]) + ci2*(u[ijk+jj1]-umean[k]) + ci3*(u[ijk+jj2]-umean[k])))
                                   + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+jj1]-umean[k]) + ci2*(u[ijk+jj2]-umean[k]) + ci3*(u[ijk+jj3]-umean[k]))) ) * cgi*dyi, 2)

                       + std::pow( ( cg0*((ci0*(u[ijk-kk3]-umean[k-3]) + ci1*(u[ijk-kk2]-umean[k-2]) + ci2*(u[ijk-kk1]-umean[k-1]) + ci3*(u[ijk    ]-umean[k  ])))
                                   + cg1*((ci0*(u[ijk-kk2]-umean[k-2]) + ci1*(u[ijk-kk1]-umean[k-1]) + ci2*(u[ijk    ]-umean[k  ]) + ci3*(u[ijk+kk1]-umean[k+1])))
                                   + cg2*((ci0*(u[ijk-kk1]-umean[k-1]) + ci1*(u[ijk    ]-umean[k  ]) + ci2*(u[ijk+kk1]-umean[k+1]) + ci3*(u[ijk+kk2]-umean[k+2])))
                                   + cg3*((ci0*(u[ijk    ]-umean[k  ]) + ci1*(u[ijk+kk1]-umean[k+1]) + ci2*(u[ijk+kk2]-umean[k+2]) + ci3*(u[ijk+kk3]-umean[k+3]))) ) * dzi4[k], 2) );

        v2_diss[k]  -= 2.*visc * (
                         std::pow( ( cg0*((ci0*(v[ijk-ii3]-vmean[k]) + ci1*(v[ijk-ii2]-vmean[k]) + ci2*(v[ijk-ii1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                                   + cg1*((ci0*(v[ijk-ii2]-vmean[k]) + ci1*(v[ijk-ii1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+ii1]-vmean[k])))
                                   + cg2*((ci0*(v[ijk-ii1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+ii1]-vmean[k]) + ci3*(v[ijk+ii2]-vmean[k])))
                                   + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+ii1]-vmean[k]) + ci2*(v[ijk+ii2]-vmean[k]) + ci3*(v[ijk+ii3]-vmean[k]))) ) * cgi*dxi, 2)

                       + std::pow( ( cg0*((ci0*(v[ijk-jj3]-vmean[k]) + ci1*(v[ijk-jj2]-vmean[k]) + ci2*(v[ijk-jj1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                                   + cg1*((ci0*(v[ijk-jj2]-vmean[k]) + ci1*(v[ijk-jj1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+jj1]-vmean[k])))
                                   + cg2*((ci0*(v[ijk-jj1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+jj1]-vmean[k]) + ci3*(v[ijk+jj2]-vmean[k])))
                                   + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+jj1]-vmean[k]) + ci2*(v[ijk+jj2]-vmean[k]) + ci3*(v[ijk+jj3]-vmean[k]))) ) * cgi*dyi, 2)

                       + std::pow( ( cg0*((ci0*(v[ijk-kk3]-vmean[k-3]) + ci1*(v[ijk-kk2]-vmean[k-2]) + ci2*(v[ijk-kk1]-vmean[k-1]) + ci3*(v[ijk    ]-vmean[k  ])))
                                   + cg1*((ci0*(v[ijk-kk2]-vmean[k-2]) + ci1*(v[ijk-kk1]-vmean[k-1]) + ci2*(v[ijk    ]-vmean[k  ]) + ci3*(v[ijk+kk1]-vmean[k+1])))
                                   + cg2*((ci0*(v[ijk-kk1]-vmean[k-1]) + ci1*(v[ijk    ]-vmean[k  ]) + ci2*(v[ijk+kk1]-vmean[k+1]) + ci3*(v[ijk+kk2]-vmean[k+2])))
                                   + cg3*((ci0*(v[ijk    ]-vmean[k  ]) + ci1*(v[ijk+kk1]-vmean[k+1]) + ci2*(v[ijk+kk2]-vmean[k+2]) + ci3*(v[ijk+kk3]-vmean[k+3]))) ) * dzi4[k], 2) );

        tke_diss[k] -= visc * (
                         std::pow( (cg0*w[ijk-ii1] + cg1*w[ijk] + cg2*w[ijk+ii1] + cg3*w[ijk+ii2]) * cgi*dxi, 2)
                       + std::pow( (cg0*w[ijk-jj1] + cg1*w[ijk] + cg2*w[ijk+jj1] + cg3*w[ijk+jj2]) * cgi*dyi, 2)
                       + std::pow( (cg0*w[ijk-kk1] + cg1*w[ijk] + cg2*w[ijk+kk1] + cg3*w[ijk+kk2]) * dzi4[k], 2) );
      }
    tke_diss[k] += 0.5*(u2_diss[k] + v2_diss[k]);
  }

  // top boundary
  k = grid->kend-1;

  u2_diss [k] = 0.;
  v2_diss [k] = 0.;
  tke_diss[k] = 0.;

  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk  = i + j*jj1 + k*kk1;
      u2_diss[k]  -= 2.*visc * (
                       std::pow( ( cg0*((ci0*(u[ijk-ii3]-umean[k]) + ci1*(u[ijk-ii2]-umean[k]) + ci2*(u[ijk-ii1]-umean[k]) + ci3*(u[ijk    ]-umean[k])))
                                 + cg1*((ci0*(u[ijk-ii2]-umean[k]) + ci1*(u[ijk-ii1]-umean[k]) + ci2*(u[ijk    ]-umean[k]) + ci3*(u[ijk+ii1]-umean[k])))
                                 + cg2*((ci0*(u[ijk-ii1]-umean[k]) + ci1*(u[ijk    ]-umean[k]) + ci2*(u[ijk+ii1]-umean[k]) + ci3*(u[ijk+ii2]-umean[k])))
                                 + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+ii1]-umean[k]) + ci2*(u[ijk+ii2]-umean[k]) + ci3*(u[ijk+ii3]-umean[k]))) ) * cgi*dxi, 2)

                     + std::pow( ( cg0*((ci0*(u[ijk-jj3]-umean[k]) + ci1*(u[ijk-jj2]-umean[k]) + ci2*(u[ijk-jj1]-umean[k]) + ci3*(u[ijk    ]-umean[k])))
                                 + cg1*((ci0*(u[ijk-jj2]-umean[k]) + ci1*(u[ijk-jj1]-umean[k]) + ci2*(u[ijk    ]-umean[k]) + ci3*(u[ijk+jj1]-umean[k])))
                                 + cg2*((ci0*(u[ijk-jj1]-umean[k]) + ci1*(u[ijk    ]-umean[k]) + ci2*(u[ijk+jj1]-umean[k]) + ci3*(u[ijk+jj2]-umean[k])))
                                 + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+jj1]-umean[k]) + ci2*(u[ijk+jj2]-umean[k]) + ci3*(u[ijk+jj3]-umean[k]))) ) * cgi*dyi, 2)

                     + std::pow( ( cg0*((ci0*(u[ijk-kk3]-umean[k-3]) + ci1*(u[ijk-kk2]-umean[k-2]) + ci2*(u[ijk-kk1]-umean[k-1]) + ci3*(u[ijk    ]-umean[k  ])))
                                 + cg1*((ci0*(u[ijk-kk2]-umean[k-2]) + ci1*(u[ijk-kk1]-umean[k-1]) + ci2*(u[ijk    ]-umean[k  ]) + ci3*(u[ijk+kk1]-umean[k+1])))
                                 + cg2*((ci0*(u[ijk-kk1]-umean[k-1]) + ci1*(u[ijk    ]-umean[k  ]) + ci2*(u[ijk+kk1]-umean[k+1]) + ci3*(u[ijk+kk2]-umean[k+2])))
                                 + cg3*((ti0*(u[ijk-kk1]-umean[k-1]) + ti1*(u[ijk    ]-umean[k  ]) + ti2*(u[ijk+kk1]-umean[k+1]) + ti3*(u[ijk+kk2]-umean[k+2]))) ) * dzi4[k], 2) );

      v2_diss[k]  -= 2.*visc * (
                       std::pow( ( cg0*((ci0*(v[ijk-ii3]-vmean[k]) + ci1*(v[ijk-ii2]-vmean[k]) + ci2*(v[ijk-ii1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                                 + cg1*((ci0*(v[ijk-ii2]-vmean[k]) + ci1*(v[ijk-ii1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+ii1]-vmean[k])))
                                 + cg2*((ci0*(v[ijk-ii1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+ii1]-vmean[k]) + ci3*(v[ijk+ii2]-vmean[k])))
                                 + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+ii1]-vmean[k]) + ci2*(v[ijk+ii2]-vmean[k]) + ci3*(v[ijk+ii3]-vmean[k]))) ) * cgi*dxi, 2)

                     + std::pow( ( cg0*((ci0*(v[ijk-jj3]-vmean[k]) + ci1*(v[ijk-jj2]-vmean[k]) + ci2*(v[ijk-jj1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                                 + cg1*((ci0*(v[ijk-jj2]-vmean[k]) + ci1*(v[ijk-jj1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+jj1]-vmean[k])))
                                 + cg2*((ci0*(v[ijk-jj1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+jj1]-vmean[k]) + ci3*(v[ijk+jj2]-vmean[k])))
                                 + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+jj1]-vmean[k]) + ci2*(v[ijk+jj2]-vmean[k]) + ci3*(v[ijk+jj3]-vmean[k]))) ) * cgi*dyi, 2)

                     + std::pow( ( cg0*((ci0*(v[ijk-kk3]-vmean[k-3]) + ci1*(v[ijk-kk2]-vmean[k-2]) + ci2*(v[ijk-kk1]-vmean[k-1]) + ci3*(v[ijk    ]-vmean[k  ])))
                                 + cg1*((ci0*(v[ijk-kk2]-vmean[k-2]) + ci1*(v[ijk-kk1]-vmean[k-1]) + ci2*(v[ijk    ]-vmean[k  ]) + ci3*(v[ijk+kk1]-vmean[k+1])))
                                 + cg2*((ci0*(v[ijk-kk1]-vmean[k-1]) + ci1*(v[ijk    ]-vmean[k  ]) + ci2*(v[ijk+kk1]-vmean[k+1]) + ci3*(v[ijk+kk2]-vmean[k+2])))
                                 + cg3*((ti0*(v[ijk-kk1]-vmean[k-1]) + ti1*(v[ijk    ]-vmean[k  ]) + ti2*(v[ijk+kk1]-vmean[k+1]) + ti3*(v[ijk+kk2]-vmean[k+2]))) ) * dzi4[k], 2) );

      tke_diss[k] -= visc * (
                       std::pow( (cg0*w[ijk-ii1] + cg1*w[ijk] + cg2*w[ijk+ii1] + cg3*w[ijk+ii2]) * cgi*dxi, 2)
                     + std::pow( (cg0*w[ijk-jj1] + cg1*w[ijk] + cg2*w[ijk+jj1] + cg3*w[ijk+jj2]) * cgi*dyi, 2)
                     + std::pow( (cg0*w[ijk-kk1] + cg1*w[ijk] + cg2*w[ijk+kk1] + cg3*w[ijk+kk2]) * dzi4[k], 2) );
    }
  tke_diss[k] += 0.5*(u2_diss[k] + v2_diss[k]);

  // calculate the w2 budget term
  for(int k=grid->kstart+1; k<grid->kend; ++k)
  {
    w2_diss[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        w2_diss[k]  -= 2.*visc * (
                         std::pow( ( cg0*(ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ])
                                   + cg1*(ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1])
                                   + cg2*(ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2])
                                   + cg3*(ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3]) ) * cgi*dxi, 2)

                       + std::pow( ( cg0*(ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ])
                                   + cg1*(ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1])
                                   + cg2*(ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2])
                                   + cg3*(ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3]) ) * cgi*dyi, 2)

                       + std::pow( ( cg0*(ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ])
                                   + cg1*(ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1])
                                   + cg2*(ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])
                                   + cg3*(ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3]) ) * dzhi4[k], 2) );
      }
  }

  master->sum(u2_diss , grid->kcells);
  master->sum(v2_diss , grid->kcells);
  master->sum(w2_diss , grid->kcells);
  master->sum(tke_diss, grid->kcells);

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    u2_diss [k] /= n;
    v2_diss [k] /= n;
    tke_diss[k] /= n;
  }

  for(int k=grid->kstart; k<grid->kend; ++k)
    w2_diss [k] /= n;

  // 7. CALCULATE THE PRESSURE REDISTRIBUTION TERM
  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    u2_rdstr [k] = 0.;
    v2_rdstr [k] = 0.;

    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        u2_rdstr [k] += 2.*(ci0*p[ijk-ii2] + ci1*p[ijk-ii1] + ci2*p[ijk] + ci3*p[ijk+ii1])*
                        ( cg0*((ci0*(u[ijk-ii3]-umean[k]) + ci1*(u[ijk-ii2]-umean[k]) + ci2*(u[ijk-ii1]-umean[k]) + ci3*(u[ijk    ]-umean[k])))
                        + cg1*((ci0*(u[ijk-ii2]-umean[k]) + ci1*(u[ijk-ii1]-umean[k]) + ci2*(u[ijk    ]-umean[k]) + ci3*(u[ijk+ii1]-umean[k])))
                        + cg2*((ci0*(u[ijk-ii1]-umean[k]) + ci1*(u[ijk    ]-umean[k]) + ci2*(u[ijk+ii1]-umean[k]) + ci3*(u[ijk+ii2]-umean[k])))
                        + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+ii1]-umean[k]) + ci2*(u[ijk+ii2]-umean[k]) + ci3*(u[ijk+ii3]-umean[k]))) ) * cgi*dxi;
        v2_rdstr [k] += 2.*(ci0*p[ijk-jj2] + ci1*p[ijk-jj1] + ci2*p[ijk] + ci3*p[ijk+jj1])*
                        ( cg0*((ci0*(v[ijk-jj3]-vmean[k]) + ci1*(v[ijk-jj2]-vmean[k]) + ci2*(v[ijk-jj1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                        + cg1*((ci0*(v[ijk-jj2]-vmean[k]) + ci1*(v[ijk-jj1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+jj1]-vmean[k])))
                        + cg2*((ci0*(v[ijk-jj1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+jj1]-vmean[k]) + ci3*(v[ijk+jj2]-vmean[k])))
                        + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+jj1]-vmean[k]) + ci2*(v[ijk+jj2]-vmean[k]) + ci3*(v[ijk+jj3]-vmean[k]))) ) * cgi*dyi;
      }
  }

  for(int k=grid->kstart+1; k<grid->kend; ++k)
  {
    w2_rdstr[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        w2_rdstr[k] += 2.*(ci0*p[ijk-kk2] + ci1*p[ijk-kk1] + ci2*p[ijk] + ci3*p[ijk+kk1])*
                       ( cg0*(ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ])
                       + cg1*(ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1])
                       + cg2*(ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])
                       + cg3*(ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3]) ) * dzhi4[k];
      }
  }

  master->sum(u2_rdstr , grid->kcells);
  master->sum(v2_rdstr , grid->kcells);
  master->sum(w2_rdstr , grid->kcells);

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    u2_rdstr[k] /= n;
    v2_rdstr[k] /= n;
    w2_rdstr[k] /= n;
  }

  return 0;
}

int Budget::calctkebudget_buoy(double * restrict w, double * restrict b,
                                double * restrict w2_buoy, double * restrict tke_buoy)
{
  int ijk,jj1,kk1,kk2;

  jj1 = 1*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  double n = grid->imax*grid->jmax;

  // calculate the buoyancy term
  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    w2_buoy [k] = 0.;
    tke_buoy[k] = 0.;

    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        tke_buoy[k] += (ci0*w[ijk-kk1] + ci1*w[ijk] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])*b[ijk];
      }
  }
  for(int k=grid->kstart+1; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        w2_buoy[k] += 2.*(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk] + ci3*b[ijk+kk1])*w[ijk];
      }

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    w2_buoy [k] /= n;
    tke_buoy[k] /= n;
  }

  grid->getProf(w2_buoy , grid->kcells);
  grid->getProf(tke_buoy, grid->kcells);

  return 0;
}

int Budget::calcpe(double * restrict b, double * restrict zsort, double * restrict zsortbot, double * restrict zsorttop,
                    double * restrict z,
                    double * restrict bsort,
                    double * restrict pe_total, double * restrict pe_avail, double * restrict pe_bg,
                    double * restrict zsortprof)
{
  int ij,ijk,jj,kk1,kk2,kstart,kend;

  jj = grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;
  kstart = grid->kstart;
  kend = grid->kend;

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    pe_total[k] = 0;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk1;
        pe_total[k] -= b[ijk] * z[k];
      }
  }

  master->sum(pe_total, grid->kcells);

  int n = grid->itot*grid->jtot;
  for(int k=grid->kstart; k<grid->kend; ++k)
    pe_total[k] /= n;

  // now find out the available potential energy
  // int ks;
  double zsortval;
  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    zsortprof[k] = 0.;
    pe_bg    [k] = 0.;
    pe_avail [k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk1;
        /*
        ks  = k;
        if(b[ijk] > bsort[k])
        {
          while(b[ijk] > bsort[ks] && ks < grid->kend-1)
            ++ks;

          // linearly interpolate the height
          zsortval = z[ks-1] + (b[ijk]-bsort[ks-1])/(bsort[ks]-bsort[ks-1]) * (z[ks]-z[ks-1]);

        }
        else if(b[ijk] < bsort[k])
        {
          while(b[ijk] < bsort[ks] && ks > grid->kstart)
            --ks;

          // linearly interpolate the height
          zsortval = z[ks] + (b[ijk]-bsort[ks])/(bsort[ks+1]-bsort[ks]) * (z[ks+1]-z[ks]);
        }
        else
          zsortval = z[ks];
          */
        zsortval = calczsort(b[ijk], bsort, z, k);
        zsort[ijk] = zsortval;

        zsortprof[k] += zsortval;
        pe_bg    [k] -= zsortval*b[ijk];
        pe_avail [k] -= (z[k]-zsortval)*b[ijk];
      }
  }

  master->sum(zsortprof, grid->kcells);
  master->sum(pe_bg    , grid->kcells);
  master->sum(pe_avail , grid->kcells);

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    zsortprof[k] /= n;
    pe_bg    [k] /= n;
    pe_avail [k] /= n;
  }

  // now, calculate the boundary conditions for zsort
  // bottom bc
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk1;
      zsortbot[ij] = calczsort(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk] + ci3*b[ijk+kk1], bsort, z, kstart);
    }

  // top bc
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + (kend-1)*kk1;
      zsorttop[ij] = calczsort(ci0*b[ijk-kk1] + ci1*b[ijk] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2], bsort, z, kend-1);
    }

  // calculate the ghost cells at the bottom
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk1;
      zsort[ijk-kk1] = (8./3.)*zsortbot[ij] - 2.*zsort[ijk] + (1./3.)*zsort[ijk+kk1];
      zsort[ijk-kk2] = 8.*zsortbot[ij] - 9.*zsort[ijk] + 2.*zsort[ijk+kk1];
    }

  // calculate the ghost cells at the top
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + (kend-1)*kk1;
      zsort[ijk+kk1] = (8./3.)*zsorttop[ij] - 2.*zsort[ijk] + (1./3.)*zsort[ijk-kk1];
      zsort[ijk+kk2] = 8.*zsorttop[ij] - 9.*zsort[ijk] + 2.*zsort[ijk-kk1];
    }

  return 0;
}

double Budget::calczsort(double b, double * restrict bsort, double * restrict z, int k)
{
  double zsortval;
  int ks = k;

  if(b > bsort[k])
  {
    while(b > bsort[ks] && ks < grid->kend-1)
      ++ks;

    // linearly interpolate the height
    zsortval = z[ks-1] + (b-bsort[ks-1])/(bsort[ks]-bsort[ks-1]) * (z[ks]-z[ks-1]);

  }
  else if(b < bsort[k])
  {
    while(b < bsort[ks] && ks > grid->kstart)
      --ks;

    // linearly interpolate the height
    zsortval = z[ks] + (b-bsort[ks])/(bsort[ks+1]-bsort[ks]) * (z[ks+1]-z[ks]);
  }
  else
    zsortval = z[ks];

  return zsortval;
}

double Budget::calcdzstardb(double b, double * restrict bsort, double * restrict z)
{
  // start the iteration below the grid to make sure not to miss values below the first full level
  int k = grid->kstart-1;
  while(bsort[k+1] < b && k < grid->kend)
    ++k;

  // our required value is in between bsort[k] and bsort[k+1]
  // calculate a spline of the form zstar(b) = a*zstar(b_k) + b*zstar(b_(k+1)) + c * zstar''(b_k) + d * zstar''(b_(k+1))
  double ca, cb;//, cc, cd;
  ca = (bsort[k+1]-b) / (bsort[k+1]-bsort[k]);
  cb = 1.-ca;
  // cc = (1./6.) * (ca*ca*ca-ca) * std::pow(bsort[k+1]-bsort[k],2);
  // cd = (1./6.) * (cb*cb*cb-cb) * std::pow(bsort[k+1]-bsort[k],2);

  // calculate the second derivatives using second order accuracy since the profile is very smooth
  double d2zstarb2k, d2zstarb2kp;
  d2zstarb2k  = 2.*((z[k+1]-z[k  ])/(bsort[k+1]-bsort[k  ]) - (z[k  ]-z[k-1])/(bsort[k  ]-bsort[k-1])) / (bsort[k+1]-bsort[k-1]);
  d2zstarb2kp = 2.*((z[k+2]-z[k+1])/(bsort[k+2]-bsort[k+1]) - (z[k+1]-z[k  ])/(bsort[k+1]-bsort[k  ])) / (bsort[k+2]-bsort[k  ]);

  // std::printf("CvH %E, %E, %E, %E, %E, %E\n", bsort[k], b, bsort[k+1], z[k], ca*z[k] + cb*z[k+1], z[k+1]);

  // the derivative is computed according to:
  // dzstar/db = (zstar[k+1]-zstar[k])/(b[k+1]-b[k]) 
  //           - (3*a^2-1)/6 * (b[k+1]-b[k])*zstar''(b_k) 
  //           + (3*b^2-1)/6 * (b[k+1]-b[k])*zstar''(b_(k+1))

  double dzstardb = (z[k+1]-z[k]) / (bsort[k+1]-bsort[k]) 
                  - (3.*ca*ca-1.)/6. * (bsort[k+1]-bsort[k])*d2zstarb2k
                  + (3.*cb*cb-1.)/6. * (bsort[k+1]-bsort[k])*d2zstarb2kp;

  return dzstardb;
}


int Budget::calcpebudget(double * restrict w, double * restrict b, double * restrict bz, double * restrict bztop,
                          double * restrict pe_turb, double * restrict pe_visc, double * restrict pe_bous,
                          double * restrict z, double * restrict zh, double * restrict dzi4, double * restrict dzhi4,
                          double visc)
{
  int ij,ijk,jj1,kk1,kk2,kk3;
  int kstart,kend;
  double zsize;

  jj1 = 1*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;
  kk3 = 3*grid->ijcells;
  kstart = grid->kstart;
  kend = grid->kend;
  zsize = grid->zsize;

  // first, calculate the Boussinesq term (kappa*db/dz). Here bz contains the buoyancy and not the PE yet
  // bottom boundary
  pe_bous[kstart] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk = i + j*jj1 + kstart*kk1;
      pe_bous[kstart] += visc * ( cg0*(bi0*b[ijk-kk2] + bi1*b[ijk-kk1] + bi2*b[ijk    ] + bi3*b[ijk+kk1])
                                + cg1*(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1])
                                + cg2*(ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2])
                                + cg3*(ci0*b[ijk    ] + ci1*b[ijk+kk1] + ci2*b[ijk+kk2] + ci3*b[ijk+kk3]) )
                                * dzi4[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; ++k)
  {
    pe_bous[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj1 + k*kk1;
        pe_bous[k] += visc * ( cg0*(ci0*b[ijk-kk3] + ci1*b[ijk-kk2] + ci2*b[ijk-kk1] + ci3*b[ijk    ])
                             + cg1*(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1])
                             + cg2*(ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2])
                             + cg3*(ci0*b[ijk    ] + ci1*b[ijk+kk1] + ci2*b[ijk+kk2] + ci3*b[ijk+kk3]) )
                             * dzi4[k];
      }
  }

  // top boundary
  pe_bous[kend-1] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      pe_bous[kend-1] += visc * ( cg0*(ci0*b[ijk-kk3] + ci1*b[ijk-kk2] + ci2*b[ijk-kk1] + ci3*b[ijk    ])
                                + cg1*(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1])
                                + cg2*(ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2])
                                + cg3*(ti0*b[ijk-kk1] + ti1*b[ijk    ] + ti2*b[ijk+kk1] + ti3*b[ijk+kk2]) )
                                * dzi4[kend-1];
    }

  // now, convert the buoyancy field into a potential energy field
  // first, before destroying the field, calculate the potential energy at the top
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj1;
      ijk = i + j*jj1 + (kend-1)*kk1;
      bztop[ij] = -zsize*(ci0*b[ijk-kk1] + ci1*b[ijk] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2]);
    }

  // calculate the potential energy
  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj1 + k*kk1;
        bz[ijk] = -b[ijk] * z[k];
      }

  // calculate the ghost cells at the bottom, making use of the fact that bz = 0
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj1;
      ijk = i + j*jj1 + kstart*kk1;
      bz[ijk-kk1] = - 2.*bz[ijk] + (1./3.)*bz[ijk+kk1];
      bz[ijk-kk2] = - 9.*bz[ijk] + 2.*bz[ijk+kk1];
    }

  // calculate the ghost cells at the top
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj1;
      ijk = i + j*jj1 + (kend-1)*kk1;
      bz[ijk+kk1] = (8./3.)*bztop[ij] - 2.*bz[ijk] + (1./3.)*bz[ijk-kk1];
      bz[ijk+kk2] = 8.*bztop[ij] - 9.*bz[ijk] + 2.*bz[ijk-kk1];
    }

  // calculate the advective transport term
  // bottom boundary
  pe_turb[kstart] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk = i + j*jj1 + kstart*kk1;
      pe_turb[kstart] -= ( cg0*(w[ijk-kk1] * (bi0*bz[ijk-kk2] + bi1*bz[ijk-kk1] + bi2*bz[ijk    ] + bi3*bz[ijk+kk1]))
                         + cg1*(w[ijk    ] * (ci0*bz[ijk-kk2] + ci1*bz[ijk-kk1] + ci2*bz[ijk    ] + ci3*bz[ijk+kk1]))
                         + cg2*(w[ijk+kk1] * (ci0*bz[ijk-kk1] + ci1*bz[ijk    ] + ci2*bz[ijk+kk1] + ci3*bz[ijk+kk2]))
                         + cg3*(w[ijk+kk2] * (ci0*bz[ijk    ] + ci1*bz[ijk+kk1] + ci2*bz[ijk+kk2] + ci3*bz[ijk+kk3])) )
                         * dzi4[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; ++k)
  {
    pe_turb[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj1 + k*kk1;
        pe_turb[k] -= ( cg0*(w[ijk-kk1] * (ci0*bz[ijk-kk3] + ci1*bz[ijk-kk2] + ci2*bz[ijk-kk1] + ci3*bz[ijk    ]))
                      + cg1*(w[ijk    ] * (ci0*bz[ijk-kk2] + ci1*bz[ijk-kk1] + ci2*bz[ijk    ] + ci3*bz[ijk+kk1]))
                      + cg2*(w[ijk+kk1] * (ci0*bz[ijk-kk1] + ci1*bz[ijk    ] + ci2*bz[ijk+kk1] + ci3*bz[ijk+kk2]))
                      + cg3*(w[ijk+kk2] * (ci0*bz[ijk    ] + ci1*bz[ijk+kk1] + ci2*bz[ijk+kk2] + ci3*bz[ijk+kk3])) )
                      * dzi4[k];
      }
  }

  // top boundary
  pe_turb[kend-1] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      pe_turb[kend-1] -= ( cg0*(w[ijk-kk1] * (ci0*bz[ijk-kk3] + ci1*bz[ijk-kk2] + ci2*bz[ijk-kk1] + ci3*bz[ijk    ]))
                         + cg1*(w[ijk    ] * (ci0*bz[ijk-kk2] + ci1*bz[ijk-kk1] + ci2*bz[ijk    ] + ci3*bz[ijk+kk1]))
                         + cg2*(w[ijk+kk1] * (ci0*bz[ijk-kk1] + ci1*bz[ijk    ] + ci2*bz[ijk+kk1] + ci3*bz[ijk+kk2]))
                         + cg3*(w[ijk+kk2] * (ti0*bz[ijk-kk1] + ti1*bz[ijk    ] + ti2*bz[ijk+kk1] + ti3*bz[ijk+kk2])) )
                         * dzi4[kend-1];
    }

  // calculate the diffusion of potential energy
  // bottom boundary
  pe_visc[kstart] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk = i + j*jj1 + kstart*kk1;
      pe_visc[kstart] -= visc * ( cg0*zh[kstart-1]*(bg0*b[ijk-kk2] + bg1*b[ijk-kk1] + bg2*b[ijk    ] + bg3*b[ijk+kk1]) * dzhi4[kstart-1]
                                + cg1*zh[kstart  ]*(cg0*b[ijk-kk2] + cg1*b[ijk-kk1] + cg2*b[ijk    ] + cg3*b[ijk+kk1]) * dzhi4[kstart  ]
                                + cg2*zh[kstart+1]*(cg0*b[ijk-kk1] + cg1*b[ijk    ] + cg2*b[ijk+kk1] + cg3*b[ijk+kk2]) * dzhi4[kstart+1]
                                + cg3*zh[kstart+2]*(cg0*b[ijk    ] + cg1*b[ijk+kk1] + cg2*b[ijk+kk2] + cg3*b[ijk+kk3]) * dzhi4[kstart+2] )
                                * dzi4[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; ++k)
  {
    pe_visc[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj1 + k*kk1;
        pe_visc[k] -= visc * ( cg0*zh[k-1]*(cg0*b[ijk-kk3] + cg1*b[ijk-kk2] + cg2*b[ijk-kk1] + cg3*b[ijk    ]) * dzhi4[k-1]
                             + cg1*zh[k  ]*(cg0*b[ijk-kk2] + cg1*b[ijk-kk1] + cg2*b[ijk    ] + cg3*b[ijk+kk1]) * dzhi4[k  ]
                             + cg2*zh[k+1]*(cg0*b[ijk-kk1] + cg1*b[ijk    ] + cg2*b[ijk+kk1] + cg3*b[ijk+kk2]) * dzhi4[k+1]
                             + cg3*zh[k+2]*(cg0*b[ijk    ] + cg1*b[ijk+kk1] + cg2*b[ijk+kk2] + cg3*b[ijk+kk3]) * dzhi4[k+2] )
                             * dzi4[k];
      }
  }

  // top boundary
  pe_visc[kend-1] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      pe_visc[kend-1] -= visc * ( cg0*zh[kend-2]*(cg0*b[ijk-kk3] + cg1*b[ijk-kk2] + cg2*b[ijk-kk1] + cg3*b[ijk    ]) * dzhi4[kend-2]
                                + cg1*zh[kend-1]*(cg0*b[ijk-kk2] + cg1*b[ijk-kk1] + cg2*b[ijk    ] + cg3*b[ijk+kk1]) * dzhi4[kend-1]
                                + cg2*zh[kend  ]*(cg0*b[ijk-kk1] + cg1*b[ijk    ] + cg2*b[ijk+kk1] + cg3*b[ijk+kk2]) * dzhi4[kend  ]
                                + cg3*zh[kend+1]*(tg0*b[ijk-kk1] + tg1*b[ijk    ] + tg2*b[ijk+kk1] + tg3*b[ijk+kk2]) * dzhi4[kend+1] )
                                * dzi4[kend-1];
    }

  master->sum(pe_turb, grid->kcells);
  master->sum(pe_visc, grid->kcells);
  master->sum(pe_bous, grid->kcells);

  int n = grid->itot*grid->jtot;
  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    pe_turb[k] /= n;
    pe_visc[k] /= n;
    pe_bous[k] /= n;
  }

  return 0;
}

int Budget::calcbpebudget(double * restrict w, double * restrict b, double * restrict bz, double * restrict bzbot, double * restrict bztop,
                           double * restrict bpe_turb, double * restrict bpe_visc, double * restrict bpe_diss,
                           double * restrict bsort,
                           double * restrict z, double * restrict dzi4, double * restrict dzhi4,
                           double visc)
{
  int ijk,ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  int kstart,kend;
  double dzstardb;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;
  kk3 = 3*grid->ijcells;
  kstart = grid->kstart;
  kend = grid->kend;

  // calculate the diffusion of potential energy
  // bottom boundary
  bpe_visc[kstart] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk = i + j*jj1 + kstart*kk1;
      bpe_visc[kstart] -= visc *
                            ( cg0*(bg0*b [ijk-kk2] + bg1*b [ijk-kk1] + bg2*b [ijk    ] + bg3*b [ijk+kk1]) * dzhi4[kstart-1]
                                 *(bi0*bz[ijk-kk2] + bi1*bz[ijk-kk1] + bi2*bz[ijk    ] + bi3*bz[ijk+kk1])
                            + cg1*(cg0*b [ijk-kk2] + cg1*b [ijk-kk1] + cg2*b [ijk    ] + cg3*b [ijk+kk1]) * dzhi4[kstart  ]
                                 *(ci0*bz[ijk-kk2] + ci1*bz[ijk-kk1] + ci2*bz[ijk    ] + ci3*bz[ijk+kk1])
                            + cg2*(cg0*b [ijk-kk1] + cg1*b [ijk    ] + cg2*b [ijk+kk1] + cg3*b [ijk+kk2]) * dzhi4[kstart+1]
                                 *(ci0*bz[ijk-kk1] + ci1*bz[ijk    ] + ci2*bz[ijk+kk1] + ci3*bz[ijk+kk2])
                            + cg3*(cg0*b [ijk    ] + cg1*b [ijk+kk1] + cg2*b [ijk+kk2] + cg3*b [ijk+kk3]) * dzhi4[kstart+2]
                                 *(ci0*bz[ijk    ] + ci1*bz[ijk+kk1] + ci2*bz[ijk+kk2] + ci3*bz[ijk+kk3]) )
                            * dzi4[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; ++k)
  {
    bpe_visc[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj1 + k*kk1;
        bpe_visc[k] -= visc *
                         ( cg0*(cg0*b [ijk-kk3] + cg1*b [ijk-kk2] + cg2*b [ijk-kk1] + cg3*b [ijk    ]) * dzhi4[k-1]
                              *(ci0*bz[ijk-kk3] + ci1*bz[ijk-kk2] + ci2*bz[ijk-kk1] + ci3*bz[ijk    ])
                         + cg1*(cg0*b [ijk-kk2] + cg1*b [ijk-kk1] + cg2*b [ijk    ] + cg3*b [ijk+kk1]) * dzhi4[k  ]
                              *(ci0*bz[ijk-kk2] + ci1*bz[ijk-kk1] + ci2*bz[ijk    ] + ci3*bz[ijk+kk1])
                         + cg2*(cg0*b [ijk-kk1] + cg1*b [ijk    ] + cg2*b [ijk+kk1] + cg3*b [ijk+kk2]) * dzhi4[k+1]
                              *(ci0*bz[ijk-kk1] + ci1*bz[ijk    ] + ci2*bz[ijk+kk1] + ci3*bz[ijk+kk2])
                         + cg3*(cg0*b [ijk    ] + cg1*b [ijk+kk1] + cg2*b [ijk+kk2] + cg3*b [ijk+kk3]) * dzhi4[k+2]
                              *(ci0*bz[ijk    ] + ci1*bz[ijk+kk1] + ci2*bz[ijk+kk2] + ci3*bz[ijk+kk3]) )
                         * dzi4[k];
      }
  }

  // top boundary
  bpe_visc[kend-1] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      bpe_visc[kend-1] -= visc *
                            ( cg0*(cg0*b [ijk-kk3] + cg1*b [ijk-kk2] + cg2*b [ijk-kk1] + cg3*b [ijk    ]) * dzhi4[kend-2]
                                 *(ci0*bz[ijk-kk3] + ci1*bz[ijk-kk2] + ci2*bz[ijk-kk1] + ci3*bz[ijk    ])
                            + cg1*(cg0*b [ijk-kk2] + cg1*b [ijk-kk1] + cg2*b [ijk    ] + cg3*b [ijk+kk1]) * dzhi4[kend-1]
                                 *(ci0*bz[ijk-kk2] + ci1*bz[ijk-kk1] + ci2*bz[ijk    ] + ci3*bz[ijk+kk1])
                            + cg2*(cg0*b [ijk-kk1] + cg1*b [ijk    ] + cg2*b [ijk+kk1] + cg3*b [ijk+kk2]) * dzhi4[kend  ]
                                 *(ci0*bz[ijk-kk1] + ci1*bz[ijk    ] + ci2*bz[ijk+kk1] + ci3*bz[ijk+kk2])
                            + cg3*(tg0*b [ijk-kk1] + tg1*b [ijk    ] + tg2*b [ijk+kk1] + tg3*b [ijk+kk2]) * dzhi4[kend+1]
                                 *(ti0*bz[ijk-kk1] + ti1*bz[ijk    ] + ti2*bz[ijk+kk1] + ti3*bz[ijk+kk2]) )
                            * dzi4[kend-1];
    }

  // calculate the dissipation term
  double dxi,dyi;
  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  bpe_diss[kstart] = 0.;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk  = i + j*jj1 + kstart*kk1;
      dzstardb = calcdzstardb(b[ijk], bsort, z);
      bpe_diss[kstart] += visc * dzstardb * (
                         std::pow( ( cg0*(ci0*b[ijk-ii3] + ci1*b[ijk-ii2] + ci2*b[ijk-ii1] + ci3*b[ijk    ])
                                   + cg1*(ci0*b[ijk-ii2] + ci1*b[ijk-ii1] + ci2*b[ijk    ] + ci3*b[ijk+ii1])
                                   + cg2*(ci0*b[ijk-ii1] + ci1*b[ijk    ] + ci2*b[ijk+ii1] + ci3*b[ijk+ii2])
                                   + cg3*(ci0*b[ijk    ] + ci1*b[ijk+ii1] + ci2*b[ijk+ii2] + ci3*b[ijk+ii3]) ) * cgi*dxi, 2)

                       + std::pow( ( cg0*(ci0*b[ijk-jj3] + ci1*b[ijk-jj2] + ci2*b[ijk-jj1] + ci3*b[ijk    ])
                                   + cg1*(ci0*b[ijk-jj2] + ci1*b[ijk-jj1] + ci2*b[ijk    ] + ci3*b[ijk+jj1])
                                   + cg2*(ci0*b[ijk-jj1] + ci1*b[ijk    ] + ci2*b[ijk+jj1] + ci3*b[ijk+jj2])
                                   + cg3*(ci0*b[ijk    ] + ci1*b[ijk+jj1] + ci2*b[ijk+jj2] + ci3*b[ijk+jj3]) ) * cgi*dyi, 2)

                       + std::pow( ( cg0*(bi0*b[ijk-kk2] + bi1*b[ijk-kk1] + bi2*b[ijk    ] + bi3*b[ijk+kk1])
                                   + cg1*(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1])
                                   + cg2*(ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2])
                                   + cg3*(ci0*b[ijk    ] + ci1*b[ijk+kk1] + ci2*b[ijk+kk2] + ci3*b[ijk+kk3]) ) * dzi4[kstart], 2) );
    }

  // interior
  for(int k=grid->kstart+1; k<grid->kend-1; k++)
  {
    bpe_diss[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj1 + k*kk1;
        dzstardb = calcdzstardb(b[ijk], bsort, z);
        bpe_diss[k] += visc * dzstardb * (
                        std::pow( ( cg0*(ci0*b[ijk-ii3] + ci1*b[ijk-ii2] + ci2*b[ijk-ii1] + ci3*b[ijk    ])
                                  + cg1*(ci0*b[ijk-ii2] + ci1*b[ijk-ii1] + ci2*b[ijk    ] + ci3*b[ijk+ii1])
                                  + cg2*(ci0*b[ijk-ii1] + ci1*b[ijk    ] + ci2*b[ijk+ii1] + ci3*b[ijk+ii2])
                                  + cg3*(ci0*b[ijk    ] + ci1*b[ijk+ii1] + ci2*b[ijk+ii2] + ci3*b[ijk+ii3]) ) * cgi*dxi, 2)

                      + std::pow( ( cg0*(ci0*b[ijk-jj3] + ci1*b[ijk-jj2] + ci2*b[ijk-jj1] + ci3*b[ijk    ])
                                  + cg1*(ci0*b[ijk-jj2] + ci1*b[ijk-jj1] + ci2*b[ijk    ] + ci3*b[ijk+jj1])
                                  + cg2*(ci0*b[ijk-jj1] + ci1*b[ijk    ] + ci2*b[ijk+jj1] + ci3*b[ijk+jj2])
                                  + cg3*(ci0*b[ijk    ] + ci1*b[ijk+jj1] + ci2*b[ijk+jj2] + ci3*b[ijk+jj3]) ) * cgi*dyi, 2)

                      + std::pow( ( cg0*(ci0*b[ijk-kk3] + ci1*b[ijk-kk2] + ci2*b[ijk-kk1] + ci3*b[ijk    ])
                                  + cg1*(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1])
                                  + cg2*(ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2])
                                  + cg3*(ci0*b[ijk    ] + ci1*b[ijk+kk1] + ci2*b[ijk+kk2] + ci3*b[ijk+kk3]) ) * dzi4[k], 2) );
      }
  }

  // top
  bpe_diss[kend-1] = 0.;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      dzstardb = calcdzstardb(b[ijk], bsort, z);
      bpe_diss[kend-1] += visc * dzstardb * (
                      std::pow( ( cg0*(ci0*b[ijk-ii3] + ci1*b[ijk-ii2] + ci2*b[ijk-ii1] + ci3*b[ijk    ])
                                + cg1*(ci0*b[ijk-ii2] + ci1*b[ijk-ii1] + ci2*b[ijk    ] + ci3*b[ijk+ii1])
                                + cg2*(ci0*b[ijk-ii1] + ci1*b[ijk    ] + ci2*b[ijk+ii1] + ci3*b[ijk+ii2])
                                + cg3*(ci0*b[ijk    ] + ci1*b[ijk+ii1] + ci2*b[ijk+ii2] + ci3*b[ijk+ii3]) ) * cgi*dxi, 2)

                    + std::pow( ( cg0*(ci0*b[ijk-jj3] + ci1*b[ijk-jj2] + ci2*b[ijk-jj1] + ci3*b[ijk    ])
                                + cg1*(ci0*b[ijk-jj2] + ci1*b[ijk-jj1] + ci2*b[ijk    ] + ci3*b[ijk+jj1])
                                + cg2*(ci0*b[ijk-jj1] + ci1*b[ijk    ] + ci2*b[ijk+jj1] + ci3*b[ijk+jj2])
                                + cg3*(ci0*b[ijk    ] + ci1*b[ijk+jj1] + ci2*b[ijk+jj2] + ci3*b[ijk+jj3]) ) * cgi*dyi, 2)

                    + std::pow( ( cg0*(ci0*b[ijk-kk3] + ci1*b[ijk-kk2] + ci2*b[ijk-kk1] + ci3*b[ijk    ])
                                + cg1*(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1])
                                + cg2*(ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2])
                                + cg3*(ti0*b[ijk-kk1] + ti1*b[ijk    ] + ti2*b[ijk+kk1] + ti3*b[ijk+kk2]) ) * dzi4[kend-1], 2) );
    }

  /*
  // CONVERT THE BZ FIELD INTO BACKGROUND POTENTIAL ENERGY
  // first, calculate the potential energy at the bottom, the bot field contains the zsort at the bottom boundary
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj1;
      ijk = i + j*jj1 + kstart*kk1;
      bzbot[ij] *= -(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk] + ci3*b[ijk+kk1]);
    }


  // calculate the potential energy at the top, the top field contains the zsort at the top boundary
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj1;
      ijk = i + j*jj1 + (kend-1)*kk1;
      bztop[ij] *= -(ci0*b[ijk-kk1] + ci1*b[ijk] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2]);
    }

  // calculate the potential energy
  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj1 + k*kk1;
        bz[ijk] = -b[ijk] * bz[ijk];
      }

  // calculate the ghost cells at the bottom
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj1;
      ijk = i + j*jj1 + kstart*kk1;
      bz[ijk-kk1] = (8./3.)*bzbot[ij] - 2.*bz[ijk] + (1./3.)*bz[ijk+kk1];
      bz[ijk-kk2] = 8.*bzbot[ij] - 9.*bz[ijk] + 2.*bz[ijk+kk1];
    }

  // calculate the ghost cells at the top
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj1;
      ijk = i + j*jj1 + (kend-1)*kk1;
      bz[ijk+kk1] = (8./3.)*bztop[ij] - 2.*bz[ijk] + (1./3.)*bz[ijk-kk1];
      bz[ijk+kk2] = 8.*bztop[ij] - 9.*bz[ijk] + 2.*bz[ijk-kk1];
    }
  */

  // calculate the advective transport term
  // bottom boundary
  bpe_turb[kstart] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk = i + j*jj1 + kstart*kk1;
      bpe_turb[kstart] += bz[ijk]*
                          ( cg0*(w[ijk-kk1] * (bi0*b[ijk-kk2] + bi1*b[ijk-kk1] + bi2*b[ijk    ] + bi3*b[ijk+kk1]))
                          + cg1*(w[ijk    ] * (ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1]))
                          + cg2*(w[ijk+kk1] * (ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2]))
                          + cg3*(w[ijk+kk2] * (ci0*b[ijk    ] + ci1*b[ijk+kk1] + ci2*b[ijk+kk2] + ci3*b[ijk+kk3])) )
                          * dzi4[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; ++k)
  {
    bpe_turb[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj1 + k*kk1;
        bpe_turb[k] += bz[ijk]*
                       ( cg0*(w[ijk-kk1] * (ci0*b[ijk-kk3] + ci1*b[ijk-kk2] + ci2*b[ijk-kk1] + ci3*b[ijk    ]))
                       + cg1*(w[ijk    ] * (ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1]))
                       + cg2*(w[ijk+kk1] * (ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2]))
                       + cg3*(w[ijk+kk2] * (ci0*b[ijk    ] + ci1*b[ijk+kk1] + ci2*b[ijk+kk2] + ci3*b[ijk+kk3])) )
                       * dzi4[k];
      }
  }

  // top boundary
  bpe_turb[kend-1] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      bpe_turb[kend-1] += bz[ijk]*
                          ( cg0*(w[ijk-kk1] * (ci0*b[ijk-kk3] + ci1*b[ijk-kk2] + ci2*b[ijk-kk1] + ci3*b[ijk    ]))
                          + cg1*(w[ijk    ] * (ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1]))
                          + cg2*(w[ijk+kk1] * (ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2]))
                          + cg3*(w[ijk+kk2] * (ti0*b[ijk-kk1] + ti1*b[ijk    ] + ti2*b[ijk+kk1] + ti3*b[ijk+kk2])) )
                          * dzi4[kend-1];
    }

  master->sum(bpe_turb, grid->kcells);
  master->sum(bpe_visc, grid->kcells);
  master->sum(bpe_diss, grid->kcells);

  int n = grid->itot*grid->jtot;
  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    bpe_turb[k] /= n;
    bpe_visc[k] /= n;
    bpe_diss[k] /= n;
  }

  return 0;
}
