/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
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
#include "model.h"
#include "thermo.h"
#include "stats.h"
#include <netcdfcpp.h>

#define NO_OFFSET 0.

cbudget::cbudget(cmodel *modelin)
{
  model = modelin;

  umodel = NULL;
  vmodel = NULL;
}

cbudget::~cbudget()
{
  delete[] umodel;
  delete[] vmodel;
}

int cbudget::readinifile(cinput *inputin)
{
  int nerror = 0;
  nerror += inputin->getItem(&swbudget, "budget", "swbudget", "", "0");

  return nerror;
}

int cbudget::init()
{
  // copy the pointers
  grid   = model->grid;
  fields = model->fields;
  stats  = model->stats;

  // if the stats is disabled, also disable the budget stats
  if(stats->getsw() == "0")
    swbudget = "0";

  umodel = new double[grid->kcells];
  vmodel = new double[grid->kcells];

  return 0;
}

int cbudget::create()
{
  if(swbudget == "0")
    return 0;

  // add the profiles to the statistics
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

  if(model->thermo->getsw() != "0")
  {
    stats->addprof("w2_buoy" , "Buoyancy production/destruction term in W2 budget" , "m2 s-3", "zh");
    stats->addprof("tke_buoy", "Buoyancy production/destruction term in TKE budget", "m2 s-3", "z" );
  }

  return 0;
}

int cbudget::execstats()
{
  if(swbudget == "0")
    return 0;

  // calculate the mean of the fields
  stats->calcmean(fields->u->data, umodel, NO_OFFSET);
  stats->calcmean(fields->v->data, vmodel, NO_OFFSET);

  if(grid->swspatialorder == "4")
  {
    // calculate the TKE budget
    calctkebudget(fields->u->data, fields->v->data, fields->w->data, fields->s["p"]->data,
                  fields->s["tmp1"]->data, fields->s["tmp2"]->data,
                  umodel, vmodel,
                  stats->profs["u2_shear"].data, stats->profs["v2_shear"].data, stats->profs["tke_shear"].data,
                  stats->profs["u2_turb"].data, stats->profs["v2_turb"].data, stats->profs["w2_turb"].data, stats->profs["tke_turb"].data,
                  stats->profs["u2_visc"].data, stats->profs["v2_visc"].data, stats->profs["w2_visc"].data, stats->profs["tke_visc"].data,
                  stats->profs["u2_diss"].data, stats->profs["v2_diss"].data, stats->profs["w2_diss"].data, stats->profs["tke_diss"].data,
                  stats->profs["w2_pres"].data, stats->profs["tke_pres"].data,
                  stats->profs["u2_rdstr"].data, stats->profs["v2_rdstr"].data, stats->profs["w2_rdstr"].data,
                  grid->dzi4, grid->dzhi4, fields->visc);

    // calculate the buoyancy term of the TKE budget
    if(model->thermo->getsw() != "0")
    {
      model->thermo->getthermofield(fields->sd["tmp1"], fields->sd["tmp2"], "b");
      calctkebudget_buoy(fields->w->data, fields->s["tmp1"]->data,
                    stats->profs["w2_buoy"].data, stats->profs["tke_buoy"].data);
    }
  }

  return 0;
}

int cbudget::calctkebudget(double * restrict u, double * restrict v, double * restrict w, double * restrict p,
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
  // get w on the x and y location
  grid->interpolatex_4th(wx, w, 0);
  grid->interpolatey_4th(wy, w, 0);

  int ijk,ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;
  kk3 = 3*grid->ijcells;

  double n = grid->imax*grid->jmax;

  // calculate the shear term u'w*dumean/dz
  for(int k=grid->kstart; k<grid->kend; ++k)
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

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    u2_shear [k] /= n;
    v2_shear [k] /= n;
    tke_shear[k] /= n;
  }

  grid->getprof(u2_shear , grid->kcells);
  grid->getprof(v2_shear , grid->kcells);
  grid->getprof(tke_shear, grid->kcells);

  // calculate the turbulent transport term
  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    u2_turb [k] = 0.;
    v2_turb [k] = 0.;
    w2_turb [k] = 0.;
    tke_turb[k] = 0.;

    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        u2_turb[k]  -= ( cg0*((ci0*std::pow(u[ijk-kk3]-umean[k-3],2.) + ci1*std::pow(u[ijk-kk2]-umean[k-2],2.) + ci2*std::pow(u[ijk-kk1]-umean[k-1],2.) + ci3*std::pow(u[ijk    ]-umean[k  ],2.))*wx[ijk-kk1])
                       + cg1*((ci0*std::pow(u[ijk-kk2]-umean[k-2],2.) + ci1*std::pow(u[ijk-kk1]-umean[k-1],2.) + ci2*std::pow(u[ijk    ]-umean[k  ],2.) + ci3*std::pow(u[ijk+kk1]-umean[k+1],2.))*wx[ijk    ])
                       + cg2*((ci0*std::pow(u[ijk-kk1]-umean[k-1],2.) + ci1*std::pow(u[ijk    ]-umean[k  ],2.) + ci2*std::pow(u[ijk+kk1]-umean[k+1],2.) + ci3*std::pow(u[ijk+kk2]-umean[k+2],2.))*wx[ijk+kk1])
                       + cg3*((ci0*std::pow(u[ijk    ]-umean[k  ],2.) + ci1*std::pow(u[ijk+kk1]-umean[k+1],2.) + ci2*std::pow(u[ijk+kk2]-umean[k+2],2.) + ci3*std::pow(u[ijk+kk3]-umean[k+3],2.))*wx[ijk+kk2]) ) * dzi4[k];

        v2_turb[k]  -= ( cg0*((ci0*std::pow(v[ijk-kk3]-vmean[k-3],2.) + ci1*std::pow(v[ijk-kk2]-vmean[k-2],2.) + ci2*std::pow(v[ijk-kk1]-vmean[k-1],2.) + ci3*std::pow(v[ijk    ]-vmean[k  ],2.))*wy[ijk-kk1]) 
                       + cg1*((ci0*std::pow(v[ijk-kk2]-vmean[k-2],2.) + ci1*std::pow(v[ijk-kk1]-vmean[k-1],2.) + ci2*std::pow(v[ijk    ]-vmean[k  ],2.) + ci3*std::pow(v[ijk+kk1]-vmean[k+1],2.))*wy[ijk    ]) 
                       + cg2*((ci0*std::pow(v[ijk-kk1]-vmean[k-1],2.) + ci1*std::pow(v[ijk    ]-vmean[k  ],2.) + ci2*std::pow(v[ijk+kk1]-vmean[k+1],2.) + ci3*std::pow(v[ijk+kk2]-vmean[k+2],2.))*wy[ijk+kk1]) 
                       + cg3*((ci0*std::pow(v[ijk    ]-vmean[k  ],2.) + ci1*std::pow(v[ijk+kk1]-vmean[k+1],2.) + ci2*std::pow(v[ijk+kk2]-vmean[k+2],2.) + ci3*std::pow(v[ijk+kk3]-vmean[k+3],2.))*wy[ijk+kk2]) ) * dzi4[k];

        tke_turb[k] -= 0.5*( cg0*std::pow(w[ijk-kk1], 3.) + cg1*std::pow(w[ijk], 3.) + cg2*std::pow(w[ijk+kk1], 3.) + cg3*std::pow(w[ijk+kk2], 3.)) * dzi4[k];
      }
    tke_turb[k] += 0.5*(u2_turb[k] + v2_turb[k]);
  }
  for(int k=grid->kstart+1; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        w2_turb[k] -= ( cg0*(ci0*std::pow(w[ijk-kk3],3.) + ci1*std::pow(w[ijk-kk2],3.) + ci2*std::pow(w[ijk-kk1],3.) + ci3*std::pow(w[ijk    ],3.))
                      + cg1*(ci0*std::pow(w[ijk-kk2],3.) + ci1*std::pow(w[ijk-kk1],3.) + ci2*std::pow(w[ijk    ],3.) + ci3*std::pow(w[ijk+kk1],3.))
                      + cg2*(ci0*std::pow(w[ijk-kk1],3.) + ci1*std::pow(w[ijk    ],3.) + ci2*std::pow(w[ijk+kk1],3.) + ci3*std::pow(w[ijk+kk2],3.))
                      + cg3*(ci0*std::pow(w[ijk    ],3.) + ci1*std::pow(w[ijk+kk1],3.) + ci2*std::pow(w[ijk+kk2],3.) + ci3*std::pow(w[ijk+kk3],3.)) ) * dzhi4[k];
      }

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    u2_turb [k] /= n;
    v2_turb [k] /= n;
    w2_turb [k] /= n;
    tke_turb[k] /= n;
  }

  grid->getprof(u2_turb , grid->kcells);
  grid->getprof(v2_turb , grid->kcells);
  grid->getprof(w2_turb , grid->kcells);
  grid->getprof(tke_turb, grid->kcells);

  // calculate the pressure transport term
  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    w2_pres [k] = 0.;
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
  for(int k=grid->kstart+1; k<grid->kend; ++k)
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

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    w2_pres [k] /= n;
    tke_pres[k] /= n;
  }

  grid->getprof(w2_pres , grid->kcells);
  grid->getprof(tke_pres, grid->kcells);

  // calculate the viscous transport term
  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    u2_visc [k] = 0.;
    v2_visc [k] = 0.;
    w2_visc [k] = 0.;
    tke_visc[k] = 0.;

    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        u2_visc[k]  += visc * ( cg0*((cg0*std::pow(u[ijk-kk3]-umean[k-3],2.) + cg1*std::pow(u[ijk-kk2]-umean[k-2],2.) + cg2*std::pow(u[ijk-kk1]-umean[k-1],2.) + cg3*std::pow(u[ijk    ]-umean[k  ],2.)) * dzhi4[k-1])
                              + cg1*((cg0*std::pow(u[ijk-kk2]-umean[k-2],2.) + cg1*std::pow(u[ijk-kk1]-umean[k-1],2.) + cg2*std::pow(u[ijk    ]-umean[k  ],2.) + cg3*std::pow(u[ijk+kk1]-umean[k+1],2.)) * dzhi4[k  ])
                              + cg2*((cg0*std::pow(u[ijk-kk1]-umean[k-1],2.) + cg1*std::pow(u[ijk    ]-umean[k  ],2.) + cg2*std::pow(u[ijk+kk1]-umean[k+1],2.) + cg3*std::pow(u[ijk+kk2]-umean[k+2],2.)) * dzhi4[k+1])
                              + cg3*((cg0*std::pow(u[ijk    ]-umean[k  ],2.) + cg1*std::pow(u[ijk+kk1]-umean[k+1],2.) + cg2*std::pow(u[ijk+kk2]-umean[k+2],2.) + cg3*std::pow(u[ijk+kk3]-umean[k+3],2.)) * dzhi4[k+2]) ) * dzi4[k];

        v2_visc[k]  += visc * ( cg0*((cg0*std::pow(v[ijk-kk3]-vmean[k-3],2.) + cg1*std::pow(v[ijk-kk2]-vmean[k-2],2.) + cg2*std::pow(v[ijk-kk1]-vmean[k-1],2.) + cg3*std::pow(v[ijk    ]-vmean[k  ],2.)) * dzhi4[k-1])
                              + cg1*((cg0*std::pow(v[ijk-kk2]-vmean[k-2],2.) + cg1*std::pow(v[ijk-kk1]-vmean[k-1],2.) + cg2*std::pow(v[ijk    ]-vmean[k  ],2.) + cg3*std::pow(v[ijk+kk1]-vmean[k+1],2.)) * dzhi4[k  ])
                              + cg2*((cg0*std::pow(v[ijk-kk1]-vmean[k-1],2.) + cg1*std::pow(v[ijk    ]-vmean[k  ],2.) + cg2*std::pow(v[ijk+kk1]-vmean[k+1],2.) + cg3*std::pow(v[ijk+kk2]-vmean[k+2],2.)) * dzhi4[k+1])
                              + cg3*((cg0*std::pow(v[ijk    ]-vmean[k  ],2.) + cg1*std::pow(v[ijk+kk1]-vmean[k+1],2.) + cg2*std::pow(v[ijk+kk2]-vmean[k+2],2.) + cg3*std::pow(v[ijk+kk3]-vmean[k+3],2.)) * dzhi4[k+2]) ) * dzi4[k];

        tke_visc[k] += 0.5*visc * ( cg0*std::pow(w[ijk-kk1],2.) + cg1*std::pow(w[ijk],2.) + cg2*std::pow(w[ijk+kk1],2.) + cg3*std::pow(w[ijk+kk2],2.)) * dzi4[k];
      }
    tke_visc[k] += 0.5*(u2_visc[k] + v2_visc[k]);
  }
  for(int k=grid->kstart+1; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        w2_visc[k] += visc * ( cg0*((cg0*std::pow(w[ijk-kk3],2.) + cg1*std::pow(w[ijk-kk2],2.) + cg2*std::pow(w[ijk-kk1],2.) + cg3*std::pow(w[ijk    ],2.)) * dzi4[k-2])
                             + cg1*((cg0*std::pow(w[ijk-kk2],2.) + cg1*std::pow(w[ijk-kk1],2.) + cg2*std::pow(w[ijk    ],2.) + cg3*std::pow(w[ijk+kk1],2.)) * dzi4[k-1])
                             + cg2*((cg0*std::pow(w[ijk-kk1],2.) + cg1*std::pow(w[ijk    ],2.) + cg2*std::pow(w[ijk+kk1],2.) + cg3*std::pow(w[ijk+kk2],2.)) * dzi4[k  ])
                             + cg3*((cg0*std::pow(w[ijk    ],2.) + cg1*std::pow(w[ijk+kk1],2.) + cg2*std::pow(w[ijk+kk2],2.) + cg3*std::pow(w[ijk+kk3],2.)) * dzi4[k+1]) ) * dzhi4[k];
      }

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    u2_visc [k] /= n;
    v2_visc [k] /= n;
    w2_visc [k] /= n;
    tke_visc[k] /= n;
  }

  grid->getprof(u2_visc , grid->kcells);
  grid->getprof(v2_visc , grid->kcells);
  grid->getprof(w2_visc , grid->kcells);
  grid->getprof(tke_visc, grid->kcells);


  // calculate the dissipation
  double dxi,dyi;
  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    u2_diss [k] = 0.;
    v2_diss [k] = 0.;
    w2_diss [k] = 0.;
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
                                   + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+ii1]-umean[k]) + ci2*(u[ijk+ii2]-umean[k]) + ci3*(u[ijk+ii3]-umean[k]))) ) * cgi*dxi, 2.)

                       + std::pow( ( cg0*((ci0*(u[ijk-jj3]-umean[k]) + ci1*(u[ijk-jj2]-umean[k]) + ci2*(u[ijk-jj1]-umean[k]) + ci3*(u[ijk    ]-umean[k])))
                                   + cg1*((ci0*(u[ijk-jj2]-umean[k]) + ci1*(u[ijk-jj1]-umean[k]) + ci2*(u[ijk    ]-umean[k]) + ci3*(u[ijk+jj1]-umean[k])))
                                   + cg2*((ci0*(u[ijk-jj1]-umean[k]) + ci1*(u[ijk    ]-umean[k]) + ci2*(u[ijk+jj1]-umean[k]) + ci3*(u[ijk+jj2]-umean[k])))
                                   + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+jj1]-umean[k]) + ci2*(u[ijk+jj2]-umean[k]) + ci3*(u[ijk+jj3]-umean[k]))) ) * cgi*dyi, 2.)

                       + std::pow( ( cg0*((ci0*(u[ijk-kk3]-umean[k-3]) + ci1*(u[ijk-kk2]-umean[k-2]) + ci2*(u[ijk-kk1]-umean[k-1]) + ci3*(u[ijk    ]-umean[k  ])))
                                   + cg1*((ci0*(u[ijk-kk2]-umean[k-2]) + ci1*(u[ijk-kk1]-umean[k-1]) + ci2*(u[ijk    ]-umean[k  ]) + ci3*(u[ijk+kk1]-umean[k+1])))
                                   + cg2*((ci0*(u[ijk-kk1]-umean[k-1]) + ci1*(u[ijk    ]-umean[k  ]) + ci2*(u[ijk+kk1]-umean[k+1]) + ci3*(u[ijk+kk2]-umean[k+2])))
                                   + cg3*((ci0*(u[ijk    ]-umean[k  ]) + ci1*(u[ijk+kk1]-umean[k+1]) + ci2*(u[ijk+kk2]-umean[k+2]) + ci3*(u[ijk+kk3]-umean[k+3]))) ) * dzi4[k], 2.) );

        v2_diss[k]  -= 2.*visc * (
                         std::pow( ( cg0*((ci0*(v[ijk-ii3]-vmean[k]) + ci1*(v[ijk-ii2]-vmean[k]) + ci2*(v[ijk-ii1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                                   + cg1*((ci0*(v[ijk-ii2]-vmean[k]) + ci1*(v[ijk-ii1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+ii1]-vmean[k])))
                                   + cg2*((ci0*(v[ijk-ii1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+ii1]-vmean[k]) + ci3*(v[ijk+ii2]-vmean[k])))
                                   + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+ii1]-vmean[k]) + ci2*(v[ijk+ii2]-vmean[k]) + ci3*(v[ijk+ii3]-vmean[k]))) ) * cgi*dxi, 2.)

                       + std::pow( ( cg0*((ci0*(v[ijk-jj3]-vmean[k]) + ci1*(v[ijk-jj2]-vmean[k]) + ci2*(v[ijk-jj1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                                   + cg1*((ci0*(v[ijk-jj2]-vmean[k]) + ci1*(v[ijk-jj1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+jj1]-vmean[k])))
                                   + cg2*((ci0*(v[ijk-jj1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+jj1]-vmean[k]) + ci3*(v[ijk+jj2]-vmean[k])))
                                   + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+jj1]-vmean[k]) + ci2*(v[ijk+jj2]-vmean[k]) + ci3*(v[ijk+jj3]-vmean[k]))) ) * cgi*dyi, 2.)

                       + std::pow( ( cg0*((ci0*(v[ijk-kk3]-vmean[k-3]) + ci1*(v[ijk-kk2]-vmean[k-2]) + ci2*(v[ijk-kk1]-vmean[k-1]) + ci3*(v[ijk    ]-vmean[k  ])))
                                   + cg1*((ci0*(v[ijk-kk2]-vmean[k-2]) + ci1*(v[ijk-kk1]-vmean[k-1]) + ci2*(v[ijk    ]-vmean[k  ]) + ci3*(v[ijk+kk1]-vmean[k+1])))
                                   + cg2*((ci0*(v[ijk-kk1]-vmean[k-1]) + ci1*(v[ijk    ]-vmean[k  ]) + ci2*(v[ijk+kk1]-vmean[k+1]) + ci3*(v[ijk+kk2]-vmean[k+2])))
                                   + cg3*((ci0*(v[ijk    ]-vmean[k  ]) + ci1*(v[ijk+kk1]-vmean[k+1]) + ci2*(v[ijk+kk2]-vmean[k+2]) + ci3*(v[ijk+kk3]-vmean[k+3]))) ) * dzi4[k], 2.) );

        tke_diss[k] -= visc * (
                         std::pow( (cg0*w[ijk-ii1] + cg1*w[ijk] + cg2*w[ijk+ii1] + cg3*w[ijk+ii2]) * cgi*dxi, 2.)
                       + std::pow( (cg0*w[ijk-jj1] + cg1*w[ijk] + cg2*w[ijk+jj1] + cg3*w[ijk+jj2]) * cgi*dyi, 2.)
                       + std::pow( (cg0*w[ijk-kk1] + cg1*w[ijk] + cg2*w[ijk+kk1] + cg3*w[ijk+kk2]) * dzi4[k], 2.) );
      }
    tke_diss[k] += 0.5*(u2_diss[k] + v2_diss[k]);
  }
  for(int k=grid->kstart+1; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj1 + k*kk1;
        w2_diss[k]  -= 2.*visc * (
                         std::pow( ( cg0*(ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ])
                                   + cg1*(ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1])
                                   + cg2*(ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2])
                                   + cg3*(ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3]) ) * cgi*dxi, 2.)

                       + std::pow( ( cg0*(ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ])
                                   + cg1*(ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1])
                                   + cg2*(ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2])
                                   + cg3*(ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3]) ) * cgi*dyi, 2.)

                       + std::pow( ( cg0*(ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ])
                                   + cg1*(ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1])
                                   + cg2*(ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])
                                   + cg3*(ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3]) ) * dzhi4[k], 2.) );
      }

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    u2_diss [k] /= n;
    v2_diss [k] /= n;
    w2_diss [k] /= n;
    tke_diss[k] /= n;
  }

  grid->getprof(u2_diss , grid->kcells);
  grid->getprof(v2_diss , grid->kcells);
  grid->getprof(w2_diss , grid->kcells);
  grid->getprof(tke_diss, grid->kcells);


  // calculate the pressure redistribution term
  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    u2_rdstr [k] = 0.;
    v2_rdstr [k] = 0.;
    w2_rdstr [k] = 0.;

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

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    u2_rdstr [k] /= n;
    v2_rdstr [k] /= n;
    w2_rdstr [k] /= n;
  }

  grid->getprof(u2_rdstr , grid->kcells);
  grid->getprof(v2_rdstr , grid->kcells);
  grid->getprof(w2_rdstr , grid->kcells);

  return 0;
}

int cbudget::calctkebudget_buoy(double * restrict w, double * restrict b,
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

  grid->getprof(w2_buoy , grid->kcells);
  grid->getprof(tke_buoy, grid->kcells);

  return 0;
}

