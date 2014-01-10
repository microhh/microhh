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
#include "grid.h"
#include "fields.h"
#include "stats_dns.h"
#include "defines.h"
#include "model.h"
#include <netcdfcpp.h>

#define NO_OFFSET 0.

cstats_dns::cstats_dns(cmodel *modelin) : cstats(modelin)
{
  allocated   = false;
  initialized = false;
}

cstats_dns::~cstats_dns()
{
  if(initialized)
    delete dataFile;

  if(allocated)
  {
    delete[] umodel;
    delete[] vmodel;

    // delete the profiles
    for(profmap::const_iterator it=profs.begin(); it!=profs.end(); ++it)
      delete[] it->second.data;
  }
}

int cstats_dns::readinifile(cinput *inputin)
{
  int nerror = 0;
  nerror += inputin->getItem(&statstime, "stats", "statstime", "");
  return nerror;
}

int cstats_dns::init(double ifactor)
{
  istatstime = (unsigned long)(ifactor * statstime);

  umodel = new double[grid->kcells];
  vmodel = new double[grid->kcells];

  allocated = true;

  // set the number of stats to zero
  nstats = 0;

  return 0;
}

int cstats_dns::create(int n)
{
  int nerror = 0;

  // create a NetCDF file for the statistics
  if(mpi->mpiid == 0)
  {
    char filename[256];
    std::sprintf(filename, "%s.%07d.nc", mpi->simname.c_str(), n);
    dataFile = new NcFile(filename, NcFile::New);
    if(!dataFile->is_valid())
    {
      std::printf("ERROR cannot write statistics file\n");
      ++nerror;
    }
    else
      initialized = true;
  }

  // crash on all processes in case the file could not be written
  mpi->broadcast(&nerror, 1);
  if(nerror)
    return 1;

  if(mpi->mpiid == 0)
  {
    // create dimensions
    z_dim  = dataFile->add_dim("z" , grid->kmax);
    zh_dim = dataFile->add_dim("zh", grid->kmax+1);
    t_dim  = dataFile->add_dim("t");

    // create variables
    iter_var = dataFile->add_var("iter", ncInt   , t_dim );
    t_var    = dataFile->add_var("t"   , ncDouble, t_dim );
    z_var    = dataFile->add_var("z"   , ncDouble, z_dim );
    zh_var   = dataFile->add_var("zh"  , ncDouble, zh_dim);
  }

  // means
  addprof("u", "z" );
  addprof("v", "z" );
  addprof("w", "zh");
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    addprof(it->first, "z");
  addprof("p", "z");

  // 2nd order
  addprof("u2", "z" );
  addprof("v2", "z" );
  addprof("w2", "zh");
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    addprof(it->first+"2", "z");

  // 3rd order
  addprof("u3", "z" );
  addprof("v3", "z" );
  addprof("w3", "zh");
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    addprof(it->first+"3", "z");

  // 4th order
  addprof("u4", "z" );
  addprof("v4", "z" );
  addprof("w4", "zh");
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    addprof(it->first+"4", "z");

  // gradients
  addprof("ugrad", "zh");
  addprof("vgrad", "zh");
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    addprof(it->first+"grad", "zh");

  // turbulent fluxes
  addprof("uw", "zh");
  addprof("vw", "zh");
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    addprof(it->first+"w", "zh");

  // diffusive fluxes
  addprof("udiff", "zh");
  addprof("vdiff", "zh");
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    addprof(it->first+"diff", "zh");

  // total fluxes
  addprof("uflux", "zh");
  addprof("vflux", "zh");
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    addprof(it->first+"flux", "zh");

  addprof("u2_shear" , "z");
  addprof("v2_shear" , "z");
  addprof("tke_shear", "z");

  addprof("u2_turb" , "z" );
  addprof("v2_turb" , "z" );
  addprof("w2_turb" , "zh");
  addprof("tke_turb", "z" );

  addprof("u2_visc" , "z" );
  addprof("v2_visc" , "z" );
  addprof("w2_visc" , "zh");
  addprof("tke_visc", "z" );

  addprof("u2_diss" , "z" );
  addprof("v2_diss" , "z" );
  addprof("w2_diss" , "zh");
  addprof("tke_diss", "z" );

  addprof("w2_pres" , "zh");
  addprof("tke_pres", "z" );

  addprof("u2_rdstr", "z" );
  addprof("v2_rdstr", "z" );
  addprof("w2_rdstr", "zh");

  addprof("w2_buoy" , "zh");
  addprof("tke_buoy", "z" );

  if(mpi->mpiid == 0)
  {
    // save the grid variables
    z_var ->put(&grid->z [grid->kstart], grid->kmax  );
    zh_var->put(&grid->zh[grid->kstart], grid->kmax+1);

    dataFile->sync();
  }

  return 0;
}

unsigned long cstats_dns::gettimelim(unsigned long itime)
{
  unsigned long idtlim = istatstime -  itime % istatstime;
  return idtlim;
}

int cstats_dns::exec(int iteration, double time, unsigned long itime)
{
  // check if time for execution
  if(itime % istatstime != 0)
    return 0;

  if(mpi->mpiid == 0) std::printf("Saving stats for time %f\n", time);

  if(mpi->mpiid == 0)
  {
    t_var   ->put_rec(&time     , nstats);
    iter_var->put_rec(&iteration, nstats);
  }

  // PROFILES
  // calculate means
  calcmean(fields->u->data, profs["u"].data, grid->utrans);
  calcmean(fields->v->data, profs["v"].data, grid->vtrans);
  calcmean(fields->w->data, profs["w"].data, NO_OFFSET);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    calcmean(it->second->data, profs[it->first].data, NO_OFFSET);

  calcmean(fields->s["p"]->data, profs["p"].data, NO_OFFSET);

  // calculate model means without correction for transformation
  calcmean(fields->u->data, umodel, NO_OFFSET);
  calcmean(fields->v->data, vmodel, NO_OFFSET);

  // 2nd order
  calcmoment(fields->u->data, umodel, profs["u2"].data, 2., 0);
  calcmoment(fields->v->data, vmodel, profs["v2"].data, 2., 0);
  calcmoment(fields->w->data, profs["w"].data, profs["w2"].data, 2., 1);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    calcmoment(it->second->data, profs[it->first].data, profs[it->first+"2"].data, 2., 0);

  // 3rd order
  calcmoment(fields->u->data, umodel, profs["u3"].data, 3., 0);
  calcmoment(fields->v->data, vmodel, profs["v3"].data, 3., 0);
  calcmoment(fields->w->data, profs["w"].data, profs["w3"].data, 3., 1);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    calcmoment(it->second->data, profs[it->first].data, profs[it->first+"3"].data, 3., 0);

  // 4th order
  calcmoment(fields->u->data, umodel, profs["u4"].data, 4., 0);
  calcmoment(fields->v->data, vmodel, profs["v4"].data, 4., 0);
  calcmoment(fields->w->data, profs["w"].data, profs["w4"].data, 4., 1);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    calcmoment(it->second->data, profs[it->first].data, profs[it->first+"4"].data, 3., 0);

  calcgrad(fields->u->data, profs["ugrad"].data, grid->dzhi4);
  calcgrad(fields->v->data, profs["vgrad"].data, grid->dzhi4);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    calcgrad(it->second->data, profs[it->first+"grad"].data, grid->dzhi4);

  // calculate turbulent fluxes
  calcflux(fields->u->data, fields->w->data, profs["uw"].data, fields->s["tmp1"]->data, 1, 0);
  calcflux(fields->v->data, fields->w->data, profs["vw"].data, fields->s["tmp1"]->data, 0, 1);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    calcflux(it->second->data, fields->w->data, profs[it->first+"w"].data, fields->s["tmp1"]->data, 0, 0);

  // calculate diffusive fluxes
  calcdiff(fields->u->data, profs["udiff"].data, grid->dzhi4, fields->visc);
  calcdiff(fields->v->data, profs["vdiff"].data, grid->dzhi4, fields->visc);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    calcdiff(it->second->data, profs[it->first+"diff"].data, grid->dzhi4, it->second->visc);

  addfluxes(profs["uflux"].data, profs["uw"].data, profs["udiff"].data);
  addfluxes(profs["vflux"].data, profs["vw"].data, profs["vdiff"].data);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    addfluxes(profs[it->first+"flux"].data, profs[it->first+"w"].data, profs[it->first+"diff"].data);

  // calculate the TKE budget
  calctkebudget(fields->u->data, fields->v->data, fields->w->data, fields->s["p"]->data, fields->s["s"]->data,
                fields->s["tmp1"]->data, fields->s["tmp2"]->data,
                umodel, vmodel,
                profs["u2_shear"].data, profs["v2_shear"].data, profs["tke_shear"].data,
                profs["u2_turb"].data, profs["v2_turb"].data, profs["w2_turb"].data, profs["tke_turb"].data,
                profs["u2_visc"].data, profs["v2_visc"].data, profs["w2_visc"].data, profs["tke_visc"].data,
                profs["u2_diss"].data, profs["v2_diss"].data, profs["w2_diss"].data, profs["tke_diss"].data,
                profs["w2_pres"].data, profs["tke_pres"].data,
                profs["u2_rdstr"].data, profs["v2_rdstr"].data, profs["w2_rdstr"].data,
                profs["w2_buoy"].data, profs["tke_buoy"].data,
                grid->dzi4, grid->dzhi4, fields->visc);

  // put the data into the NetCDF file
  if(mpi->mpiid == 0)
  {
    for(profmap::const_iterator it=profs.begin(); it!=profs.end(); ++it)
      profs[it->first].ncvar->put_rec(&profs[it->first].data[grid->kstart], nstats);

    dataFile->sync();
  }

  ++nstats;

  return 0;
}

int cstats_dns::addprof(std::string name, std::string zloc)
{
  // create the NetCDF variable
  if(mpi->mpiid == 0)
  {
    if(zloc == "z")
      profs[name].ncvar = dataFile->add_var(name.c_str(), ncDouble, t_dim, z_dim );
    else if(zloc == "zh")
      profs[name].ncvar = dataFile->add_var(name.c_str(), ncDouble, t_dim, zh_dim);
  }

  // and allocate the memory
  profs[name].data = new double[grid->kcells];

  return 0;
}

// COMPUTATIONAL KERNELS
int cstats_dns::calcmean(double * restrict data, double * restrict prof, double offset)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  
  for(int k=0; k<grid->kcells; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += data[ijk] + offset;
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=0; k<grid->kcells; k++)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats_dns::calcmoment(double * restrict data, double * restrict datamean, double * restrict prof, double power, int a)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  
  for(int k=grid->kstart; k<grid->kend+a; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += std::pow(data[ijk]-datamean[k], power);
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+a; k++)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats_dns::calcflux(double * restrict data, double * restrict w, double * restrict prof, double * restrict tmp1, int locx, int locy)
{
  int ijk,jj,kk1,kk2;

  jj  = 1*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  // set a pointer to the field that contains w, either interpolated or the original
  double * restrict calcw = w;
  if(locx == 1)
  {
    grid->interpolatex_4th(tmp1, w, 0);
    calcw = tmp1;
  }
  else if(locy == 1)
  {
    grid->interpolatey_4th(tmp1, w, 0);
    calcw = tmp1;
  }
  
  for(int k=grid->kstart; k<grid->kend+1; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk1;
        prof[k] += (ci0*data[ijk-kk2] + ci1*data[ijk-kk1] + ci2*data[ijk] + ci3*data[ijk+kk1])*calcw[ijk];
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+1; k++)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats_dns::calcgrad(double * restrict data, double * restrict prof, double * restrict dzhi4)
{
  int ijk,jj,kk1,kk2;

  jj  = 1*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  
  for(int k=grid->kstart; k<grid->kend+1; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk1;
        prof[k] += (cg0*data[ijk-kk2] + cg1*data[ijk-kk1] + cg2*data[ijk] + cg3*data[ijk+kk1])*dzhi4[k];
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+1; k++)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats_dns::calcdiff(double * restrict data, double * restrict prof, double * restrict dzhi4, double visc)
{
  int ijk,jj,kk1,kk2;

  jj  = 1*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  
  for(int k=grid->kstart; k<grid->kend+1; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk1;
        prof[k] += -visc*(cg0*data[ijk-kk2] + cg1*data[ijk-kk1] + cg2*data[ijk] + cg3*data[ijk+kk1])*dzhi4[k];
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+1; k++)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats_dns::calctkebudget(double * restrict u, double * restrict v, double * restrict w, double * restrict p, double * restrict b,
                              double * restrict wx, double * restrict wy,
                              double * restrict umean, double * restrict vmean,
                              double * restrict u2_shear, double * restrict v2_shear, double * restrict tke_shear,
                              double * restrict u2_turb, double * restrict v2_turb, double * restrict w2_turb, double * restrict tke_turb,
                              double * restrict u2_visc, double * restrict v2_visc, double * restrict w2_visc, double * restrict tke_visc,
                              double * restrict u2_diss, double * restrict v2_diss, double * restrict w2_diss, double * restrict tke_diss,
                              double * restrict w2_pres, double * restrict tke_pres,
                              double * restrict u2_rdstr, double * restrict v2_rdstr, double * restrict w2_rdstr,
                              double * restrict w2_buoy, double * restrict tke_buoy,
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
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  double n = grid->imax*grid->jmax;


  // calculate the shear term u'w*dumean/dz
  for(int k=grid->kstart; k<grid->kend; k++)
  {
    u2_shear [k] = 0.;
    v2_shear [k] = 0.;
    tke_shear[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
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

  for(int k=grid->kstart; k<grid->kend; k++)
  {
    u2_shear [k] /= n;
    v2_shear [k] /= n;
    tke_shear[k] /= n;
  }

  grid->getprof(u2_shear , grid->kcells);
  grid->getprof(v2_shear , grid->kcells);
  grid->getprof(tke_shear, grid->kcells);


  // calculate the turbulent transport term
  for(int k=grid->kstart; k<grid->kend; k++)
  {
    u2_turb [k] = 0.;
    v2_turb [k] = 0.;
    w2_turb [k] = 0.;
    tke_turb[k] = 0.;

    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
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
  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj1 + k*kk1;
        w2_turb[k] -= ( cg0*(ci0*std::pow(w[ijk-kk3],3.) + ci1*std::pow(w[ijk-kk2],3.) + ci2*std::pow(w[ijk-kk1],3.) + ci3*std::pow(w[ijk    ],3.))
                      + cg1*(ci0*std::pow(w[ijk-kk2],3.) + ci1*std::pow(w[ijk-kk1],3.) + ci2*std::pow(w[ijk    ],3.) + ci3*std::pow(w[ijk+kk1],3.))
                      + cg2*(ci0*std::pow(w[ijk-kk1],3.) + ci1*std::pow(w[ijk    ],3.) + ci2*std::pow(w[ijk+kk1],3.) + ci3*std::pow(w[ijk+kk2],3.))
                      + cg3*(ci0*std::pow(w[ijk    ],3.) + ci1*std::pow(w[ijk+kk1],3.) + ci2*std::pow(w[ijk+kk2],3.) + ci3*std::pow(w[ijk+kk3],3.)) ) * dzhi4[k];
      }

  for(int k=grid->kstart; k<grid->kend; k++)
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
  for(int k=grid->kstart; k<grid->kend; k++)
  {
    w2_pres [k] = 0.;
    tke_pres[k] = 0.;

    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj1 + k*kk1;
        tke_pres[k] -= ( cg0*((ci0*p[ijk-kk3] + ci1*p[ijk-kk2] + ci2*p[ijk-kk1] + ci3*p[ijk    ])*w[ijk-kk1])
                       + cg1*((ci0*p[ijk-kk2] + ci1*p[ijk-kk1] + ci2*p[ijk    ] + ci3*p[ijk+kk1])*w[ijk    ])
                       + cg2*((ci0*p[ijk-kk1] + ci1*p[ijk    ] + ci2*p[ijk+kk1] + ci3*p[ijk+kk2])*w[ijk+kk1])
                       + cg3*((ci0*p[ijk    ] + ci1*p[ijk+kk1] + ci2*p[ijk+kk2] + ci3*p[ijk+kk3])*w[ijk+kk2]) ) * dzi4[k];
      }
  }
  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj1 + k*kk1;
        w2_pres[k] -= 2.*( cg0*((ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ])*p[ijk-kk2])
                         + cg1*((ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1])*p[ijk-kk1])
                         + cg2*((ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])*p[ijk    ])
                         + cg3*((ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3])*p[ijk+kk1]) ) * dzhi4[k];
      }

  for(int k=grid->kstart; k<grid->kend; k++)
  {
    w2_pres [k] /= n;
    tke_pres[k] /= n;
  }

  grid->getprof(w2_pres , grid->kcells);
  grid->getprof(tke_pres, grid->kcells);


  // calculate the viscous transport term
  for(int k=grid->kstart; k<grid->kend; k++)
  {
    u2_visc [k] = 0.;
    v2_visc [k] = 0.;
    w2_visc [k] = 0.;
    tke_visc[k] = 0.;

    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
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
  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj1 + k*kk1;
        w2_visc[k] += visc * ( cg0*((cg0*std::pow(w[ijk-kk3],2.) + cg1*std::pow(w[ijk-kk2],2.) + cg2*std::pow(w[ijk-kk1],2.) + cg3*std::pow(w[ijk    ],2.)) * dzi4[k-2])
                             + cg1*((cg0*std::pow(w[ijk-kk2],2.) + cg1*std::pow(w[ijk-kk1],2.) + cg2*std::pow(w[ijk    ],2.) + cg3*std::pow(w[ijk+kk1],2.)) * dzi4[k-1])
                             + cg2*((cg0*std::pow(w[ijk-kk1],2.) + cg1*std::pow(w[ijk    ],2.) + cg2*std::pow(w[ijk+kk1],2.) + cg3*std::pow(w[ijk+kk2],2.)) * dzi4[k  ])
                             + cg3*((cg0*std::pow(w[ijk    ],2.) + cg1*std::pow(w[ijk+kk1],2.) + cg2*std::pow(w[ijk+kk2],2.) + cg3*std::pow(w[ijk+kk3],2.)) * dzi4[k+1]) ) * dzhi4[k];
      }

  for(int k=grid->kstart; k<grid->kend; k++)
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

  for(int k=grid->kstart; k<grid->kend; k++)
  {
    u2_diss [k] = 0.;
    v2_diss [k] = 0.;
    w2_diss [k] = 0.;
    tke_diss[k] = 0.;

    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
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
  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
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

  for(int k=grid->kstart; k<grid->kend; k++)
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
  for(int k=grid->kstart; k<grid->kend; k++)
  {
    u2_rdstr [k] = 0.;
    v2_rdstr [k] = 0.;
    w2_rdstr [k] = 0.;

    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
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
  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj1 + k*kk1;
        w2_rdstr[k] += 2.*(ci0*p[ijk-kk2] + ci1*p[ijk-kk1] + ci2*p[ijk] + ci3*p[ijk+kk1])*
                       ( cg0*(ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ])
                       + cg1*(ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1])
                       + cg2*(ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])
                       + cg3*(ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3]) ) * dzhi4[k];
      }

  for(int k=grid->kstart; k<grid->kend; k++)
  {
    u2_rdstr [k] /= n;
    v2_rdstr [k] /= n;
    w2_rdstr [k] /= n;
  }

  grid->getprof(u2_rdstr , grid->kcells);
  grid->getprof(v2_rdstr , grid->kcells);
  grid->getprof(w2_rdstr , grid->kcells);


  // calculate the buoyancy term
  // CvH, check the correct usage of the gravity term later! check also the reference b!
  for(int k=grid->kstart; k<grid->kend; k++)
  {
    w2_buoy [k] = 0.;
    tke_buoy[k] = 0.;

    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj1 + k*kk1;
        tke_buoy[k] += (ci0*w[ijk-kk1] + ci1*w[ijk] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])*b[ijk];
      }
  }
  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj1 + k*kk1;
        w2_buoy[k] += 2.*(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk] + ci3*b[ijk+kk1])*w[ijk];
      }

  for(int k=grid->kstart; k<grid->kend; k++)
  {
    w2_buoy [k] /= n;
    tke_buoy[k] /= n;
  }

  grid->getprof(w2_buoy , grid->kcells);
  grid->getprof(tke_buoy, grid->kcells);

  return 0;
}

int cstats_dns::addfluxes(double * restrict flux, double * restrict turb, double * restrict diff)
{
  for(int k=grid->kstart; k<grid->kend+1; ++k)
    flux[k] = turb[k] + diff[k];

  return 0;
}

