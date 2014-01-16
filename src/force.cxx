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
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "force.h"
#include "defines.h"
#include "model.h"

cforce::cforce(cmodel *modelin)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  master = model->master;

  allocated = false;
}

cforce::~cforce()
{
  if(allocated)
  {
    if(swlspres == "geo")
    {
      delete[] ug;
      delete[] vg;
    }

    if(swls == "1")
    {
      for(std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
        delete[] lsprofs[*it];
    }

    if(swwls == "1")
      delete[] wls;
  }
}

int cforce::readinifile(cinput *inputin)
{
  int nerror = 0;

  nerror += inputin->getItem(&swlspres, "force", "swlspres", "", "0");
  nerror += inputin->getItem(&swls    , "force", "swls"    , "", "0");
  nerror += inputin->getItem(&swwls   , "force", "swwls"   , "", "0");
 
  if(swlspres != "0")
  {
    if(swlspres == "uflux")
      nerror += inputin->getItem(&uflux, "force", "uflux", "");
    else if(swlspres == "geo")
      nerror += inputin->getItem(&fc, "force", "fc", "");
    else
    {
      ++nerror;
      if(master->mpiid == 0) std::printf("ERROR \"%s\" is an illegal option for swlspres\n", swlspres.c_str());
    }
  }

  if(swls == "1")
    nerror += inputin->getList(&lslist, "force", "lslist", "");
  else if(swls != "0")
  {
    ++nerror;
    if(master->mpiid == 0) std::printf("ERROR \"%s\" is an illegal option for swls\n", swls.c_str());
  }

  if(swwls == "1")
    fields->setcalcprofs(true);
  else if(swwls != "0")
  {
    ++nerror;
    if(master->mpiid == 0) std::printf("ERROR \"%s\" is an illegal option for swwls\n", swwls.c_str());
  }

  return nerror;
}

int cforce::init()
{
  if(swlspres == "geo")
  {
    ug = new double[grid->kcells];
    vg = new double[grid->kcells];
  }

  if(swls == "1")
  {
    for(std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
      lsprofs[*it] = new double[grid->kcells];
  }

  if(swwls == "1")
    wls = new double[grid->kcells];

  allocated = true;

  return 0;
}

int cforce::create(cinput *inputin)
{
  int nerror = 0;

  if(swlspres == "geo")
  {
    nerror += inputin->getProf(&ug[grid->kstart], "ug", grid->kmax);
    nerror += inputin->getProf(&vg[grid->kstart], "vg", grid->kmax);
  }

  if(swls == "1")
  {
    // read the large scale sources, which are the variable names with a "ls" suffix
    for(std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
      nerror += inputin->getProf(&lsprofs[*it][grid->kstart], *it+"ls", grid->kmax);
  }

  if(swwls == "1")
    nerror += inputin->getProf(&wls[grid->kstart], "wls", grid->kmax);

  if(nerror > 0)
    return 1;

  return 0;
}

int cforce::exec(double dt)
{
  if(swlspres == "uflux")
    flux(fields->ut->data, fields->u->data, grid->dz, dt);

  else if(swlspres == "geo")
  {
    if(grid->swspatialorder == "2")
      coriolis_2nd(fields->ut->data, fields->vt->data, fields->u->data, fields->v->data, ug, vg);
    else if(grid->swspatialorder == "4")
      coriolis_4th(fields->ut->data, fields->vt->data, fields->u->data, fields->v->data, ug, vg);
  }

  if(swls == "1")
  {
    for(std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
      lssource(fields->st[*it]->data, lsprofs[*it]);
  }

  if(swwls == "1")
  {
    for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
      advecwls_2nd(it->second->data, fields->s[it->first]->datamean, wls, grid->dzhi);
  }

  return 0;
}

int cforce::flux(double * const restrict ut, const double * const restrict u, 
                 const double * const restrict dz, const double dt)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  
  double uavg, utavg, ugrid;

  uavg  = 0.;
  utavg = 0.;
  ugrid = grid->utrans;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        uavg  = uavg  + u [ijk]*dz[k];
        utavg = utavg + ut[ijk]*dz[k];
      }

  grid->getsum(&uavg);
  grid->getsum(&utavg);

  uavg  = uavg  / (grid->itot*grid->jtot*grid->zsize);
  utavg = utavg / (grid->itot*grid->jtot*grid->zsize);

  double fbody; 
  fbody = (uflux - uavg - ugrid) / dt - utavg;

  for(int n=0; n<grid->ncells; n++)
    ut[n] += fbody;

  return 0;
}

int cforce::coriolis_2nd(double * const restrict ut, double * const restrict vt,
                         const double * const restrict u , const double * const restrict v ,
                         const double * const restrict ug, const double * const restrict vg)
{
  int ijk,ii,jj,kk;
  double ugrid, vgrid;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  ugrid = grid->utrans;
  vgrid = grid->vtrans;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        ut[ijk] += fc * (0.25*(v[ijk-ii] + v[ijk] + v[ijk-ii+jj] + v[ijk+jj]) + vgrid - vg[k]);
      }

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        vt[ijk] -= fc * (0.25*(u[ijk-jj] + u[ijk] + u[ijk+ii-jj] + u[ijk+ii]) + ugrid - ug[k]);
      }

  return 0;
}

int cforce::coriolis_4th(double * const restrict ut, double * const restrict vt,
                         const double * const restrict u , const double * const restrict v ,
                         const double * const restrict ug, const double * const restrict vg)
{
  int ijk,ii1,ii2,jj1,jj2,kk1;
  double ugrid, vgrid;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;

  ugrid = grid->utrans;
  vgrid = grid->vtrans;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj1 + k*kk1;
        ut[ijk] += fc * ( ( ci0*(ci0*v[ijk-ii2-jj1] + ci1*v[ijk-ii1-jj1] + ci2*v[ijk-jj1] + ci3*v[ijk+ii1-jj1])
                          + ci1*(ci0*v[ijk-ii2    ] + ci1*v[ijk-ii1    ] + ci2*v[ijk    ] + ci3*v[ijk+ii1    ])
                          + ci2*(ci0*v[ijk-ii2+jj1] + ci1*v[ijk-ii1+jj1] + ci2*v[ijk+jj1] + ci3*v[ijk+ii1+jj1])
                          + ci3*(ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii1+jj2] + ci2*v[ijk+jj2] + ci3*v[ijk+ii1+jj2]) )
                        + vgrid - vg[k] );
      }

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj1 + k*kk1;
        vt[ijk] -= fc * ( ( ci0*(ci0*u[ijk-ii1-jj2] + ci1*u[ijk-jj2] + ci2*u[ijk+ii1-jj2] + ci3*u[ijk+ii2-jj2])
                          + ci1*(ci0*u[ijk-ii1-jj1] + ci1*u[ijk-jj1] + ci2*u[ijk+ii1-jj1] + ci3*u[ijk+ii2-jj1])
                          + ci2*(ci0*u[ijk-ii1    ] + ci1*u[ijk    ] + ci2*u[ijk+ii1    ] + ci3*u[ijk+ii2    ])
                          + ci3*(ci0*u[ijk-ii1+jj1] + ci1*u[ijk+jj1] + ci2*u[ijk+ii1+jj1] + ci3*u[ijk+ii2+jj1]) )
                        + ugrid - ug[k]);
      }

  return 0;
}

int cforce::lssource(double * const restrict st, const double * const restrict sls)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        st[ijk] += sls[k];
      }

  return 0;
}

int cforce::advecwls_2nd(double * const restrict st, const double * const restrict s,
                         const double * const restrict wls, const double * const dzhi)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->ijcells;

  // use an upwind differentiation
  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    if(wls[k] > 0.)
    {
      for(int j=grid->jstart; j<grid->jend; ++j)
        for(int i=grid->istart; i<grid->iend; ++i)
        {
          ijk = i + j*jj + k*kk;
          st[ijk] -=  wls[k] * (s[k]-s[k-1])*dzhi[k];
        }
    }
    else
    {
      for(int j=grid->jstart; j<grid->jend; ++j)
        for(int i=grid->istart; i<grid->iend; ++i)
        {
          ijk = i + j*jj + k*kk;
          st[ijk] -=  wls[k] * (s[k+1]-s[k])*dzhi[k+1];
        }
    }
  }

  return 0;
}

inline double cforce::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}
