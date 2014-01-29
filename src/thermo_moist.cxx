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
#include "grid.h"
#include "fields.h"
#include "thermo_moist.h"
#include "defines.h"

cthermo_moist::cthermo_moist(cgrid *gridin, cfields *fieldsin, cmaster *mpiin) : cthermo(gridin, fieldsin, mpiin)
{
  swthermo = "moist";

  allocated = false;
}

cthermo_moist::~cthermo_moist()
{
  if (allocated)
  {
    delete[] pmn;
  }
}

int cthermo_moist::readinifile(cinput *inputin)
{
  int nerror = 0;
  nerror += inputin->getItem(&ps, "thermo", "ps", "");

  nerror += fields->initpfld("s");
  nerror += inputin->getItem(&fields->sp["s"]->visc, "fields", "svisc", "s");
  nerror += fields->initpfld("qt");
  nerror += inputin->getItem(&fields->sp["qt"]->visc, "fields", "svisc", "qt");
  // nerror += fields->initdfld("ql");

  return (nerror > 0);
}

int cthermo_moist::create()
{
  int nerror = 0;
  double qtsurf, ssurf;
  
  // Create hydrostatic profile
  // Doesnt work; just put a ref value in Namelist for both gravity and 1/beta
//   nerror += grid->calcmean(&ssurf, fields->s["s"]->databot,1);
//   nerror += grid->calcmean(&qtsurf, fields->s["qt"]->databot,1);

  thvs = 303.2;//ssurf * (1. - (1. - rv/rd)*qtsurf);
  double tvs  = exner(ps) * thvs;
  rhos = ps / (rd * tvs);

  pmn = new double[grid->kcells];
  for(int k=0; k<grid->kcells; k++)
  {
    pmn[k] = ps - rhos * grav * grid->z[k];
  }

  allocated = true;
  return nerror;
}

int cthermo_moist::exec()
{
  int ijk,kk,nerror;
  kk = grid->icells*grid->jcells;

  nerror = 0;

  calcpres(fields->s["tmp1"]->data, fields->s["p"]->data, pmn);

  // nerror += calcqlfield(fields->s["ql"]->data, fields->s["s"]->data, fields->s["qt"]->data, fields->s["tmp1"]->data);

  // extend later for gravity vector not normal to surface
  if(grid->swspatialorder == "2")
  {
    calcbuoyancytend_2nd(fields->wt->data, fields->s["s"]->data, fields->s["qt"]->data, fields->s["tmp1"]->data);
  }
  else if(grid->swspatialorder == "4")
  {
    calcbuoyancytend_4th(fields->wt->data, fields->s["s"]->data, fields->s["qt"]->data, fields->s["tmp1"]->data);
  }


  return (nerror>0);
}

int cthermo_moist::getql(cfield3d *qlfield, cfield3d *pfield)
{
  int ijk,kk;
  kk = grid->icells*grid->jcells;

  // calculate the hydrostatic pressure
  calcpres(pfield->data, fields->s["p"]->data, this->pmn);

  // calculate the ql field
  calcqlfield(qlfield->data, fields->s["s"]->data, fields->s["qt"]->data, pfield->data);

  return 0;
}

int cthermo_moist::getbuoyancy(cfield3d *bfield, cfield3d *tmp)
{
  // first calculate the pressure
  calcpres(tmp->data, fields->s["p"]->data, this->pmn);
  // then calculate the buoyancy at the cell centers
  calcbuoyancy(bfield->data, fields->s["s"]->data, fields->s["qt"]->data, tmp->data);

  return 0;
}

int cthermo_moist::getbuoyancysurf(cfield3d *bfield)
{
  calcbuoyancybot(bfield->data         , bfield->databot,
                  fields->s["s" ]->data, fields->s["s" ]->databot,
                  fields->s["qt"]->data, fields->s["qt"]->databot);
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->s["s"]->databot, fields->s["s"]->datafluxbot, fields->s["qt"]->databot, fields->s["qt"]->datafluxbot);
  return 0;
}

int cthermo_moist::getbuoyancyfluxbot(cfield3d *bfield)
{
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->s["s"]->databot, fields->s["s"]->datafluxbot, fields->s["qt"]->databot, fields->s["qt"]->datafluxbot);
  return 0;
}

int cthermo_moist::calcpres(double * restrict p, double * restrict pi, double * restrict pmn)
{
  int ijk,jj,kk;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double rhos = this->rhos;

  for(int k=0; k<grid->kcells; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        p[ijk] = pi[ijk]*rhos + pmn[k]; // minus trace
      }

  return 0;
}

int cthermo_moist::calcbuoyancytend_2nd(double * restrict wt, double * restrict s, double * restrict qt, double * restrict p)
{
  int ijk,jj,kk;
  double sh, qth, ph, ql;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double thvref = thvs;

  // CvH check the usage of the gravity term here, in case of scaled DNS we use one. But thermal expansion coeff??
  for(int k=grid->kstart+1; k<grid->kend; k++)
  {
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        sh  = interp2(s[ijk-kk], s[ijk]);
        qth = interp2(qt[ijk-kk], qt[ijk]);
        ph  = interp2(p[ijk-kk], p[ijk]);
        ql  = calcql(sh, qth, ph);
        wt[ijk] += bu(ph, sh, qth, ql, thvref);
        // if(ql>0)
        //  std::printf("%d qt %f rs %f ql %f bu %f bu0 %f %f\n", k, qth, rslf(ph, sh*exner(ph)),ql, thvs/grav*(bu(ph, sh,qth,ql)), thvs/grav*(bu(ph, sh,qth,0.)), ph);
      }
//           printf("%d qt %f rs %f ql %f bu %f bu0 %f %f\n", k, qth, qth/rslf(ph, sh*exner(ph)),ql, thvs/grav*(bu(ph, sh,qth,ql)), thvs/grav*(bu(ph, sh,qth,0.)), ph);

  }
  return 0;
}

int cthermo_moist::calcbuoyancytend_4th(double * restrict wt, double * restrict s, double * restrict qt, double * restrict p)
{
  int ijk,jj;
  int kk1,kk2;
  double sh, qth, ph, ql;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  double thvref = thvs;

  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk1;
        sh  = interp4(s[ijk-kk2] , s[ijk-kk1] , s[ijk] , s[ijk+kk1]);
        qth = interp4(qt[ijk-kk2], qt[ijk-kk1], qt[ijk], qt[ijk+kk1]);
        ph  = interp4(p[ijk-kk2] , p[ijk-kk1] , p[ijk] , p[ijk+kk1]);
        ql  = calcql(sh, qth, ph);
        wt[ijk] += bu(ph, sh, qth, ql, thvref);
      }

  return 0;
}

int cthermo_moist::calcbuoyancy(double * restrict b, double * restrict s, double * restrict qt, double * restrict p)
{
  int ijk,jj,kk;
  double ql;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double thvref = thvs;

  for(int k=0; k<grid->kcells; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ql = calcql(s[ijk], qt[ijk], p[ijk]);
        b[ijk] = bu(p[ijk], s[ijk], qt[ijk], ql, thvref);
      }

  return 0;
}

int cthermo_moist::calcqlfield(double * restrict ql, double * restrict s, double * restrict qt, double * restrict p)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ql[ijk] = calcql(s[ijk], qt[ijk], p[ijk]);
      }

  return 0;
}

int cthermo_moist::calcbuoyancybot(double * restrict b , double * restrict bbot,
                                   double * restrict s , double * restrict sbot,
                                   double * restrict qt, double * restrict qtbot)
{
  int ij,ijk,jj,kk,kstart;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  kstart = grid->kstart;

  double thvref = thvs;

  // assume no liquid water at the lowest model level
  for(int j=0; j<grid->jcells; j++)
#pragma ivdep
    for(int i=0; i<grid->icells; i++)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      bbot[ij ] = bunoql(sbot[ij], qtbot[ij], thvref);
      b   [ijk] = bunoql(s[ijk], qt[ijk], thvref);
    }

  return 0;
}

int cthermo_moist::calcbuoyancyfluxbot(double * restrict bfluxbot, double * restrict sbot, double * restrict sfluxbot, double * restrict qtbot, double * restrict qtfluxbot)
{
  int ij,jj,kk;
  jj = grid->icells;

  double thvref = thvs;

  // assume no liquid water at the lowest model level
  for(int j=0; j<grid->jcells; j++)
#pragma ivdep
    for(int i=0; i<grid->icells; i++)
    {
      ij  = i + j*jj;
      bfluxbot[ij] = bufluxnoql(sbot[ij], sfluxbot[ij], qtbot[ij], qtfluxbot[ij], thvref);
    }

  return 0;
}

// INLINE FUNCTIONS
inline double cthermo_moist::bu(const double p, const double s, const double qt, const double ql, const double thvref)
{
  return grav * ((s + lv*ql/(cp*exner(p))) * (1. - (1. - rv/rd)*qt - rv/rd*ql) - thvref) / thvref;
}

inline double cthermo_moist::bunoql(const double s, const double qt, const double thvref)
{
  return grav * (s * (1. - (1. - rv/rd)*qt) - thvref) / thvref;
}

inline double cthermo_moist::bufluxnoql(const double s, const double sflux, const double qt, const double qtflux, const double thvref)
{
  return grav/thvref * (sflux * (1. - (1.-rv/rd)*qt) - (1.-rv/rd)*s*qtflux);
}

inline double cthermo_moist::calcql(const double s, const double qt, const double p)
{
  int niter = 0, nitermax = 5;
//   double sabs, sguess = 1.e9, t, qs, ql, dtldt;
  double ql, tl, tnr_old = 1.e9, tnr, qs;
  tl = s * exner(p);
  tnr = tl;
  while (std::fabs(tnr-tnr_old)/tnr_old> 1e-5)// && niter < nitermax)
  {
    ++niter;
    tnr_old = tnr;
//     es      = es0*exp(at*(Tnr-tmelt)/(Tnr-bt));
//     qs   = rd/rv*es/(pressure(k)-(1-rd/rv)*es);
    qs = rslf(p,tnr);
//     printf("%f, %f %f %f\n",qs,qt, tnr, p);
    tnr     = tnr - (tnr+(lv/cp)*qs-tl-(lv/cp)*qt)/(1+(std::pow(lv,2)*qs)/ (rv*cp*std::pow(tnr,2)));

//     sguess = t*std::exp(-lv*ql/(cp*t));
//     sguess = t + lv*ql/cp;
//     dtldt = sabs/t*(1. - lv * ql / (pow(cp,2.)* t));
//     t += (sguess-sabs) / dtldt;
  }
  ql = std::max(0.,qt - qs);
  return ql;
}

inline double cthermo_moist::exner(const double p)
{
  return pow((p/p0),(rd/cp));
}
inline double cthermo_moist::rslf(const double p, const double t)
{
  return ep*esl(t)/(p-(1-ep)*esl(t));
}

inline double cthermo_moist::esl(const double t)
{
  const double c0=0.6105851e+03;
  const double c1=0.4440316e+02;
  const double c2=0.1430341e+01;
  const double c3=0.2641412e-01;
  const double c4=0.2995057e-03;
  const double c5=0.2031998e-05;
  const double c6=0.6936113e-08;
  const double c7=0.2564861e-11;
  const double c8=-.3704404e-13;
//   const double x=std::max(-80.,t-tmelt);
  const double at= 17.27;
  const double bt= 35.86;
  const double es0 = 610.78;
//   return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))));
  return es0*std::exp(at*(t-tmelt)/(t-bt));
}

inline double cthermo_moist::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}

inline double cthermo_moist::interp4(const double a, const double b, const double c, const double d)
{
  return (-a + 9.*b + 9.*c - d) / 16.;
}

