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
#include "stdlib.h"
#include <cmath>
#include "grid.h"
#include "fields.h"
#include "thermo_moist.h"
#include "defines.h"

cthermo_moist::cthermo_moist(cmodel *modelin) : cthermo(modelin)
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

  return (nerror > 0);
}

int cthermo_moist::create()
{
  int nerror = 0;
  
  thvs = 303.2;  //ssurf * (1. - (1. - rv/rd)*qtsurf);

  pmn  = new double[grid->kcells];  // Hydrostatic pressure (full levels)
 
  allocated = true;
  return nerror;
}

int cthermo_moist::exec()
{
  int kk,nerror;
  kk = grid->icells*grid->jcells;

  nerror = 0;

  // extend later for gravity vector not normal to surface
  if(grid->swspatialorder == "2")
  {
    calchydropres_2nd(pmn,fields->s["s"]->data,fields->s["s"]->datamean,fields->s["qt"]->data,fields->s["qt"]->datamean);
    calcbuoyancytend_2nd(fields->wt->data, fields->s["s"]->data, fields->s["qt"]->data, pmn,
                         &fields->s["tmp2"]->data[0*kk], &fields->s["tmp2"]->data[1*kk], &fields->s["tmp2"]->data[2*kk]);
  }
  else if(grid->swspatialorder == "4")
  {
    calchydropres_4th(pmn,fields->s["s"]->data,fields->s["s"]->datamean,fields->s["qt"]->data,fields->s["qt"]->datamean);
    calcbuoyancytend_4th(fields->wt->data, fields->s["s"]->data, fields->s["qt"]->data, pmn,
                         &fields->s["tmp2"]->data[0*kk], &fields->s["tmp2"]->data[1*kk], &fields->s["tmp2"]->data[2*kk]);
  }

  return (nerror>0);
}

int cthermo_moist::checkthermofield(std::string name)
{
  if(name == "b" || name == "ql")
    return 0;
  else
    return 1;
}

int cthermo_moist::getthermofield(cfield3d *field, cfield3d *tmp, std::string name)
{
  // calculate the hydrostatic pressure
  if(grid->swspatialorder == "2")
    calchydropres_2nd(pmn,fields->s["s"]->data,fields->s["s"]->datamean,fields->s["qt"]->data,fields->s["qt"]->datamean);
  else if(grid->swspatialorder == "4")
    calchydropres_4th(pmn,fields->s["s"]->data,fields->s["s"]->datamean,fields->s["qt"]->data,fields->s["qt"]->datamean);

  if(name == "b")
    calcbuoyancy(field->data, fields->s["s"]->data, fields->s["qt"]->data, pmn, tmp->data);
  else if(name == "ql")
    calcqlfield(field->data, fields->s["s"]->data, fields->s["qt"]->data, pmn);

  return 0;
}

//int cthermo_moist::getql(cfield3d *qlfield, cfield3d *pfield)
//{
//  // calculate the hydrostatic pressure
//  if(grid->swspatialorder == "2")
//  {
//    calchydropres_2nd(pmn,fields->s["s"]->data,fields->s["s"]->datamean,fields->s["qt"]->data,fields->s["qt"]->datamean);
//  }
//  else if(grid->swspatialorder == "4")
//  {
//    calchydropres_4th(pmn,fields->s["s"]->data,fields->s["s"]->datamean,fields->s["qt"]->data,fields->s["qt"]->datamean);
//  }
//
//  // calculate the ql field
//  calcqlfield(qlfield->data, fields->s["s"]->data, fields->s["qt"]->data, pmn);
//
//  return 0;
//}
//
//int cthermo_moist::getbuoyancy(cfield3d *bfield, cfield3d *tmp)
//{
//  // calculate the hydrostatic pressure
//  if(grid->swspatialorder == "2")
//  {
//    calchydropres_2nd(pmn,fields->s["s"]->data,fields->s["s"]->datamean,fields->s["qt"]->data,fields->s["qt"]->datamean);
//  }
//  else if(grid->swspatialorder == "4")
//  {
//    calchydropres_4th(pmn,fields->s["s"]->data,fields->s["s"]->datamean,fields->s["qt"]->data,fields->s["qt"]->datamean);
//  }
//
//  // calculate the buoyancy at the cell centers
//  calcbuoyancy(bfield->data, fields->s["s"]->data, fields->s["qt"]->data, pmn, fields->s["tmp2"]->data);
//
//  return 0;
//}

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

/**
 * This function calculates the hydrostatic pressure
 * Solves: dpi/dz=-g/thv with pi=cp*(p/p0)**(rd/cp)
 * @param pmn Pointer to hydrostatic pressure array
 * @param s,smean,qt,qtmean .... 
 * @return Returns 1 on error, 0 otherwise.
 */
int cthermo_moist::calchydropres_2nd(double * restrict pmn, double * restrict s, double * restrict smean,
                                 double * restrict qt, double * restrict qtmean)
{
  int kstart,kend;
  double thv,ssurf,qtsurf,stop,qttop,ptop;
  double rdcp = rd/cp;

  kstart = grid->kstart;
  kend = grid->kend;

  // Calculate horizontal mean profiles, and interpolate surface and model top values 
  grid->calcmean(smean,s,grid->kcells);
  grid->calcmean(qtmean,qt,grid->kcells);

  ssurf  = interp2(smean[kstart-1], smean[kstart]);
  stop   = interp2(smean[kend-1],   smean[kend]);
  qtsurf = interp2(qtmean[kstart-1],qtmean[kstart]);
  qttop  = interp2(qtmean[kend-1],  qtmean[kend]);

  // Calculate lowest full level (kstart) from surface values p,s,qt
  thv = ssurf*(1.+(rv/rd-1)*qtsurf);
  pmn[kstart] = pow((pow(ps,rdcp) - grav * pow(p0,rdcp) * grid->z[kstart] / (cp * thv)),(1./rdcp)); 

  for(int k=kstart+1; k<kend; k++)
  {
    thv = interp2(smean[k-1],smean[k])*(1.+(rv/rd-1.)*interp2(qtmean[k-1],qtmean[k]));   // BvS: assume no ql for now..
    pmn[k] = pow((pow(pmn[k-1],rdcp) - grav * pow(p0,rdcp) * grid->dzh[k] / (cp * thv)),(1./rdcp)); 
  }

  // Calculate pressure at top of domain, needed to fill ghost cells
  thv = stop*(1.+(rv/rd-1)*qttop);
  ptop = pow((pow(pmn[kend-1],rdcp) - grav * pow(p0,rdcp) * (grid->zh[kend]-grid->z[kend-1]) / (cp * thv)),(1./rdcp));

  // Fill bottom and top ghost cells 
  pmn[kstart-1] = 2.*ps - pmn[kstart];
  pmn[kend] = 2.*ptop - pmn[kend-1];

  return 0;
}

/**
 * This function calculates the hydrostatic pressure
 * Solves: dpi/dz=-g/thv with pi=cp*(p/p0)**(rd/cp)
 * @param pmn Pointer to hydrostatic pressure array
 * @param s,smean,qt,qtmean .... 
 * @return Returns 1 on error, 0 otherwise.
 */
int cthermo_moist::calchydropres_4th(double * restrict pmn, double * restrict s, double * restrict smean,
                                 double * restrict qt, double * restrict qtmean)
{
  int kstart,kend;
  double thv,ssurf,qtsurf,stop,qttop,ptop;
  double rdcp = rd/cp;

  kstart = grid->kstart;
  kend = grid->kend;

  // Calculate horizontal mean profiles, and interpolate surface and model top values 
  grid->calcmean(smean,s,grid->kcells);
  grid->calcmean(qtmean,qt,grid->kcells);

  ssurf  = interp4(smean[kstart-2], smean[kstart-1], smean[kstart], smean[kstart+1]);
  stop   = interp4(smean[kend-2],   smean[kend-1],   smean[kend],   smean[kend+1]);
  qtsurf = interp4(qtmean[kstart-2],qtmean[kstart-1],qtmean[kstart],qtmean[kstart+1]);
  qttop  = interp4(qtmean[kend-2],  qtmean[kend-1],  qtmean[kend],  qtmean[kend+1]);

  // Calculate lowest full level (kstart) from surface values p,s,qt
  thv = ssurf*(1.+(rv/rd-1)*qtsurf);
  pmn[kstart] = pow((pow(ps,rdcp) - grav * pow(p0,rdcp) * grid->z[kstart] / (cp * thv)),(1./rdcp)); 

  for(int k=kstart+1; k<kend; k++)
  {
    thv = interp4(smean[k-2],smean[k-1],smean[k],smean[k+1])*(1.+(rv/rd-1.)*interp4(qtmean[k-2],qtmean[k-1],qtmean[k],qtmean[k+1]));   // BvS: assume no ql for now..
    pmn[k] = pow((pow(pmn[k-1],rdcp) - grav * pow(p0,rdcp) * grid->dzh[k] / (cp * thv)),(1./rdcp)); 
  }

  // Calculate pressure at top of domain, needed to fill ghost cells
  thv = stop*(1.+(rv/rd-1)*qttop);
  ptop = pow((pow(pmn[kend-1],rdcp) - grav * pow(p0,rdcp) * (grid->zh[kend]-grid->z[kend-1]) / (cp * thv)),(1./rdcp));

  // Fill bottom and top ghost cells 
  pmn[kstart-1] = (8./3.)*ps - 2.*pmn[kstart] + (1./3.)*pmn[kstart+1];
  pmn[kstart-2] = 8.*ps - 9.*pmn[kstart] + 2.*pmn[kstart+1];
  pmn[kend] = (8./3.)*ptop - 2.*pmn[kend-1] + (1./3.)*pmn[kend-2];
  pmn[kend+1] = 8.*ptop - 9.*pmn[kend-1] + 2.*pmn[kend-2];

  return 0;
}


int cthermo_moist::calcbuoyancytend_2nd(double * restrict wt, double * restrict s, double * restrict qt, double * restrict p,
                                        double * restrict sh, double * restrict qth, double * restrict ql)
{
  int ijk,jj,kk,ij;
  double tl, ph, exnh;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double thvref = thvs;

  // CvH check the usage of the gravity term here, in case of scaled DNS we use one. But thermal expansion coeff??
  for(int k=grid->kstart+1; k<grid->kend; k++)
  {
    ph   = interp2(p[k-1],p[k]);   // BvS To-do: calculate pressure at full and half levels
    exnh = exner2(ph);
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ij  = i + j*jj;
        sh[ij]  = interp2(s[ijk-kk], s[ijk]);
        qth[ij] = interp2(qt[ijk-kk], qt[ijk]);
        tl      = sh[ij] * exnh;
        // Calculate first estimate of ql using Tl
        // if ql(Tl)>0, saturation adjustment routine needed
        ql[ij]  = qth[ij]-rslf(ph,tl);
      }
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ij  = i + j*jj;
        if(ql[ij]>0)   // already doesn't vectorize because of iteration in calcql()
          ql[ij] = calcql(sh[ij], qth[ij], ph, exnh);
        else
          ql[ij] = 0.;
      }
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ij  = i + j*jj;
        wt[ijk] += bu(ph, sh[ij], qth[ij], ql[ij], thvref);
      }
  }
  return 0;
}

int cthermo_moist::calcbuoyancytend_4th(double * restrict wt, double * restrict s, double * restrict qt, double * restrict p,
                                        double * restrict sh, double * restrict qth, double * restrict ql)
{
  int ijk,jj,ij;
  int kk1,kk2;
  double tl, ph, exnh;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  double thvref = thvs;

  for(int k=grid->kstart+1; k<grid->kend; k++)
  {
    ph  = interp4(p[k-2] , p[k-1] , p[k] , p[k+1]); // BvS To-do: calculate pressure at full and half levels
    exnh = exner2(ph);
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk1;
        ij  = i + j*jj;
        sh[ij]  = interp4(s[ijk-kk2] , s[ijk-kk1] , s[ijk] , s[ijk+kk1]);
        qth[ij] = interp4(qt[ijk-kk2], qt[ijk-kk1], qt[ijk], qt[ijk+kk1]);
        tl      = sh[ij] * exnh;
        // Calculate first estimate of ql using Tl
        // if ql(Tl)>0, saturation adjustment routine needed
        ql[ij]  = qth[ij]-rslf(ph,tl);   
      }
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ij  = i + j*jj;
        if(ql[ij]>0)   // already doesn't vectorize because of iteration in calcql()
          ql[ij] = calcql(sh[ij], qth[ij], ph, exnh);
        else
          ql[ij] = 0.;
      }
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk1;
        ij  = i + j*jj;
        wt[ijk] += bu(ph, sh[ij], qth[ij], ql[ij], thvref);
      }
  }
  return 0;
}

int cthermo_moist::calcbuoyancy(double * restrict b, double * restrict s, double * restrict qt, double * restrict p, double * restrict ql)
{
  int ijk,jj,kk,ij;
  double tl, exn;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double thvref = thvs;

  for(int k=0; k<grid->kcells; k++)
  {
    exn = exner2(p[k]);
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ij  = i + j*jj;
        tl  = s[ijk] * exn;
        ql[ij]  = qt[ijk]-rslf(p[k],tl);   // not real ql, just estimate
      }
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ij  = i + j*jj;
        if(ql[ij] > 0)
          ql[ij] = calcql(s[ijk], qt[ijk], p[k], exn);
        else
          ql[ij] = 0.;
      }
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ij  = i + j*jj;
        b[ijk] = bu(p[k], s[ijk], qt[ijk], ql[ij], thvref);
      }
  }

  return 0;
}

int cthermo_moist::calcqlfield(double * restrict ql, double * restrict s, double * restrict qt, double * restrict p)
{
  int ijk,jj,kk;
  double exn;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int k=grid->kstart; k<grid->kend; k++)
  {
    exn = exner2(p[k]);
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ql[ijk] = calcql(s[ijk], qt[ijk], p[k], exn);
      }
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
  int ij,jj;
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
  return grav * ((s + lv*ql/(cp*exner2(p))) * (1. - (1. - rv/rd)*qt - rv/rd*ql) - thvref) / thvref;
}

inline double cthermo_moist::bunoql(const double s, const double qt, const double thvref)
{
  return grav * (s * (1. - (1. - rv/rd)*qt) - thvref) / thvref;
}

inline double cthermo_moist::bufluxnoql(const double s, const double sflux, const double qt, const double qtflux, const double thvref)
{
  return grav/thvref * (sflux * (1. - (1.-rv/rd)*qt) - (1.-rv/rd)*s*qtflux);
}

inline double cthermo_moist::calcql(const double s, const double qt, const double p, const double exn)
{
  int niter = 0; //, nitermax = 5;
  double ql, tl, tnr_old = 1.e9, tnr, qs;
  tl = s * exn;
  tnr = tl;
  while (std::fabs(tnr-tnr_old)/tnr_old> 1e-5)// && niter < nitermax)
  {
    ++niter;
    tnr_old = tnr;
    qs = rslf(p,tnr);
    tnr = tnr - (tnr+(lv/cp)*qs-tl-(lv/cp)*qt)/(1+(std::pow(lv,2)*qs)/ (rv*cp*std::pow(tnr,2)));
  }
  ql = std::max(0.,qt - qs);
  return ql;
}

inline double cthermo_moist::exner(const double p)
{
  return pow((p/p0),(rd/cp));
}

inline double cthermo_moist::exner2(const double p)
{
  double dp=p-p0;
  return (1+(dp*(ex1+dp*(ex2+dp*(ex3+dp*(ex4+dp*(ex5+dp*(ex6+ex7*dp)))))))); 
}

inline double cthermo_moist::rslf(const double p, const double t)
{
  return ep*esl(t)/(p-(1-ep)*esl(t));
}

inline double cthermo_moist::esl(const double t)
{
  const double x=std::max(-80.,t-tmelt);
  return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))));

  //return es0*std::exp(at*(t-tmelt)/(t-bt));
}

inline double cthermo_moist::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}

inline double cthermo_moist::interp4(const double a, const double b, const double c, const double d)
{
  return (-a + 9.*b + 9.*c - d) / 16.;
}


