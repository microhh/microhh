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
#include <algorithm>
#include <cmath>
#include "grid.h"
#include "fields.h"
#include "thermo_dry.h"
#include "defines.h"
#include "model.h"
#include "stats.h"
#include "diff_les2s.h"
#include "master.h"
#include "cross.h"

#define grav 9.81
#define Rd 287.04
#define cp 1005.

#define NO_OFFSET 0.

cthermo_dry::cthermo_dry(cmodel *modelin, cinput *inputin) : cthermo(modelin, inputin)
{
  swthermo = "dry";

  thref  = 0;
  threfh = 0;
  pref   = 0;
  prefh  = 0;
  exner  = 0;
  exnerh = 0;

  int nerror = 0;
  nerror += inputin->getItem(&pbot, "thermo", "pbot", "");

  nerror += fields->initpfld("th", "Potential Temperature", "K");
  nerror += inputin->getItem(&fields->sp["th"]->visc, "fields", "svisc", "th");

  // Read list of cross sections
  nerror += inputin->getList(&crosslist , "thermo", "crosslist" , "");

  if(nerror)
    throw 1;
}

cthermo_dry::~cthermo_dry()
{
  delete[] this->thref;
  delete[] this->threfh;
  delete[] this->pref;
  delete[] this->prefh;
  delete[] this->exner;
  delete[] this->exnerh;
}

int cthermo_dry::init()
{
  // copy pointers
  stats = model->stats;

  // fields for Boussinesq and anelastic solver
  thref  = new double[grid->kcells];
  threfh = new double[grid->kcells];
  pref   = new double[grid->kcells];
  prefh  = new double[grid->kcells];
  exner  = new double[grid->kcells];
  exnerh = new double[grid->kcells];

  return 0;
}

int cthermo_dry::create(cinput *inputin)
{
  // Only in case of Boussinesq, read in reference potential temperature
  int nerror = 0;
  if(model->swbasestate == "boussinesq")
    nerror += inputin->getItem(&thref0, "thermo", "thref0", "");
  if(nerror)
    throw 1;

  // Setup base state for anelastic solver
  if(model->swbasestate == "anelastic")
  {
    // take the initial profile as the reference
    if(inputin->getProf(&thref[grid->kstart], "th", grid->kmax))
      return 1;

    int kstart = grid->kstart;
    int kend   = grid->kend;

    // extrapolate the profile to get the bottom value
    threfh[kstart] = thref[kstart] - grid->z[kstart]*(thref[kstart+1]-thref[kstart])*grid->dzhi[kstart+1];

    // extrapolate the profile to get the top value
    threfh[kend] = thref[kend-1] + (grid->zh[kend]-grid->z[kend-1])*(thref[kend-1]-thref[kend-2])*grid->dzhi[kend-1];

    // set the ghost cells for the reference temperature
    thref[kstart-1] = 2.*threfh[kstart] - thref[kstart];
    thref[kend]     = 2.*threfh[kend]   - thref[kend-1];

    // interpolate the reference temperature profile
    for(int k=grid->kstart+1; k<grid->kend; ++k)
      threfh[k] = 0.5*(thref[k-1] + thref[k]);

    // ANELASTIC
    // calculate the base state pressure and density
    for(int k=grid->kstart; k<grid->kend; ++k)
    {
      pref [k] = pbot*std::exp(-grav/(Rd*thref[k])*grid->z[k]);
      exner[k] = std::pow(pref[k]/pbot, Rd/cp);

      // set the base density for the entire model
      fields->rhoref[k] = pref[k] / (Rd*exner[k]*thref[k]);
    }

    for(int k=grid->kstart; k<grid->kend+1; ++k)
    {
      prefh [k] = pbot*std::exp(-grav/(Rd*threfh[k])*grid->zh[k]);
      exnerh[k] = std::pow(prefh[k]/pbot, Rd/cp);

      // set the base density for the entire model
      fields->rhorefh[k] = prefh[k] / (Rd*exnerh[k]*threfh[k]);
    }

    // set the ghost cells for the reference variables
    // CvH for now in 2nd order
    pref [kstart-1] = 2.*prefh [kstart] - pref [kstart];
    exner[kstart-1] = 2.*exnerh[kstart] - exner[kstart];
    fields->rhoref[kstart-1] = 2.*fields->rhorefh[kstart] - fields->rhoref[kstart];

    pref [kend] = 2.*prefh [kend] - pref [kend-1];
    exner[kend] = 2.*exnerh[kend] - exner[kend-1];
    fields->rhoref[kend] = 2.*fields->rhorefh[kend] - fields->rhoref[kend-1];
  }
  else
  {
    // Set entire column to reference value
    for(int k=0; k<grid->kcells; ++k)
    {
      thref[k]  = thref0;
      threfh[k] = thref0;
    }
  }
  
  // add variables to the statistics
  if(stats->getsw() == "1")
  {
    // Add base state profiles to statistics -> needed/wanted for Boussinesq? Or write as 0D var?
    //stats->addfixedprof("pref",    "Full level basic state pressure", "Pa",     "z",  pref);
    //stats->addfixedprof("prefh",   "Half level basic state pressure", "Pa",     "zh", prefh);
    stats->addfixedprof("rhoref",  "Full level basic state density",  "kg m-3", "z",  fields->rhoref);
    stats->addfixedprof("rhorefh", "Half level basic state density",  "kg m-3", "zh", fields->rhorefh);
    stats->addfixedprof("thref",   "Full level reference potential temperature", "K", "z",thref);
    stats->addfixedprof("threfh",  "Half level reference potential temperature", "K", "zh",thref);

    stats->addprof("b", "Buoyancy", "m s-2", "z");
    for(int n=2; n<5; ++n)
    {
      std::stringstream ss;
      ss << n;
      std::string sn = ss.str();
      stats->addprof("b"+sn, "Moment " +sn+" of the buoyancy", "(m s-2)"+sn,"z");
    }

    stats->addprof("bgrad", "Gradient of the buoyancy", "s-2", "zh");
    stats->addprof("bw"   , "Turbulent flux of the buoyancy", "m2 s-3", "zh");
    stats->addprof("bdiff", "Diffusive flux of the buoyancy", "m2 s-3", "zh");
    stats->addprof("bflux", "Total flux of the buoyancy", "m2 s-3", "zh");

    stats->addprof("bsort", "Sorted buoyancy", "m s-2", "z");
  }

  // Cross sections (isn't there an easier way to populate this list?)
  allowedcrossvars.push_back("b");
  allowedcrossvars.push_back("bbot");
  allowedcrossvars.push_back("bfluxbot");
  if(grid->swspatialorder == "4")
    allowedcrossvars.push_back("blngrad");

  // Check input list of cross variables (crosslist)
  std::vector<std::string>::iterator it=crosslist.begin();
  while(it != crosslist.end())
  {
    if(!std::count(allowedcrossvars.begin(),allowedcrossvars.end(),*it))
    {
      if(master->mpiid == 0) std::printf("WARNING field %s in [thermo][crosslist] is illegal\n", it->c_str());
      it = crosslist.erase(it);  // erase() returns iterator of next element..
    }
    else
      ++it;
  }

  // Sort crosslist to group ql and b variables
  std::sort(crosslist.begin(),crosslist.end());

  return 0;
}

int cthermo_dry::exec()
{
  if(grid->swspatialorder== "2")
    calcbuoyancytend_2nd(fields->wt->data, fields->s["th"]->data, threfh);
  else if(grid->swspatialorder == "4")
    calcbuoyancytend_4th(fields->wt->data, fields->s["th"]->data, threfh);

  return 0;
}

int cthermo_dry::execstats(mask *m)
{
  // calculate the buoyancy and its surface flux for the profiles
  calcbuoyancy(fields->s["tmp1"]->data, fields->s["th"]->data, thref);
  calcbuoyancyfluxbot(fields->s["tmp1"]->datafluxbot, fields->s["th"]->datafluxbot, threfh);

  // define the location
  const int sloc[] = {0,0,0};

  // calculate the mean
  stats->calcmean(fields->s["tmp1"]->data, m->profs["b"].data, NO_OFFSET, sloc,
                  fields->s["tmp3"]->data, stats->nmask);

  // calculate the moments
  for(int n=2; n<5; ++n)
  {
    std::stringstream ss;
    ss << n;
    std::string sn = ss.str();
    stats->calcmoment(fields->s["tmp1"]->data, m->profs["b"].data, m->profs["b"+sn].data, n, sloc,
                      fields->s["tmp3"]->data, stats->nmask);
  }

  // calculate the gradients
  if(grid->swspatialorder == "2")
    stats->calcgrad_2nd(fields->s["tmp1"]->data, m->profs["bgrad"].data, grid->dzhi, sloc,
                        fields->s["tmp4"]->data, stats->nmaskh);
  if(grid->swspatialorder == "4")
    stats->calcgrad_4th(fields->s["tmp1"]->data, m->profs["bgrad"].data, grid->dzhi4, sloc,
                        fields->s["tmp4"]->data, stats->nmaskh);

  // calculate turbulent fluxes
  if(grid->swspatialorder == "2")
    stats->calcflux_2nd(fields->s["tmp1"]->data, m->profs["b"].data, fields->w->data, m->profs["w"].data,
                        m->profs["bw"].data, fields->s["tmp2"]->data, sloc,
                        fields->s["tmp4"]->data, stats->nmaskh);
  if(grid->swspatialorder == "4")
    stats->calcflux_4th(fields->s["tmp1"]->data, fields->w->data, m->profs["bw"].data, fields->s["tmp2"]->data, sloc,
                        fields->s["tmp4"]->data, stats->nmaskh);

  // calculate diffusive fluxes
  if(model->diff->getname() == "les2s")
  {
    cdiff_les2s *diffptr = static_cast<cdiff_les2s *>(model->diff);
    stats->calcdiff_2nd(fields->s["tmp1"]->data, fields->w->data, fields->s["evisc"]->data,
                        m->profs["bdiff"].data, grid->dzhi,
                        fields->s["tmp1"]->datafluxbot, fields->s["tmp1"]->datafluxtop, diffptr->tPr, sloc,
                        fields->s["tmp4"]->data, stats->nmaskh);
  }
  else
    stats->calcdiff_4th(fields->s["tmp1"]->data, m->profs["bdiff"].data, grid->dzhi4, fields->s["th"]->visc, sloc,
                        fields->s["tmp4"]->data, stats->nmaskh);

  // calculate the total fluxes
  stats->addfluxes(m->profs["bflux"].data, m->profs["bw"].data, m->profs["bdiff"].data);

  // calculate the sorted buoyancy profile
  stats->calcsortprof(fields->sd["tmp1"]->data, fields->sd["tmp2"]->data, m->profs["bsort"].data);

  return 0;
}

int cthermo_dry::execcross()
{
  int nerror = 0;

  // With one additional temp field, we wouldn't have to re-calculate the ql or b field for simple,lngrad,path, etc.
  for(std::vector<std::string>::iterator it=crosslist.begin(); it<crosslist.end(); ++it)
  {
    if(*it == "b")
    {
      getthermofield(fields->s["tmp1"], fields->s["tmp2"], *it);
      nerror += model->cross->crosssimple(fields->s["tmp1"]->data, fields->s["tmp2"]->data, *it);
    }
    else if(*it == "blngrad")
    {
      getthermofield(fields->s["tmp1"], fields->s["tmp2"], "b");
      // Note: tmp1 twice used as argument -> overwritten in crosspath()
      nerror += model->cross->crosslngrad(fields->s["tmp1"]->data, fields->s["tmp2"]->data, fields->s["tmp1"]->data, grid->dzi4, *it);
    }
    else if(*it == "bbot" or *it == "bfluxbot")
    {
      getbuoyancysurf(fields->s["tmp1"]);
      if(*it == "bbot")
        nerror += model->cross->crossplane(fields->s["tmp1"]->databot, fields->s["tmp1"]->data, "bbot");
      else if(*it == "bfluxbot")
        nerror += model->cross->crossplane(fields->s["tmp1"]->datafluxbot, fields->s["tmp1"]->data, "bfluxbot");
    }
  }  

  return nerror; 
}

int cthermo_dry::checkthermofield(std::string name)
{
  if(name == "b")
    return 0;
  else
    return 1;
}

int cthermo_dry::getthermofield(cfield3d *fld, cfield3d *tmp, std::string name)
{
  if(name == "b")
    calcbuoyancy(fld->data, fields->s["th"]->data, thref);
  else if(name == "N2")
    calcN2(fld->data, fields->s["th"]->data, grid->dzi, thref);
  else
    return 1;

  return 0;
}

int cthermo_dry::getbuoyancyfluxbot(cfield3d *bfield)
{
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->s["th"]->datafluxbot, threfh);

  return 0;
}

int cthermo_dry::getbuoyancysurf(cfield3d *bfield)
{
  calcbuoyancybot(bfield->data, bfield->databot,
                  fields->s["th"]->data, fields->s["th"]->databot, thref, threfh);
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->s["th"]->datafluxbot, threfh);
  return 0;
}

int cthermo_dry::getprogvars(std::vector<std::string> *list)
{
  list->push_back("th");
  return 0;
}

int cthermo_dry::calcbuoyancy(double * restrict b, double * restrict th, double * restrict thref)
{
  int ijk,jj,kk;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int k=0; k<grid->kcells; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        b[ijk] = grav/thref[k] * (th[ijk] - thref[k]);
      }

  return 0;
}

int cthermo_dry::calcN2(double * restrict N2, double * restrict th, double * restrict dzi, double * restrict thref)
{
  int ijk,jj,kk;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int k=0; k<grid->kcells; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        N2[ijk] = grav/thref[k]*0.5*(th[ijk+kk] - th[ijk-kk])*dzi[k];
      }

  return 0;
}

int cthermo_dry::calcbuoyancybot(double * restrict b , double * restrict bbot,
                                 double * restrict th, double * restrict thbot,
                                 double * restrict thref, double * restrict threfh)
{
  int ij,ijk,jj,kk,kstart;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  kstart = grid->kstart;

  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      bbot[ij] = grav/threfh[kstart] * (thbot[ij] - threfh[kstart]);
      b[ijk]   = grav/thref [kstart] * (th[ijk]   - thref [kstart]);
    }

  return 0;
}

int cthermo_dry::calcbuoyancyfluxbot(double * restrict bfluxbot, double * restrict thfluxbot, double * restrict threfh)
{
  int ij,jj;
  jj = grid->icells;

  int kstart = grid->kstart;

  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ij  = i + j*jj;
      bfluxbot[ij] = grav/threfh[kstart]*thfluxbot[ij];
    }

  return 0;
}

int cthermo_dry::calcbuoyancytend_2nd(double * restrict wt, double * restrict th, double * restrict threfh)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int k=grid->kstart+1; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        wt[ijk] += grav/threfh[k] * (interp2(th[ijk-kk], th[ijk]) - threfh[k]);
      }

  return 0;
}

int cthermo_dry::calcbuoyancytend_4th(double * restrict wt, double * restrict th, double * restrict threfh)
{
  int ijk,jj;
  int kk1,kk2;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  for(int k=grid->kstart+1; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk1;
        wt[ijk] += grav/threfh[k] * (interp4(th[ijk-kk2], th[ijk-kk1], th[ijk], th[ijk+kk1]) - threfh[k]);
      }

  return 0;
}

inline double cthermo_dry::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}

inline double cthermo_dry::interp4(const double a, const double b, const double c, const double d)
{
  return (-a + 9.*b + 9.*c - d) / 16.;
}

