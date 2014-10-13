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
#include <algorithm>
#include <cmath>
#include "grid.h"
#include "fields.h"
#include "thermo_dry.h"
#include "defines.h"
#include "constants.h"
#include "fd.h"
#include "model.h"
#include "stats.h"
#include "diff_smag2.h"
#include "master.h"
#include "cross.h"
#include "dump.h"

using fd::o2::interp2;
using fd::o4::interp4;
using namespace constants;

ThermoDry::ThermoDry(Model *modelin, Input *inputin) : Thermo(modelin, inputin)
{
  swthermo = "dry";

  thref  = 0;
  threfh = 0;
  pref   = 0;
  prefh  = 0;
  exner  = 0;
  exnerh = 0;

  fields->initPrognosticField("th", "Potential Temperature", "K");

  int nerror = 0;
  nerror += inputin->getItem(&fields->sp["th"]->visc, "fields", "svisc", "th");
  nerror += inputin->getList(&crosslist , "thermo", "crosslist" , "");
  nerror += inputin->getList(&dumplist ,  "thermo", "dumplist" ,  "");

  if(nerror)
    throw 1;
}

ThermoDry::~ThermoDry()
{
  delete[] this->thref;
  delete[] this->threfh;
  delete[] this->pref;
  delete[] this->prefh;
  delete[] this->exner;
  delete[] this->exnerh;

  #ifdef USECUDA
  clearDevice();
  #endif
}

void ThermoDry::init()
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
}

void ThermoDry::create(Input *inputin)
{
  // Only in case of Boussinesq, read in reference potential temperature
  if(model->swbasestate == "boussinesq")
  {
    if(inputin->getItem(&thref0, "thermo", "thref0", ""))
      throw 1;
  }

  // For dry thermo, only anelastic needs surface pressure
  if(model->swbasestate == "anelastic")
  {
    if(inputin->getItem(&pbot, "thermo", "pbot", ""))
      throw 1;
  }

  // Setup base state for anelastic solver
  if(model->swbasestate == "anelastic")
  {
    // take the initial profile as the reference
    if(inputin->getProf(&thref[grid->kstart], "th", grid->kmax))
      throw 1;

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
  if(stats->getSwitch() == "1")
  {
    // Add base state profiles to statistics -> needed/wanted for Boussinesq? Or write as 0D var?
    //stats->addfixedprof("pref",    "Full level basic state pressure", "Pa",     "z",  pref);
    //stats->addfixedprof("prefh",   "Half level basic state pressure", "Pa",     "zh", prefh);
    stats->addFixedProf("rhoref",  "Full level basic state density",  "kg m-3", "z",  fields->rhoref);
    stats->addFixedProf("rhorefh", "Half level basic state density",  "kg m-3", "zh", fields->rhorefh);
    stats->addFixedProf("thref",   "Full level reference potential temperature", "K", "z",thref);
    stats->addFixedProf("threfh",  "Half level reference potential temperature", "K", "zh",thref);

    stats->addProf("b", "Buoyancy", "m s-2", "z");
    for(int n=2; n<5; ++n)
    {
      std::stringstream ss;
      ss << n;
      std::string sn = ss.str();
      stats->addProf("b"+sn, "Moment " +sn+" of the buoyancy", "(m s-2)"+sn,"z");
    }

    stats->addProf("bgrad", "Gradient of the buoyancy", "s-2", "zh");
    stats->addProf("bw"   , "Turbulent flux of the buoyancy", "m2 s-3", "zh");
    stats->addProf("bdiff", "usive flux of the buoyancy", "m2 s-3", "zh");
    stats->addProf("bflux", "Total flux of the buoyancy", "m2 s-3", "zh");

    stats->addProf("bsort", "Sorted buoyancy", "m s-2", "z");
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

  // Check if fields in dumplist are retrievable thermo fields, if not delete them and print warning
  std::vector<std::string>::iterator dumpvar=dumplist.begin();
  while(dumpvar != dumplist.end())
  {
    if(checkThermoField(*dumpvar))
    {
      master->printWarning("field %s in [thermo][dumplist] is not a thermo field\n", dumpvar->c_str());
      dumpvar = dumplist.erase(dumpvar);  // erase() returns iterator of next element
    }
    else
      ++dumpvar;
  }
}

#ifndef USECUDA
void ThermoDry::exec()
{
  if(grid->swspatialorder== "2")
    calcbuoyancytend_2nd(fields->wt->data, fields->sp["th"]->data, threfh);
  else if(grid->swspatialorder == "4")
    calcbuoyancytend_4th(fields->wt->data, fields->sp["th"]->data, threfh);
}
#endif

void ThermoDry::execStats(Mask *m)
{
  const double NoOffset = 0.;

  // calculate the buoyancy and its surface flux for the profiles
  calcbuoyancy(fields->atmp["tmp1"]->data, fields->sp["th"]->data, thref);
  calcbuoyancyfluxbot(fields->atmp["tmp1"]->datafluxbot, fields->sp["th"]->datafluxbot, threfh);

  // define the location
  const int sloc[] = {0,0,0};

  // calculate the mean
  stats->calcMean(m->profs["b"].data, fields->atmp["tmp1"]->data, NoOffset, sloc,
                  fields->atmp["tmp3"]->data, stats->nmask);

  // calculate the moments
  for(int n=2; n<5; ++n)
  {
    std::stringstream ss;
    ss << n;
    std::string sn = ss.str();
    stats->calcMoment(fields->atmp["tmp1"]->data, m->profs["b"].data, m->profs["b"+sn].data, n, sloc,
                      fields->atmp["tmp3"]->data, stats->nmask);
  }

  // calculate the gradients
  if(grid->swspatialorder == "2")
    stats->calcGrad_2nd(fields->atmp["tmp1"]->data, m->profs["bgrad"].data, grid->dzhi, sloc,
                        fields->atmp["tmp4"]->data, stats->nmaskh);
  else if(grid->swspatialorder == "4")
    stats->calcGrad_4th(fields->atmp["tmp1"]->data, m->profs["bgrad"].data, grid->dzhi4, sloc,
                        fields->atmp["tmp4"]->data, stats->nmaskh);

  // calculate turbulent fluxes
  if(grid->swspatialorder == "2")
    stats->calcFlux_2nd(fields->atmp["tmp1"]->data, m->profs["b"].data, fields->w->data, m->profs["w"].data,
                        m->profs["bw"].data, fields->atmp["tmp2"]->data, sloc,
                        fields->atmp["tmp4"]->data, stats->nmaskh);
  else if(grid->swspatialorder == "4")
    stats->calcFlux_4th(fields->atmp["tmp1"]->data, fields->w->data, m->profs["bw"].data, fields->atmp["tmp2"]->data, sloc,
                        fields->atmp["tmp4"]->data, stats->nmaskh);

  // calculate diffusive fluxes
  if(grid->swspatialorder == "2")
  {
    if(model->diff->getName() == "smag2")
    {
      DiffSmag2 *diffptr = static_cast<DiffSmag2 *>(model->diff);
      stats->calcDiff_2nd(fields->atmp["tmp1"]->data, fields->w->data, fields->sd["evisc"]->data,
                          m->profs["bdiff"].data, grid->dzhi,
                          fields->atmp["tmp1"]->datafluxbot, fields->atmp["tmp1"]->datafluxtop, diffptr->tPr, sloc,
                          fields->atmp["tmp4"]->data, stats->nmaskh);
    }
    else
      stats->calcDiff_2nd(fields->atmp["tmp1"]->data, m->profs["bdiff"].data, grid->dzhi, fields->sp["th"]->visc, sloc,
                          fields->atmp["tmp4"]->data, stats->nmaskh);
  }
  else if(grid->swspatialorder == "4")
  {
    stats->calcDiff_4th(fields->atmp["tmp1"]->data, m->profs["bdiff"].data, grid->dzhi4, fields->sp["th"]->visc, sloc,
                        fields->atmp["tmp4"]->data, stats->nmaskh);
  }

  // calculate the total fluxes
  stats->addFluxes(m->profs["bflux"].data, m->profs["bw"].data, m->profs["bdiff"].data);

  // calculate the sorted buoyancy profile
  //stats->calcSortedProf(fields->sd["tmp1"]->data, fields->sd["tmp2"]->data, m->profs["bsort"].data);
}

void ThermoDry::execCross()
{
  int nerror = 0;

  Cross *cross = model->cross;

  // With one additional temp field, we wouldn't have to re-calculate the ql or b field for simple,lngrad,path, etc.
  for(std::vector<std::string>::iterator it=crosslist.begin(); it<crosslist.end(); ++it)
  {
    /* BvS: for now, don't call getThermoField() or getBuoyancySurf(), but directly the function itself. With CUDA enabled, 
       statistics etc. is done on the host, while getThermoField() is executed on the GPU */ 
    
    if(*it == "b")
    {
      //getThermoField(fields->s["tmp1"], fields->s["tmp2"], *it);
      calcbuoyancy(fields->atmp["tmp1"]->data, fields->sp["th"]->data, thref);
      nerror += cross->crossSimple(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, *it);
    }
    else if(*it == "blngrad")
    {
      //getThermoField(fields->s["tmp1"], fields->s["tmp2"], "b");
      calcbuoyancy(fields->atmp["tmp1"]->data, fields->sp["th"]->data, thref);
      // Note: tmp1 twice used as argument -> overwritten in crosspath()
      nerror += cross->crossLngrad(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, grid->dzi4, *it);
    }
    else if(*it == "bbot" or *it == "bfluxbot")
    {
      //getBuoyancySurf(fields->s["tmp1"]);
      calcbuoyancybot(fields->atmp["tmp1"]->data, fields->atmp["tmp1"]->databot, fields->sp["th"]->data, fields->sp["th"]->databot, thref, threfh);
      calcbuoyancyfluxbot(fields->atmp["tmp1"]->datafluxbot, fields->sp["th"]->datafluxbot, threfh);

      if(*it == "bbot")
        nerror += cross->crossPlane(fields->atmp["tmp1"]->databot, fields->atmp["tmp1"]->data, "bbot");
      else if(*it == "bfluxbot")
        nerror += cross->crossPlane(fields->atmp["tmp1"]->datafluxbot, fields->atmp["tmp1"]->data, "bfluxbot");
    }
  }

  if(nerror)
    throw 1;
}

void ThermoDry::execDump()
{
  for(std::vector<std::string>::const_iterator it=dumplist.begin(); it<dumplist.end(); ++it)
  {
    // TODO BvS restore getThermoField(), the combination of checkThermoField with getThermoField is more elegant... 
    if(*it == "b")
      calcbuoyancy(fields->atmp["tmp2"]->data, fields->sp["th"]->data, thref);
    else
      throw 1;

    model->dump->saveDump(fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, *it);
  }
}

bool ThermoDry::checkThermoField(std::string name)
{
  if(name == "b")
    return false;
  else
    return true;
}

#ifndef USECUDA
void ThermoDry::getThermoField(Field3d *fld, Field3d *tmp, std::string name)
{
  if(name == "b")
    calcbuoyancy(fld->data, fields->sp["th"]->data, thref);
  else if(name == "N2")
    calcN2(fld->data, fields->sp["th"]->data, grid->dzi, thref);
  else
    throw 1;
}
#endif

#ifndef USECUDA
void ThermoDry::getBuoyancyFluxbot(Field3d *bfield)
{
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->sp["th"]->datafluxbot, threfh);
}
#endif

#ifndef USECUDA
void ThermoDry::getBuoyancySurf(Field3d *bfield)
{
  calcbuoyancybot(bfield->data, bfield->databot,
                  fields->sp["th"]->data, fields->sp["th"]->databot, thref, threfh);
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->sp["th"]->datafluxbot, threfh);
}
#endif

void ThermoDry::getProgVars(std::vector<std::string> *list)
{
  list->push_back("th");
}

int ThermoDry::calcbuoyancy(double * restrict b, double * restrict th, double * restrict thref)
{
  int ijk,jj,kk;
  jj = grid->icells;
  kk = grid->ijcells;

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

int ThermoDry::calcN2(double * restrict N2, double * restrict th, double * restrict dzi, double * restrict thref)
{
  int ijk,jj,kk;
  jj = grid->icells;
  kk = grid->ijcells;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        N2[ijk] = grav/thref[k]*0.5*(th[ijk+kk] - th[ijk-kk])*dzi[k];
      }

  return 0;
}

int ThermoDry::calcbuoyancybot(double * restrict b , double * restrict bbot,
                                 double * restrict th, double * restrict thbot,
                                 double * restrict thref, double * restrict threfh)
{
  int ij,ijk,jj,kk,kstart;
  jj = grid->icells;
  kk = grid->ijcells;
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

int ThermoDry::calcbuoyancyfluxbot(double * restrict bfluxbot, double * restrict thfluxbot, double * restrict threfh)
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

int ThermoDry::calcbuoyancytend_2nd(double * restrict wt, double * restrict th, double * restrict threfh)
{
  using namespace fd::o2;

  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->ijcells;

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

int ThermoDry::calcbuoyancytend_4th(double * restrict wt, double * restrict th, double * restrict threfh)
{
  int ijk,jj;
  int kk1,kk2;

  jj  = grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

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
