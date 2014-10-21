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
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "thermo_moist.h"
#include "diff_smag2.h"
#include "defines.h"
#include "constants.h"
#include "fd.h"
#include "model.h"
#include "stats.h"
#include "master.h"
#include "cross.h"
#include "dump.h"

using fd::o2::interp2;
using fd::o4::interp4;
using namespace constants;

ThermoMoist::ThermoMoist(Model *modelin, Input *inputin) : Thermo(modelin, inputin)
{
  swthermo = "moist";

  thl0 = 0;
  qt0  = 0;

  thvref  = 0;
  thvrefh = 0;
  exnref  = 0;
  exnrefh = 0;
  pref    = 0;
  prefh   = 0;

  thvref_g  = 0;
  thvrefh_g = 0;
  exnref_g  = 0;
  exnrefh_g = 0;
  pref_g    = 0;
  prefh_g   = 0;

  // Initialize the prognostic fields
  fields->initPrognosticField("thl", "Liquid water potential temperature", "K");
  fields->initPrognosticField("qt", "Total water mixing ratio", "kg kg-1");

  int nerror = 0;
  nerror += inputin->getItem(&fields->sp["thl" ]->visc, "fields", "svisc", "thl" );
  nerror += inputin->getItem(&fields->sp["qt"]->visc, "fields", "svisc", "qt");
  nerror += inputin->getItem(&pbot, "thermo", "pbot", "");

  // Read list of cross sections
  nerror += inputin->getList(&crosslist, "thermo", "crosslist", "");
  // Read list of 3d dumps
  nerror += inputin->getList(&dumplist,  "thermo", "dumplist",  "");
  
  // BvS test for updating hydrostatic prssure during run
  // swupdate..=0 -> initial base state pressure used in saturation calculation
  // swupdate..=1 -> base state pressure updated before saturation calculation
  nerror += inputin->getItem(&swupdatebasestate, "thermo", "swupdatebasestate", ""); 

  if(nerror)
    throw 1;
}

ThermoMoist::~ThermoMoist()
{
  delete[] thl0;
  delete[] qt0;
  delete[] thvref;
  delete[] thvrefh;
  delete[] exnref;
  delete[] exnrefh;
  delete[] pref;
  delete[] prefh;

  #ifdef USECUDA
  clearDevice();
  #endif
}

void ThermoMoist::init()
{
  stats = model->stats;

  thl0    = new double[grid->kcells];
  qt0     = new double[grid->kcells];
  thvref  = new double[grid->kcells];
  thvrefh = new double[grid->kcells];
  exnref  = new double[grid->kcells];
  exnrefh = new double[grid->kcells];
  pref    = new double[grid->kcells];
  prefh   = new double[grid->kcells];

  for(int k=0; k<grid->kcells; ++k)
  {
    thl0   [k] = 0.;
    qt0    [k] = 0.;
    thvref [k] = 0.;
    thvrefh[k] = 0.;
    exnref [k] = 0.;
    exnrefh[k] = 0.;
    pref   [k] = 0.;
    prefh  [k] = 0.;
  }
}

void ThermoMoist::create(Input *inputin)
{
  int kstart = grid->kstart;
  int kend   = grid->kend;

  // Enable automated calculation of horizontally averaged fields
  if(swupdatebasestate)
    fields->set_calcMeanProfs(true);

  // Calculate the base state profiles. With swupdatebasestate=1, these profiles are updated on every iteration. 
  // 1. Take the initial profile as the reference
  if(inputin->getProf(&thl0[grid->kstart], "thl", grid->kmax))
    throw 1;
  if(inputin->getProf(&qt0[grid->kstart], "qt", grid->kmax))
    throw 1;

  // 2. Calculate surface and model top values thl and qt
  double thl0s, qt0s, thl0t, qt0t;
  thl0s = thl0[kstart] - grid->z[kstart]*(thl0[kstart+1]-thl0[kstart])*grid->dzhi[kstart+1];
  qt0s  = qt0[kstart]  - grid->z[kstart]*(qt0[kstart+1] -qt0[kstart] )*grid->dzhi[kstart+1];
  thl0t = thl0[kend-1] + (grid->zh[kend]-grid->z[kend-1])*(thl0[kend-1]-thl0[kend-2])*grid->dzhi[kend-1];
  qt0t  = qt0[kend-1]  + (grid->zh[kend]-grid->z[kend-1])*(qt0[kend-1]- qt0[kend-2] )*grid->dzhi[kend-1];

  // 3. Set the ghost cells for the reference temperature and moisture
  thl0[kstart-1]  = 2.*thl0s - thl0[kstart];
  thl0[kend]      = 2.*thl0t - thl0[kend-1];
  qt0[kstart-1]   = 2.*qt0s  - qt0[kstart];
  qt0[kend]       = 2.*qt0t  - qt0[kend-1];

  // 4. Calculate the initial/reference base state
  calcBaseState(pref, prefh, fields->rhoref, fields->rhorefh, thvref, thvrefh, exnref, exnrefh, thl0, qt0);

  // 5. In Boussinesq mode, overwrite reference temperature and density
  if(model->swbasestate == "boussinesq")
  {
    if(inputin->getItem(&thvref0, "thermo", "thvref0", ""))
      throw 1;

    for(int k=0; k<grid->kcells; ++k)
    {
      fields->rhoref[k]  = 1.;
      fields->rhorefh[k] = 1.;
      thvref[k]          = thvref0;
      thvrefh[k]         = thvref0;
    }
  }

  initStat();
  initCross();
  initDump();
}

#ifndef USECUDA
void ThermoMoist::exec()
{
  const int kk = grid->ijcells;
  const int kcells = grid->kcells;

  // Re-calculate hydrostatic pressure and exner, pass dummy as rhoref,thvref to prevent overwriting base state 
  double *tmp2 = fields->atmp["tmp2"]->data;
  if(swupdatebasestate)
    calcBaseState(pref, prefh, &tmp2[0*kcells], &tmp2[1*kcells], &tmp2[2*kcells], &tmp2[3*kcells], exnref, exnrefh, 
                  fields->sp["thl"]->datamean, fields->sp["qt"]->datamean);
  
  // extend later for gravity vector not normal to surface
  if(grid->swspatialorder == "2")
  {
    calcBuoyancyTend_2nd(fields->wt->data, fields->sp["thl"]->data, fields->sp["qt"]->data, prefh,
                         &fields->atmp["tmp2"]->data[0*kk], &fields->atmp["tmp2"]->data[1*kk], &fields->atmp["tmp2"]->data[2*kk],
                         thvrefh);
  }
  else if(grid->swspatialorder == "4")
  {
    calcBuoyancyTend_4th(fields->wt->data, fields->sp["thl"]->data, fields->sp["qt"]->data, prefh,
                         &fields->atmp["tmp2"]->data[0*kk], &fields->atmp["tmp2"]->data[1*kk], &fields->atmp["tmp2"]->data[2*kk],
                         thvrefh);
  }
}
#endif

void ThermoMoist::getMask(Field3d *mfield, Field3d *mfieldh, Mask *m)
{
  if(m->name == "ql")
  {
    calcLiquidWater(fields->atmp["tmp1"]->data, fields->sp["thl"]->data, fields->sp["qt"]->data, pref);
    calcMask_ql(mfield->data, mfieldh->data, mfieldh->databot,
                stats->nmask, stats->nmaskh, &stats->nmaskbot,
                fields->atmp["tmp1"]->data);
  }
  else if(m->name == "qlcore")
  {
    calcBuoyancy(fields->atmp["tmp2"]->data, fields->sp["thl"]->data, fields->sp["qt"]->data, pref, fields->atmp["tmp1"]->data,thvref);
    // calculate the mean buoyancy to determine positive buoyancy
    grid->calcMean(fields->atmp["tmp2"]->datamean, fields->atmp["tmp2"]->data, grid->kcells);
    calcLiquidWater(fields->atmp["tmp1"]->data, fields->sp["thl"]->data, fields->sp["qt"]->data, pref);
    calcMask_qlcore(mfield->data, mfieldh->data, mfieldh->databot,
                    stats->nmask, stats->nmaskh, &stats->nmaskbot,
                    fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp2"]->datamean);
  }
}

void ThermoMoist::calcMask_ql(double * restrict mask, double * restrict maskh, double * restrict maskbot,
                              int * restrict nmask, int * restrict nmaskh, int * restrict nmaskbot,
                              double * restrict ql)
{
  int ijk,ij,jj,kk;
  int kstart;

  jj = grid->icells;
  kk = grid->ijcells;
  kstart = grid->kstart;

  int ntmp;

  for(int k=grid->kstart; k<grid->kend; k++)
  {
    nmask[k] = 0;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ntmp = ql[ijk] > 0.;
        nmask[k] += ntmp;
        mask[ijk] = (double)ntmp;
      }
  }

  for(int k=grid->kstart; k<grid->kend+1; k++)
  {
    nmaskh[k] = 0;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ntmp = (ql[ijk-kk] + ql[ijk]) > 0.;

        nmaskh[k] += ntmp;
        maskh[ijk] = (double)ntmp;
      }
  }

  // Set the mask for surface projected quantities
  // In this case: ql at surface
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      maskbot[ij] = maskh[ijk];
    }

  grid->boundaryCyclic(mask);
  grid->boundaryCyclic(maskh);
  grid->boundaryCyclic2d(maskbot);

  master->sum(nmask , grid->kcells);
  master->sum(nmaskh, grid->kcells);
  *nmaskbot = nmaskh[grid->kstart];

  // BvS: should no longer be necessary now that the ql ghost cells are set to zero
  //nmaskh[kstart] = 0;
  //nmaskh[kend  ] = 0;
}

void ThermoMoist::calcMask_qlcore(double * restrict mask, double * restrict maskh, double * restrict maskbot,
                                  int * restrict nmask, int * restrict nmaskh, int * restrict nmaskbot,
                                  double * restrict ql, double * restrict b, double * restrict bmean)
{
  int ijk,ij,jj,kk;
  int kstart;

  jj = grid->icells;
  kk = grid->ijcells;
  kstart = grid->kstart;

  int ntmp;

  for(int k=grid->kstart; k<grid->kend; k++)
  {
    nmask[k] = 0;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ntmp = (ql[ijk] > 0.)*(b[ijk]-bmean[k] > 0.);
        nmask[k] += ntmp;
        mask[ijk] = (double)ntmp;
      }
  }

  for(int k=grid->kstart; k<grid->kend+1; k++)
  {
    nmaskh[k] = 0;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ntmp = (ql[ijk-kk]+ql[ijk] > 0.)*(b[ijk-kk]+b[ijk]-bmean[k-1]-bmean[k] > 0.);
        nmaskh[k] += ntmp;
        maskh[ijk] = (double)ntmp;
      }
  }

  // Set the mask for surface projected quantities
  // In this case: qlcore at surface
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      maskbot[ij] = maskh[ijk];
    }

  grid->boundaryCyclic(mask);
  grid->boundaryCyclic(maskh);
  grid->boundaryCyclic2d(maskbot);

  master->sum(nmask , grid->kcells);
  master->sum(nmaskh, grid->kcells);
  *nmaskbot = nmaskh[grid->kstart];

  // BvS: should no longer be necessary now that the ql ghost cells are set to zero
  //nmaskh[kstart] = 0;
  //nmaskh[kend  ] = 0;
}

void ThermoMoist::execStats(Mask *m)
{
  const double NoOffset = 0.;

  // calc the buoyancy and its surface flux for the profiles
  calcBuoyancy(fields->atmp["tmp1"]->data, fields->sp["thl"]->data, fields->sp["qt"]->data, pref, fields->atmp["tmp2"]->data, thvref);
  calcBuoyancyFluxBot(fields->atmp["tmp1"]->datafluxbot, fields->sp["thl"]->databot, fields->sp["thl"]->datafluxbot, fields->sp["qt"]->databot, fields->sp["qt"]->datafluxbot, thvrefh);

  // define location
  const int sloc[] = {0,0,0};

  // mean
  stats->calcMean(m->profs["b"].data, fields->atmp["tmp1"]->data, NoOffset, sloc,
                  fields->atmp["tmp3"]->data, stats->nmask);

  // moments
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
    {
      stats->calcDiff_2nd(fields->atmp["tmp1"]->data, m->profs["bdiff"].data, grid->dzhi, fields->sp["th"]->visc, sloc,
                          fields->atmp["tmp4"]->data, stats->nmaskh);
    }
  }
  else if(grid->swspatialorder == "4")
  {
    // take the diffusivity of temperature for that of buoyancy
    stats->calcDiff_4th(fields->atmp["tmp1"]->data, m->profs["bdiff"].data, grid->dzhi4, fields->sp["th"]->visc, sloc,
                        fields->atmp["tmp4"]->data, stats->nmaskh);
  }

  // calculate the total fluxes
  stats->addFluxes(m->profs["bflux"].data, m->profs["bw"].data, m->profs["bdiff"].data);

  // calculate the liquid water stats
  calcLiquidWater(fields->atmp["tmp1"]->data, fields->sp["thl"]->data, fields->sp["qt"]->data, pref);
  stats->calcMean(m->profs["ql"].data, fields->atmp["tmp1"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);
  stats->calcCount(fields->atmp["tmp1"]->data, m->profs["cfrac"].data, 0.,
                   fields->atmp["tmp3"]->data, stats->nmask);

  stats->calcCover(fields->atmp["tmp1"]->data, fields->atmp["tmp4"]->databot, &stats->nmaskbot, &m->tseries["ccover"].data, 0.);
  stats->calcPath (fields->atmp["tmp1"]->data, fields->atmp["tmp4"]->databot, &stats->nmaskbot, &m->tseries["lwp"].data);

  // Calculate base state in tmp array
  if(swupdatebasestate == 1)
  {
    const int kcells = grid->kcells;
    double * restrict tmp2 = fields->atmp["tmp2"]->data;
    calcBaseState(&tmp2[0*kcells], &tmp2[1*kcells], &tmp2[2*kcells], &tmp2[3*kcells], 
                  &tmp2[4*kcells], &tmp2[5*kcells], &tmp2[6*kcells], &tmp2[7*kcells], 
                  fields->sp["thl"]->datamean, fields->sp["qt"]->datamean);

    for(int k=0; k<kcells; ++k)
    {
      m->profs["ph"  ].data[k] = tmp2[0*kcells+k];
      m->profs["phh" ].data[k] = tmp2[1*kcells+k];
      m->profs["rho" ].data[k] = tmp2[2*kcells+k];
      m->profs["rhoh"].data[k] = tmp2[3*kcells+k];
    } 
  }
}

void ThermoMoist::execCross()
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
      calcBuoyancy(fields->atmp["tmp1"]->data, fields->sp["thl"]->data, fields->sp["qt"]->data, pref, fields->atmp["tmp2"]->data, thvref);
      nerror += cross->crossSimple(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, *it);
    }
    else if(*it == "ql")
    {
      calcLiquidWater(fields->atmp["tmp1"]->data, fields->sp["thl"]->data, fields->sp["qt"]->data, pref);
      nerror += cross->crossSimple(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, *it);
    }
    else if(*it == "blngrad")
    {
      calcBuoyancy(fields->atmp["tmp1"]->data, fields->sp["thl"]->data, fields->sp["qt"]->data, pref, fields->atmp["tmp2"]->data, thvref);
      // Note: tmp1 twice used as argument -> overwritten in crosspath()
      nerror += cross->crossLngrad(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, grid->dzi4, *it);
    }
    else if(*it == "qlpath")
    {
      calcLiquidWater(fields->atmp["tmp1"]->data, fields->sp["thl"]->data, fields->sp["qt"]->data, pref);
      // Note: tmp1 twice used as argument -> overwritten in crosspath()
      nerror += cross->crossPath(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, "qlpath");
    }
    else if(*it == "bbot" or *it == "bfluxbot")
    {
      calcBuoyancyBot(fields->atmp["tmp1"]->data, fields->atmp["tmp1"]->databot, fields->sp["thl" ]->data, fields->sp["thl"]->databot, fields->sp["qt"]->data, fields->sp["qt"]->databot, thvref, thvrefh);
      calcBuoyancyFluxBot(fields->atmp["tmp1"]->datafluxbot, fields->sp["thl"]->databot, fields->sp["thl"]->datafluxbot, fields->sp["qt"]->databot, fields->sp["qt"]->datafluxbot, thvrefh);

      if(*it == "bbot")
        nerror += cross->crossPlane(fields->atmp["tmp1"]->databot, fields->atmp["tmp1"]->data, "bbot");
      else if(*it == "bfluxbot")
        nerror += cross->crossPlane(fields->atmp["tmp1"]->datafluxbot, fields->atmp["tmp1"]->data, "bfluxbot");
    }
  }

  if(nerror)
    throw 1;
}

void ThermoMoist::execDump()
{
  for(std::vector<std::string>::const_iterator it=dumplist.begin(); it<dumplist.end(); ++it)
  {
    // TODO BvS restore getThermoField(), the combination of checkThermoField with getThermoField is more elegant... 
    if(*it == "b")
      calcBuoyancy(fields->atmp["tmp2"]->data, fields->sp["thl"]->data, fields->sp["qt"]->data, pref, fields->atmp["tmp1"]->data, thvref);
    else if(*it == "ql")
      calcLiquidWater(fields->atmp["tmp2"]->data, fields->sp["thl"]->data, fields->sp["qt"]->data, pref);
    else
      throw 1;

    model->dump->saveDump(fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, *it);
  }
}

bool ThermoMoist::checkThermoField(std::string name)
{
  if(name == "b" || name == "ql")
    return false;
  else
    return true;
}

#ifndef USECUDA
void ThermoMoist::getThermoField(Field3d *fld, Field3d *tmp, std::string name)
{
  const int kcells = grid->kcells;

  // BvS: getThermoField() is called from subgrid-model, before thermo(), so re-calculate the hydrostatic pressure
  // Pass dummy as rhoref,thvref to prevent overwriting base state 
  double * restrict tmp2 = fields->atmp["tmp2"]->data;
  if(swupdatebasestate)
    calcBaseState(pref, prefh, &tmp2[0*kcells], &tmp2[1*kcells], &tmp2[2*kcells], &tmp2[3*kcells], exnref, exnrefh, 
                  fields->sp["thl"]->datamean, fields->sp["qt"]->datamean);

  if(name == "b")
    calcBuoyancy(fld->data, fields->sp["thl"]->data, fields->sp["qt"]->data, pref, tmp->data, thvref);
  else if(name == "ql")
    calcLiquidWater(fld->data, fields->sp["thl"]->data, fields->sp["qt"]->data, pref);
  else if(name == "N2")
    calcN2(fld->data, fields->sp["thl"]->data, grid->dzi, thvref);
  else
    throw 1;
}
#endif

#ifndef USECUDA
void ThermoMoist::getBuoyancySurf(Field3d *bfield)
{
  calcBuoyancyBot(bfield->data         , bfield->databot,
                  fields->sp["thl" ]->data, fields->sp["thl" ]->databot,
                  fields->sp["qt"]->data, fields->sp["qt"]->databot,
                  thvref, thvrefh);
  calcBuoyancyFluxBot(bfield->datafluxbot, fields->sp["thl"]->databot, fields->sp["thl"]->datafluxbot, fields->sp["qt"]->databot, fields->sp["qt"]->datafluxbot, thvrefh);
}
#endif

#ifndef USECUDA
void ThermoMoist::getBuoyancyFluxbot(Field3d *bfield)
{
  calcBuoyancyFluxBot(bfield->datafluxbot, fields->sp["thl"]->databot, fields->sp["thl"]->datafluxbot, fields->sp["qt"]->databot, fields->sp["qt"]->datafluxbot, thvrefh);
}
#endif

void ThermoMoist::getProgVars(std::vector<std::string> *list)
{
  list->push_back("thl");
  list->push_back("qt");
}

/**
 * This function calculates the hydrostatic pressure at full and half levels, 
 * with option to return base state profiles like reference density and temperature
 * Solves: dpi/dz=-g/thv with pi=cp*(p/p0)**(rd/cp)
 * @param pref Pointer to output hydrostatic pressure array (full level) 
 * @param prefh Pointer to output hydrostatic pressure array (half level) 
 * @param dn Pointer to output density array (full level) 
 * @param dnh Pointer to output density array (half level) 
 * @param thv Pointer to output virtual potential temperature array (full level) 
 * @param thvh Pointer to output virtual potential temperature array (half level) 
 * @param ex Pointer to output exner array (full level) 
 * @param exh Pointer to output exner array (half level) 
 * @param thlmean Pointer to input liq. water potential temperature array (horizontal mean, full level) 
 * @param qtmean Pointer to input tot. moisture mix. ratio  array (horizontal mean, full level) 
 */
void  ThermoMoist::calcBaseState(double * restrict pref,     double * restrict prefh,
                                 double * restrict rho,      double * restrict rhoh,
                                 double * restrict thv,      double * restrict thvh,
                                 double * restrict ex,       double * restrict exh,
                                 double * restrict thlmean,  double * restrict qtmean)
{
  const double rdcp = Rd/cp;

  const int kstart = grid->kstart;
  const int kend = grid->kend;

  double ssurf = constants::dhuge;
  double qtsurf = constants::dhuge;

  if(grid->swspatialorder == "2")
  {
    ssurf  = interp2(thlmean[kstart-1], thlmean[kstart]);
    qtsurf = interp2(qtmean[kstart-1],  qtmean[kstart]);
  }
  else if(grid->swspatialorder == "4")
  {
    ssurf  = interp4(thlmean[kstart-2], thlmean[kstart-1], thlmean[kstart], thlmean[kstart+1]);
    qtsurf = interp4(qtmean[kstart-2],  qtmean[kstart-1],  qtmean[kstart],  qtmean[kstart+1]);
  }

  double ql,si,qti,qli;

  // Calculate surface (half=kstart) values
  exh[kstart]   = exner(pbot);
  ql            = satAdjust(ssurf,qtsurf,pbot,exh[kstart]); 
  thvh[kstart]  = (ssurf + Lv*ql/(cp*exh[kstart])) * (1. - (1. - Rv/Rd)*qtsurf - Rv/Rd*ql);
  prefh[kstart] = pbot;
  rhoh[kstart]  = pbot / (Rd * exh[kstart] * thvh[kstart]);

  // First full grid level pressure
  pref[kstart] = pow((pow(pbot,rdcp) - grav * pow(p0,rdcp) * grid->z[kstart] / (cp * thvh[kstart])),(1./rdcp)); 

  for(int k=kstart+1; k<kend+1; k++)
  {
    // 1. Calculate values at full level below zh[k] 
    ex[k-1]  = exner(pref[k-1]);
    ql       = satAdjust(thlmean[k-1],qtmean[k-1],pref[k-1],ex[k-1]); 
    thv[k-1] = (thlmean[k-1] + Lv*ql/(cp*ex[k-1])) * (1. - (1. - Rv/Rd)*qtmean[k-1] - Rv/Rd*ql); 
    rho[k-1] = pref[k-1] / (Rd * ex[k-1] * thv[k-1]);
 
    // 2. Calculate half level pressure at zh[k] using values at z[k-1]
    prefh[k] = pow((pow(prefh[k-1],rdcp) - grav * pow(p0,rdcp) * grid->dz[k-1] / (cp * thv[k-1])),(1./rdcp));

    // 3. Interpolate conserved variables to zh[k] and calculate virtual temp and ql
    if(grid->swspatialorder == "2")
    {
      si     = interp2(thlmean[k-1],thlmean[k]);
      qti    = interp2(qtmean[k-1],qtmean[k]);
    }
    else if(grid->swspatialorder == "4")
    {
      si     = interp4(thlmean[k-2],thlmean[k-1],thlmean[k],thlmean[k+1]);
      qti    = interp4(qtmean[k-2],qtmean[k-1],qtmean[k],qtmean[k+1]);
    }

    exh[k]   = exner(prefh[k]);
    qli      = satAdjust(si,qti,prefh[k],exh[k]);
    thvh[k]  = (si + Lv*qli/(cp*exh[k])) * (1. - (1. - Rv/Rd)*qti - Rv/Rd*qli); 
    rhoh[k]  = prefh[k] / (Rd * exh[k] * thvh[k]); 

    // 4. Calculate full level pressure at z[k]
    pref[k]  = pow((pow(pref[k-1],rdcp) - grav * pow(p0,rdcp) * grid->dzh[k] / (cp * thvh[k])),(1./rdcp)); 
  }

  // Fill bottom and top full level ghost cells 
  if(grid->swspatialorder == "2")
  {
    pref[kstart-1] = 2.*prefh[kstart] - pref[kstart];
    pref[kend]     = 2.*prefh[kend]   - pref[kend-1];
  }
  else if(grid->swspatialorder == "4")
  {
    pref[kstart-1] = (8./3.)*prefh[kstart] - 2.*pref[kstart] + (1./3.)*pref[kstart+1];
    pref[kstart-2] = 8.*prefh[kstart]      - 9.*pref[kstart] + 2.*pref[kstart+1];
    pref[kend]     = (8./3.)*prefh[kend]   - 2.*pref[kend-1] + (1./3.)*pref[kend-2];
    pref[kend+1]   = 8.*prefh[kend]        - 9.*pref[kend-1] + 2.*pref[kend-2];
  }
}

void ThermoMoist::calcBuoyancyTend_2nd(double * restrict wt, double * restrict s, double * restrict qt, 
                                       double * restrict ph, double * restrict sh, double * restrict qth, double * restrict ql,
                                       double * restrict thvrefh)
{
  int ijk,jj,kk,ij;
  double tl, exnh;
  jj = grid->icells;
  kk = grid->ijcells;

  // CvH check the usage of the gravity term here, in case of scaled DNS we use one. But thermal expansion coeff??
  for(int k=grid->kstart+1; k<grid->kend; k++)
  {
    exnh = exner(ph[k]);
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
        ql[ij]  = qth[ij]-qsat(ph[k],tl);
      }
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ij  = i + j*jj;
        if(ql[ij]>0)   // already doesn't vectorize because of iteration in satAdjust()
        {
          ql[ij] = satAdjust(sh[ij], qth[ij], ph[k], exnh);
        }
        else
          ql[ij] = 0.;
      }
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ij  = i + j*jj;
        wt[ijk] += buoyancy(ph[k], exnh, sh[ij], qth[ij], ql[ij], thvrefh[k]);
      }
  }
}

void ThermoMoist::calcBuoyancyTend_4th(double * restrict wt, double * restrict s, double * restrict qt, 
                                       double * restrict ph, double * restrict sh, double * restrict qth, double * restrict ql,
                                       double * restrict thvrefh)
{
  int ijk,jj,ij;
  int kk1,kk2;
  double tl, exnh;

  jj  = grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  for(int k=grid->kstart+1; k<grid->kend; k++)
  {
    exnh = exner(ph[k]);
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
        ql[ij]  = qth[ij]-qsat(ph[k],tl);   
      }
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ij  = i + j*jj;
        if(ql[ij]>0)   // already doesn't vectorize because of iteration in satAdjust()
          ql[ij] = satAdjust(sh[ij], qth[ij], ph[k], exnh);
        else
          ql[ij] = 0.;
      }
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk1;
        ij  = i + j*jj;
        wt[ijk] += buoyancy(ph[k], exnh, sh[ij], qth[ij], ql[ij], thvrefh[k]);
      }
  }
}

void ThermoMoist::calcBuoyancy(double * restrict b, double * restrict s, double * restrict qt, double * restrict p, double * restrict ql,
                               double * restrict thvref)
{
  int ijk,jj,kk,ij;
  double tl, ex;
  jj = grid->icells;
  kk = grid->ijcells;

  for(int k=0; k<grid->kcells; k++)
  {
    ex = exner(p[k]);
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ij  = i + j*jj;
        tl  = s[ijk] * ex;
        ql[ij]  = qt[ijk]-qsat(p[k],tl);   // not real ql, just estimate
      }
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ij  = i + j*jj;
        if(ql[ij] > 0)
          ql[ij] = satAdjust(s[ijk], qt[ijk], p[k], ex);
        else
          ql[ij] = 0.;
      }
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ij  = i + j*jj;
        b[ijk] = buoyancy(p[k], ex, s[ijk], qt[ijk], ql[ij], thvref[k]);
      }
  }
}

void ThermoMoist::calcLiquidWater(double * restrict ql, double * restrict s, double * restrict qt, double * restrict p)
{
  int ijk,jj,kk;
  double ex;

  jj = grid->icells;
  kk = grid->ijcells;

  // Fill ghost cells with zeros to prevent problems in calculating ql or qlcore masks 
  for(int k=0; k<grid->kstart; k++)
  {
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ql[ijk] = 0.;
      }
  }

  for(int k=grid->kend; k<grid->kcells; k++)
  {
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ql[ijk] = 0.;
      }
  }

  // Calculate the ql field
  for(int k=grid->kstart; k<grid->kend; k++)
  {
    ex = exner(p[k]);
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ql[ijk] = satAdjust(s[ijk], qt[ijk], p[k], ex);
      }
  }
}

void ThermoMoist::calcN2(double * restrict N2, double * restrict s, double * restrict dzi, double * restrict thvref)
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
        N2[ijk] = grav/thvref[k]*0.5*(s[ijk+kk] - s[ijk-kk])*dzi[k];
      }
}

void ThermoMoist::calcBuoyancyBot(double * restrict b , double * restrict bbot,
                                  double * restrict s , double * restrict sbot,
                                  double * restrict qt, double * restrict qtbot,
                                  double * restrict thvref, double * restrict thvrefh)
{
  int ij,ijk,jj,kk,kstart;
  jj = grid->icells;
  kk = grid->ijcells;
  kstart = grid->kstart;

  // double thvref = thvs;

  // assume no liquid water at the lowest model level
  for(int j=0; j<grid->jcells; j++)
#pragma ivdep
    for(int i=0; i<grid->icells; i++)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      bbot[ij ] = buoyancyNoql(sbot[ij], qtbot[ij], thvrefh[kstart]);
      b   [ijk] = buoyancyNoql(s[ijk], qt[ijk], thvref[kstart]);
    }
}

void ThermoMoist::calcBuoyancyFluxBot(double * restrict bfluxbot, double * restrict sbot, double * restrict sfluxbot, double * restrict qtbot, double * restrict qtfluxbot,
                                      double * restrict thvrefh)
{
  int ij,jj,kstart;
  jj = grid->icells;
  kstart = grid->kstart;

  // double thvref = thvs;

  // assume no liquid water at the lowest model level
  for(int j=0; j<grid->jcells; j++)
#pragma ivdep
    for(int i=0; i<grid->icells; i++)
    {
      ij  = i + j*jj;
      bfluxbot[ij] = buoyancyFluxNoql(sbot[ij], sfluxbot[ij], qtbot[ij], qtfluxbot[ij], thvrefh[kstart]);
    }
}

void ThermoMoist::initStat()
{
  // Add variables to the statistics
  if(stats->getSwitch() == "1")
  {
    /* Add fixed base-state density and temperature profiles. Density should probably be in fields (?), but
       there the statistics are initialized before thermo->create() is called */
    stats->addFixedProf("rhoref",  "Full level basic state density", "kg m-3", "z",  fields->rhoref );
    stats->addFixedProf("rhorefh", "Half level basic state density", "kg m-3", "zh", fields->rhorefh);
    stats->addFixedProf("thvref",  "Full level basic state virtual potential temperature", "K", "z", thvref );
    stats->addFixedProf("thvrefh", "Half level basic state virtual potential temperature", "K", "zh",thvrefh);

    if(swupdatebasestate == 1)
    {
      stats->addProf("ph",   "Full level hydrostatic pressure", "Pa",     "z" );
      stats->addProf("phh",  "Half level hydrostatic pressure", "Pa",     "zh");
      stats->addProf("rho",  "Full level density",  "kg m-3", "z" );
      stats->addProf("rhoh", "Half level density",  "kg m-3", "zh");
    }
    else
    {
      stats->addFixedProf("ph",  "Full level hydrostatic pressure", "Pa", "z",  pref );
      stats->addFixedProf("phh", "Half level hydrostatic pressure", "Pa", "zh", prefh);
    }

    stats->addProf("b", "Buoyancy", "m s-2", "z");
    for(int n=2; n<5; ++n)
    {
      std::stringstream ss;
      ss << n;
      std::string sn = ss.str();
      stats->addProf("b"+sn, "Moment " +sn+" of the buoyancy", "(m s-2)"+sn,"z");
    }

    stats->addProf("bgrad", "Gradient of the buoyancy", "m s-3", "zh");
    stats->addProf("bw"   , "Turbulent flux of the buoyancy", "m2 s-3", "zh");
    stats->addProf("bdiff", "Diffusive flux of the buoyancy", "m2 s-3", "zh");
    stats->addProf("bflux", "Total flux of the buoyancy", "m2 s-3", "zh");

    stats->addProf("ql", "Liquid water mixing ratio", "kg kg-1", "z");
    stats->addProf("cfrac", "Cloud fraction", "-","z");

    stats->addTimeSeries("lwp", "Liquid water path", "kg m-2");
    stats->addTimeSeries("ccover", "Projected cloud cover", "-");
  }
}

void ThermoMoist::initCross()
{
  if(model->cross->getSwitch() == "1")
  {
    allowedcrossvars.push_back("b");
    allowedcrossvars.push_back("bbot");
    allowedcrossvars.push_back("bfluxbot");
    if(grid->swspatialorder == "4")
      allowedcrossvars.push_back("blngrad");
    allowedcrossvars.push_back("ql");
    allowedcrossvars.push_back("qlpath");

    // Check input list of cross variables (crosslist)
    std::vector<std::string>::iterator it=crosslist.begin();
    while(it != crosslist.end())
    {
      if(!std::count(allowedcrossvars.begin(),allowedcrossvars.end(),*it))
      {
        master->printWarning("WARNING field %s in [thermo][crosslist] is illegal\n", it->c_str());
        it = crosslist.erase(it);  // erase() returns iterator of next element..
      }
      else
        ++it;
    }

    // Sort crosslist to group ql and b variables
    std::sort(crosslist.begin(),crosslist.end());
  }
}

void ThermoMoist::initDump()
{
  if(model->dump->getSwitch() == "1")
  {
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
}

// INLINE FUNCTIONS
inline double ThermoMoist::buoyancy(const double p, const double exn, const double s, const double qt, const double ql, const double thvref)
{
  return grav * ((s + Lv*ql/(cp*exn)) * (1. - (1. - Rv/Rd)*qt - Rv/Rd*ql) - thvref) / thvref;
}

inline double ThermoMoist::buoyancyNoql(const double s, const double qt, const double thvref)
{
  return grav * (s * (1. - (1. - Rv/Rd)*qt) - thvref) / thvref;
}

inline double ThermoMoist::buoyancyFluxNoql(const double s, const double sflux, const double qt, const double qtflux, const double thvref)
{
  return grav/thvref * (sflux * (1. - (1.-Rv/Rd)*qt) - (1.-Rv/Rd)*s*qtflux);
}

inline double ThermoMoist::satAdjust(const double s, const double qt, const double p, const double exn)
{
  int niter = 0; //, nitermax = 5;
  double ql, tl, tnr_old = 1.e9, tnr, qs;
  tl = s * exn;
  tnr = tl;
  while (std::fabs(tnr-tnr_old)/tnr_old> 1e-5)// && niter < nitermax)
  {
    ++niter;
    tnr_old = tnr;
    qs = qsat(p,tnr);
    tnr = tnr - (tnr+(Lv/cp)*qs-tl-(Lv/cp)*qt)/(1+(std::pow(Lv,2)*qs)/ (Rv*cp*std::pow(tnr,2)));
  }
  ql = std::max(0.,qt - qs);
  return ql;
}

inline double ThermoMoist::exner(const double p)
{
  return pow((p/p0),(Rd/cp));
}

inline double ThermoMoist::exn2(const double p)
{
  double dp=p-p0;
  return (1+(dp*(ex1+dp*(ex2+dp*(ex3+dp*(ex4+dp*(ex5+dp*(ex6+ex7*dp)))))))); 
}

inline double ThermoMoist::qsat(const double p, const double T)
{
  return ep*esat(T)/(p-(1-ep)*esat(T));
}

inline double ThermoMoist::esat(const double T)
{
  const double x=std::max(-80.,T-T0);
  return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))));
}
