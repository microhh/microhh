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

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "fd.h"
#include "model.h"
#include "stats.h"
#include "cross.h"
#include "diff_smag2.h"

Fields::Fields(Model *modelin, Input *inputin)
{
  model  = modelin;
  grid   = model->grid;
  master = model->master;

  calcMeanProfs = false;

  // Initialize the pointers.
  rhoref  = 0;
  rhorefh = 0;
  umodel  = 0;
  vmodel  = 0;

  // Initialize GPU pointers
  rhoref_g  = 0;
  rhorefh_g = 0;

  // input parameters
  int nerror = 0;

  // obligatory parameters
  nerror += inputin->getItem(&visc, "fields", "visc", "");

  // read the name of the passive scalars
  std::vector<std::string> slist;
  nerror += inputin->getList(&slist, "fields", "slist", "");

  // initialize the scalars
  for(std::vector<std::string>::const_iterator it=slist.begin(); it!=slist.end(); ++it)
  {
    initPrognosticField(*it, *it, "-");
    nerror += inputin->getItem(&sp[*it]->visc, "fields", "svisc", *it);
  }

  // Read list of cross sections
  nerror += inputin->getList(&crosslist , "fields", "crosslist" , "");

  // Read list of 3d fiels to save
  nerror += inputin->getList(&dumplist , "fields", "dumplist" , "");

  if(nerror)
    throw 1;

  // initialize the basic set of fields
  initMomentumField(u, ut, "u", "U velocity", "m s-1");
  initMomentumField(v, vt, "v", "V velocity", "m s-1");
  initMomentumField(w, wt, "w", "Vertical velocity", "m s-1");
  initDiagnosticField("p", "Pressure", "Pa");
  initTmpField("tmp1", "", "");
  initTmpField("tmp2", "", "");
  initTmpField("tmp3", "", "");
  initTmpField("tmp4", "", "");
}

Fields::~Fields()
{
  // DEALLOCATE ALL THE FIELDS
  // deallocate the prognostic velocity fields
  for(FieldMap::iterator it=mp.begin(); it!=mp.end(); ++it)
    delete it->second;

  // deallocate the velocity tendency fields
  for(FieldMap::iterator it=mt.begin(); it!=mt.end(); ++it)
    delete it->second;

  // deallocate the prognostic scalar fields
  for(FieldMap::iterator it=sp.begin(); it!=sp.end(); ++it)
    delete it->second;

  // deallocate the scalar tendency fields
  for(FieldMap::iterator it=st.begin(); it!=st.end(); ++it)
    delete it->second;

  // deallocate the diagnostic scalars
  for(FieldMap::iterator it=sd.begin(); it!=sd.end(); ++it)
    delete it->second;

  // deallocate the tmp fields
  for(FieldMap::iterator it=atmp.begin(); it!=atmp.end(); ++it)
    delete it->second;

  // delete the arrays
  delete[] rhoref;
  delete[] rhorefh;
  delete[] umodel;
  delete[] vmodel;

#ifdef USECUDA
  clearDevice();
#endif
}

void Fields::init()
{
  // set the convenience pointers
  stats = model->stats;

  master->printMessage("Initializing fields\n");

  int nerror = 0;

  // ALLOCATE ALL THE FIELDS
  // allocate the prognostic velocity fields
  for(FieldMap::iterator it=mp.begin(); it!=mp.end(); ++it)
    nerror += it->second->init();

  // allocate the velocity tendency fields
  for(FieldMap::iterator it=mt.begin(); it!=mt.end(); ++it)
    nerror += it->second->init();

  // allocate the prognostic scalar fields
  for(FieldMap::iterator it=sp.begin(); it!=sp.end(); ++it)
    nerror += it->second->init();

  // allocate the scalar tendency fields
  for(FieldMap::iterator it=st.begin(); it!=st.end(); ++it)
    nerror += it->second->init();

  // allocate the diagnostic scalars
  for(FieldMap::iterator it=sd.begin(); it!=sd.end(); ++it)
    nerror += it->second->init();

  // allocate the tmp fields
  for(FieldMap::iterator it=atmp.begin(); it!=atmp.end(); ++it)
    nerror += it->second->init();

  if(nerror > 0)
    throw 1;

  // allocate the base density profiles
  rhoref  = new double[grid->kcells];
  rhorefh = new double[grid->kcells];

  // \TODO Define a reference density. Needs to be replaced once anelastic is there
  // BvS: Always init rhoref at 1 for situation with e.g. thermo=0? For anelastic, overwrite it.
  for(int k=0; k<grid->kcells; ++k)
  {
    rhoref[k] = 1.;
    rhorefh[k] = 1.; 
  }

  // allocate help arrays for statistics;
  umodel = new double[grid->kcells];
  vmodel = new double[grid->kcells];

  // Initialize at zero
  for(int k=0; k<grid->kcells; ++k)
  {
    umodel[k] = 0.;
    vmodel[k] = 0.; 
  }

  // Check different type of crosses and put them in their respective lists 
  for(FieldMap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    checkAddedCross(it->first, "",        &crosslist, &crosssimple);
    checkAddedCross(it->first, "lngrad",  &crosslist, &crosslngrad);
    checkAddedCross(it->first, "bot",     &crosslist, &crossbot);
    checkAddedCross(it->first, "top",     &crosslist, &crosstop);
    checkAddedCross(it->first, "fluxbot", &crosslist, &crossfluxbot);
    checkAddedCross(it->first, "fluxtop", &crosslist, &crossfluxtop);
  }

  for(FieldMap::const_iterator it=sd.begin(); it!=sd.end(); ++it)
  {
    checkAddedCross(it->first, "",       &crosslist, &crosssimple);
    checkAddedCross(it->first, "lngrad", &crosslist, &crosslngrad);
  }

  // If crosslist not empty, illegal variables or cross types were selected
  if(crosslist.size() > 0)
  {
    for(std::vector<std::string>::const_iterator it=crosslist.begin(); it!=crosslist.end(); ++it)
      master->printWarning("field %s in [fields][crosslist] is illegal\n", it->c_str());
  } 

  // Check if fields in dumplist are diagnostic fields, if not delete them and print warning
  std::vector<std::string>::iterator it=dumplist.begin();
  while(it != dumplist.end())
  {
    if(sd.count(*it) == 0)
    {
      master->printWarning("field %s in [fields][dumplist] is illegal\n", it->c_str());
      it = dumplist.erase(it);  // erase() returns iterator of next element..
    }
    else
      ++it;
  }
}

void Fields::checkAddedCross(std::string var, std::string type, std::vector<std::string> *crosslist, std::vector<std::string> *typelist)
{
  std::vector<std::string>::iterator position;
  
  position = std::find(crosslist->begin(), crosslist->end(), var + type);
  if(position != crosslist->end()) 
  {
    // don't allow lngrad in 2nd order mode
    if(!(type == "lngrad" && grid->swspatialorder == "2"))
    {
      typelist->push_back(var);
      crosslist->erase(position);
    }
  }
}

#ifndef USECUDA
void Fields::exec()
{
  // calculate the means for the prognostic scalars
  if(calcMeanProfs)
  {
    for(FieldMap::iterator it=sp.begin(); it!=sp.end(); ++it)
      grid->calcMean(it->second->datamean, it->second->data, grid->kcells);
  }
}
#endif

void Fields::getMask(Field3d *mfield, Field3d *mfieldh, Mask *m)
{
  if(m->name == "wplus")
    calcMask_wplus(mfield->data, mfieldh->data, mfieldh->databot, 
                   stats->nmask, stats->nmaskh, &stats->nmaskbot, w->data);
  else if(m->name == "wmin")                                                  
    calcMask_wmin(mfield->data, mfieldh->data, mfieldh->databot,
                  stats->nmask, stats->nmaskh, &stats->nmaskbot, w->data);
}

void Fields::calcMask_wplus(double * restrict mask, double * restrict maskh, double * restrict maskbot,
                            int * restrict nmask, int * restrict nmaskh, int * restrict nmaskbot,
                            double * restrict w)
{
  int ijk,ij,jj,kk,kstart;

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
        ntmp = (w[ijk] + w[ijk+kk]) > 0.;
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
        ntmp = w[ijk] > 0.;
        nmaskh[k] += ntmp;
        maskh[ijk] = (double)ntmp;
      }
  }

  // Set the mask for surface projected quantities
  // In this case: velocity at surface, so zero
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
}

void Fields::calcMask_wmin(double * restrict mask, double * restrict maskh, double * restrict maskbot,
                           int * restrict nmask, int * restrict nmaskh, int * restrict nmaskbot,
                           double * restrict w)
{
  int ijk,ij,jj,kk,kstart;

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
        ntmp = (w[ijk] + w[ijk+kk]) <= 0.;
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
        ntmp = w[ijk] <= 0.;
        nmaskh[k] += ntmp;
        maskh[ijk] = (double)ntmp;
      }
  }

  // Set the mask for surface projected quantities
  // In this case: velocity at surface, so zero
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
}

void Fields::execStats(Mask *m)
{
  // define locations
  const int uloc[] = {1,0,0};
  const int vloc[] = {0,1,0};
  const int wloc[] = {0,0,1};
  const int sloc[] = {0,0,0};

  const int uwloc[] = {1,0,1};
  const int vwloc[] = {0,1,1};

  const double NoOffset = 0.;

  // save the area coverage of the mask
  stats->calcArea(m->profs["area" ].data, sloc, stats->nmask );
  stats->calcArea(m->profs["areah"].data, wloc, stats->nmaskh);

  // start with the stats on the w location, to make the wmean known for the flux calculations
  stats->calcMean(m->profs["w"].data, w->data, NoOffset, wloc, atmp["tmp4"]->data, stats->nmaskh);
  for(int n=2; n<5; ++n)
  {
    std::stringstream ss;
    ss << n;
    std::string sn = ss.str();
    stats->calcMoment(w->data, m->profs["w"].data, m->profs["w"+sn].data, n, wloc,
                      atmp["tmp4"]->data, stats->nmaskh);
  }

  // calculate the stats on the u location
  // interpolate the mask horizontally onto the u coordinate
  grid->interpolate_2nd(atmp["tmp1"]->data, atmp["tmp3"]->data, sloc, uloc);
  stats->calcMean(m->profs["u"].data, u->data, grid->utrans, uloc, atmp["tmp1"]->data, stats->nmask);
  stats->calcMean(umodel            , u->data, NoOffset   , uloc, atmp["tmp1"]->data, stats->nmask);
  for(int n=2; n<5; ++n)
  {
    std::stringstream ss;
    ss << n;
    std::string sn = ss.str();
    stats->calcMoment(u->data, umodel, m->profs["u"+sn].data, n, uloc,
                      atmp["tmp1"]->data, stats->nmask);
  }

  // interpolate the mask on half level horizontally onto the u coordinate
  grid->interpolate_2nd(atmp["tmp1"]->data, atmp["tmp4"]->data, wloc, uwloc);
  if(grid->swspatialorder == "2")
  {
    stats->calcGrad_2nd(u->data, m->profs["ugrad"].data, grid->dzhi, uloc,
                        atmp["tmp1"]->data, stats->nmaskh);
    stats->calcFlux_2nd(u->data, umodel, w->data, m->profs["w"].data,
                        m->profs["uw"].data, atmp["tmp2"]->data, uloc,
                        atmp["tmp1"]->data, stats->nmaskh);
    if(model->diff->getName() == "smag2")
      stats->calcDiff_2nd(u->data, w->data, sd["evisc"]->data,
                          m->profs["udiff"].data, grid->dzhi,
                          u->datafluxbot, u->datafluxtop, 1., uloc,
                          atmp["tmp1"]->data, stats->nmaskh);
    else
      stats->calcDiff_2nd(u->data, m->profs["udiff"].data, grid->dzhi, visc, uloc,
                          atmp["tmp1"]->data, stats->nmaskh);

  }
  else if(grid->swspatialorder == "4")
  {
    stats->calcGrad_4th(u->data, m->profs["ugrad"].data, grid->dzhi4, uloc,
                        atmp["tmp1"]->data, stats->nmaskh);
    stats->calcFlux_4th(u->data, w->data, m->profs["uw"].data, atmp["tmp2"]->data, uloc,
                        atmp["tmp1"]->data, stats->nmaskh);
    stats->calcDiff_4th(u->data, m->profs["udiff"].data, grid->dzhi4, visc, uloc,
                        atmp["tmp1"]->data, stats->nmaskh);
  }

  // calculate the stats on the v location
  grid->interpolate_2nd(atmp["tmp1"]->data, atmp["tmp3"]->data, sloc, vloc);
  stats->calcMean(m->profs["v"].data, v->data, grid->vtrans, vloc, atmp["tmp1"]->data, stats->nmask);
  stats->calcMean(vmodel            , v->data, NoOffset   , vloc, atmp["tmp1"]->data, stats->nmask);
  for(int n=2; n<5; ++n)
  {
    std::stringstream ss;
    ss << n;
    std::string sn = ss.str();
    stats->calcMoment(v->data, vmodel, m->profs["v"+sn].data, n, vloc,
                      atmp["tmp1"]->data, stats->nmask);
  }

  // interpolate the mask on half level horizontally onto the u coordinate
  grid->interpolate_2nd(atmp["tmp1"]->data, atmp["tmp4"]->data, wloc, vwloc);
  if(grid->swspatialorder == "2")
  {
    stats->calcGrad_2nd(v->data, m->profs["vgrad"].data, grid->dzhi, vloc,
                        atmp["tmp1"]->data, stats->nmaskh);
    stats->calcFlux_2nd(v->data, vmodel, w->data, m->profs["w"].data,
                        m->profs["vw"].data, atmp["tmp2"]->data, vloc,
                        atmp["tmp1"]->data, stats->nmaskh);
    if(model->diff->getName() == "smag2")
      stats->calcDiff_2nd(v->data, w->data, sd["evisc"]->data,
                          m->profs["vdiff"].data, grid->dzhi,
                          v->datafluxbot, v->datafluxtop, 1., vloc,
                          atmp["tmp1"]->data, stats->nmaskh);
    else
      stats->calcDiff_2nd(v->data, m->profs["vdiff"].data, grid->dzhi, visc, vloc,
                          atmp["tmp1"]->data, stats->nmaskh);

  }
  else if(grid->swspatialorder == "4")
  {
    stats->calcGrad_4th(v->data, m->profs["vgrad"].data, grid->dzhi4, vloc,
                        atmp["tmp1"]->data, stats->nmaskh);
    stats->calcFlux_4th(v->data, w->data, m->profs["vw"].data, atmp["tmp2"]->data, vloc,
                        atmp["tmp1"]->data, stats->nmaskh);
    stats->calcDiff_4th(v->data, m->profs["vdiff"].data, grid->dzhi4, visc, vloc,
                        atmp["tmp1"]->data, stats->nmaskh);
  }

  // calculate stats for the prognostic scalars
  DiffSmag2 *diffptr = static_cast<DiffSmag2 *>(model->diff);
  for(FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
  {
    stats->calcMean(m->profs[it->first].data, it->second->data, NoOffset, sloc, atmp["tmp3"]->data, stats->nmask);
    for(int n=2; n<5; ++n)
    {
      std::stringstream ss;
      ss << n;
      std::string sn = ss.str();
      stats->calcMoment(it->second->data, m->profs[it->first].data, m->profs[it->first+sn].data, n, sloc,
                        atmp["tmp3"]->data, stats->nmask);
    }
    if(grid->swspatialorder == "2")
    {
      stats->calcGrad_2nd(it->second->data, m->profs[it->first+"grad"].data, grid->dzhi, sloc,
                          atmp["tmp4"]->data, stats->nmaskh);
      stats->calcFlux_2nd(it->second->data, m->profs[it->first].data, w->data, m->profs["w"].data,
                          m->profs[it->first+"w"].data, atmp["tmp1"]->data, sloc,
                          atmp["tmp4"]->data, stats->nmaskh);
      if(model->diff->getName() == "smag2")
        stats->calcDiff_2nd(it->second->data, w->data, sd["evisc"]->data,
                            m->profs[it->first+"diff"].data, grid->dzhi,
                            it->second->datafluxbot, it->second->datafluxtop, diffptr->tPr, sloc,
                            atmp["tmp4"]->data, stats->nmaskh);
      else
        stats->calcDiff_2nd(it->second->data, m->profs[it->first+"diff"].data, grid->dzhi, it->second->visc, sloc,
                            atmp["tmp4"]->data, stats->nmaskh);
    }
    else if(grid->swspatialorder == "4")
    {
      stats->calcGrad_4th(it->second->data, m->profs[it->first+"grad"].data, grid->dzhi4, sloc,
                          atmp["tmp4"]->data, stats->nmaskh);
      stats->calcFlux_4th(it->second->data, w->data, m->profs[it->first+"w"].data, atmp["tmp1"]->data, sloc,
                          atmp["tmp4"]->data, stats->nmaskh);
      stats->calcDiff_4th(it->second->data, m->profs[it->first+"diff"].data, grid->dzhi4, it->second->visc, sloc,
                          atmp["tmp4"]->data, stats->nmaskh);
    }
  }

  // calculate the total fluxes
  stats->addFluxes(m->profs["uflux"].data, m->profs["uw"].data, m->profs["udiff"].data);
  stats->addFluxes(m->profs["vflux"].data, m->profs["vw"].data, m->profs["vdiff"].data);
  for(FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
    stats->addFluxes(m->profs[it->first+"flux"].data, m->profs[it->first+"w"].data, m->profs[it->first+"diff"].data);

  // other statistics
  stats->calcMean(m->profs["p"].data, sd["p"]->data, NoOffset, sloc, atmp["tmp3"]->data, stats->nmask);

  if(model->diff->getName() == "smag2")
    stats->calcMean(m->profs["evisc"].data, sd["evisc"]->data, NoOffset, sloc, atmp["tmp3"]->data, stats->nmask);
}

void Fields::set_calcMeanProfs(bool sw)
{
  calcMeanProfs = sw;
}

void Fields::initMomentumField(Field3d *&fld, Field3d *&fldt, std::string fldname, std::string longname, std::string unit)
{
  if(mp.find(fldname)!=mp.end())
  {
    master->printError("\"%s\" already exists\n", fldname.c_str());
    throw 1;
  }

  // add a new prognostic momentum variable
  mp[fldname] = new Field3d(grid, master, fldname, longname, unit);

  // add a new tendency for momentum variable
  std::string fldtname  = fldname + "t";
  std::string tunit     = unit + "s-1";
  std::string tlongname = "Tendency of " + longname;
  mt[fldname] = new Field3d(grid, master, fldtname, tlongname, tunit);

  // TODO remove these from the model?
  fld  = mp[fldname];
  fldt = mt[fldname];

  // add the prognostic variable and its tendency to the collection
  // of all fields and tendencies
  a [fldname] = mp[fldname];
  ap[fldname] = mp[fldname];
  at[fldname] = mt[fldname];
}

void Fields::initPrognosticField(std::string fldname, std::string longname, std::string unit)
{
  if(sp.find(fldname)!=sp.end())
  {
    master->printError("\"%s\" already exists\n", fldname.c_str());
    throw 1;
  }
  
  // add a new scalar variable
  sp[fldname] = new Field3d(grid, master, fldname,longname, unit);

  // add a new tendency for scalar variable
  std::string fldtname  = fldname + "t";
  std::string tlongname = "Tendency of " + longname;
  std::string tunit     = unit + "s-1";
  st[fldname] = new Field3d(grid, master, fldtname,tlongname, tunit);

  // add the prognostic variable and its tendency to the collection
  // of all fields and tendencies
  a [fldname] = sp[fldname];
  ap[fldname] = sp[fldname];
  at[fldname] = st[fldname];
}

void Fields::initDiagnosticField(std::string fldname,std::string longname, std::string unit)
{
  if(sd.find(fldname)!=sd.end())
  {
    master->printError("\"%s\" already exists\n", fldname.c_str());
    throw 1;
  }

  sd[fldname] = new Field3d(grid, master, fldname, longname, unit);
  a [fldname] = sd[fldname];
}

void Fields::initTmpField(std::string fldname,std::string longname, std::string unit)
{
  if(atmp.find(fldname)!=atmp.end())
  {
    master->printError("\"%s\" already exists\n", fldname.c_str());
    throw 1;
  }

  atmp[fldname] = new Field3d(grid, master, fldname, longname, unit);
}

void Fields::create(Input *inputin)
{
  master->printMessage("Creating fields\n");
  
  int nerror = 0;
  
  // Randomnize the momentum
  for(FieldMap::iterator it=mp.begin(); it!=mp.end(); ++it)
    nerror += randomnize(inputin, it->first, it->second->data);
  
  // Randomnize the scalars
  for(FieldMap::iterator it=sp.begin(); it!=sp.end(); ++it)
    nerror += randomnize(inputin, it->first, it->second->data);
  
  // Add Vortices
  nerror += addVortexPair(inputin);
  
  // Add the mean profiles to the fields
  nerror += addMeanProf(inputin, "u", mp["u"]->data, grid->utrans);
  nerror += addMeanProf(inputin, "v", mp["v"]->data, grid->vtrans);
 
  for(FieldMap::iterator it=sp.begin(); it!=sp.end(); ++it)
    nerror += addMeanProf(inputin, it->first, it->second->data, 0.);
  
  // set w equal to zero at the boundaries, just to be sure
  int lbot = grid->kstart*grid->ijcells;
  int ltop = grid->kend  *grid->ijcells;
  for(int l=0; l<grid->ijcells; ++l)
  {
    w->data[lbot+l] = 0.;
    w->data[ltop+l] = 0.;
  }

  if(nerror)
    throw 1;
}

int Fields::randomnize(Input *inputin, std::string fld, double * restrict data)
{
  int nerror = 0;

  // set mpiid as random seed to avoid having the same field at all procs
  int static seed = 0;

  if(!seed)
  {
    nerror += inputin->getItem(&seed, "fields", "rndseed", "", 0);
    seed += master->mpiid + 2;
    std::srand(seed);
  }
  
  int ijk,jj,kk;
  int kendrnd;
  double rndfac;

  jj = grid->icells;
  kk = grid->ijcells;
  
  // look up the specific randomnizer variables
  nerror += inputin->getItem(&rndamp, "fields", "rndamp", fld, 0.);
  nerror += inputin->getItem(&rndz  , "fields", "rndz"  , fld, 0.);
  nerror += inputin->getItem(&rndexp, "fields", "rndexp", fld, 0.);

  // find the location of the randomizer height
  kendrnd = grid->kstart;
  while(grid->zh[kendrnd+1] < rndz)
    ++kendrnd;

  if(kendrnd > grid->kend)
  {
    master->printError("randomnizer height rndz (%f) higher than domain top (%f)\n", grid->z[kendrnd],grid->zsize);
    return 1;
  }
  
  if(kendrnd == grid->kstart)
    kendrnd = grid->kend;

  for(int k=grid->kstart; k<kendrnd; ++k)
  {
    rndfac = std::pow((rndz-grid->z [k])/rndz, rndexp);
    for(int j=grid->jstart; j<grid->jend; ++j)
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        data[ijk] = rndfac * rndamp * ((double) std::rand() / (double) RAND_MAX - 0.5);
      }
  }

  return nerror;
}

int Fields::addVortexPair(Input *inputin)
{
  int nerror = 0;

  // add a double vortex to the initial conditions
  const double pi = std::acos((double)-1.);
  int ijk, jj, kk;

  jj = grid->icells;
  kk = grid->ijcells;
    
  // optional parameters
  nerror += inputin->getItem(&vortexnpair, "fields", "vortexnpair", "", 0    );
  nerror += inputin->getItem(&vortexamp  , "fields", "vortexamp"  , "", 1.e-3);
  nerror += inputin->getItem(&vortexaxis , "fields", "vortexaxis" , "", "y"  );

  if(vortexnpair > 0)
  {
    if(vortexaxis == "y")
      for(int k=grid->kstart; k<grid->kend; ++k)
        for(int j=grid->jstart; j<grid->jend; ++j)
          for(int i=grid->istart; i<grid->iend; ++i)
          {
            ijk = i + j*jj + k*kk;
            u->data[ijk] +=  vortexamp*std::sin(vortexnpair*2.*pi*(grid->xh[i])/grid->xsize)*std::cos(pi*grid->z [k]/grid->zsize);
            w->data[ijk] += -vortexamp*std::cos(vortexnpair*2.*pi*(grid->x [i])/grid->xsize)*std::sin(pi*grid->zh[k]/grid->zsize);
          }
    else if(vortexaxis == "x")
      for(int k=grid->kstart; k<grid->kend; ++k)
        for(int j=grid->jstart; j<grid->jend; ++j)
          for(int i=grid->istart; i<grid->iend; ++i)
          {
            ijk = i + j*jj + k*kk;
            v->data[ijk] +=  vortexamp*std::sin(vortexnpair*2.*pi*(grid->yh[j])/grid->ysize)*std::cos(pi*grid->z [k]/grid->zsize);
            w->data[ijk] += -vortexamp*std::cos(vortexnpair*2.*pi*(grid->y [j])/grid->ysize)*std::sin(pi*grid->zh[k]/grid->zsize);
          }
  }
  
  return nerror;
}

int Fields::addMeanProf(Input *inputin, std::string fld, double * restrict data, double offset)
{
  int ijk, jj, kk;
  double proftemp[grid->kmax];

  jj = grid->icells;
  kk = grid->ijcells;
  
  if(inputin->getProf(proftemp, fld, grid->kmax))
    return 1;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        data[ijk] += proftemp[k-grid->kstart] - offset;
      }
      
  return 0;
}

void Fields::load(int n)
{
  const double NoOffset = 0.;

  int nerror = 0;

  for(FieldMap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    // the offset is kept at zero, otherwise bitwise identical restarts is not possible
    char filename[256];
    std::sprintf(filename, "%s.%07d", it->second->name.c_str(), n);
    master->printMessage("Loading \"%s\" ... ", filename);
    if(grid->loadField3d(it->second->data, atmp["tmp1"]->data, atmp["tmp2"]->data, filename, NoOffset))
    {
      master->printMessage("FAILED\n");
      ++nerror;
    }
    else
    {
      master->printMessage("OK\n");
    }  
  }

  // add the profiles to te statistics
  if(stats->getSwitch() == "1")
  {
    // add variables to the statistics
    stats->addProf(u->name, u->longname, u->unit, "z" );
    stats->addProf(v->name, v->longname, v->unit, "z" );
    stats->addProf(w->name, w->longname, w->unit, "zh" );
  
    for(FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->addProf(it->first,it->second->longname, it->second->unit, "z");
    stats->addProf(sd["p"]->name, sd["p"]->longname, sd["p"]->unit, "z");
 
    // CvH, shouldn't this call be in the diffusion class?
    if(model->diff->getName() == "smag2")
      stats->addProf(sd["evisc"]->name, sd["evisc"]->longname, sd["evisc"]->unit, "z");
  
    // moments
    for(int n=2; n<5; ++n)
    {
      std::stringstream ss;
      ss << n;
      std::string sn = ss.str();
      stats->addProf(u->name + sn,"Moment "+ sn + " of the " + u->longname,"(" + u->unit + ")"+sn, "z" );
      stats->addProf(v->name + sn,"Moment "+ sn + " of the " + v->longname,"(" + v->unit + ")"+sn, "z" );
      stats->addProf(w->name + sn,"Moment "+ sn + " of the " + w->longname,"(" + w->unit + ")"+sn, "zh" );
      for(FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
        stats->addProf(it->first + sn,"Moment "+ sn + " of the " + it->second->longname,"(" + it->second->unit + ")"+sn, "z" );
    }
  
    // gradients
    stats->addProf(u->name + "grad", "Gradient of the " + u->longname,"s-1","zh");
    stats->addProf(v->name + "grad", "Gradient of the " + v->longname,"s-1","zh");
    for(FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->addProf(it->first+"grad", "Gradient of the " + it->second->longname, it->second->unit + " m-1", "zh");
  
    // turbulent fluxes
    stats->addProf("uw", "Turbulent flux of the " + u->longname, "m2 s-2", "zh");
    stats->addProf("vw", "Turbulent flux of the " + v->longname, "m2 s-2", "zh");
    for(FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->addProf(it->first+"w", "Turbulent flux of the " + it->second->longname, it->second->unit + " m s-1", "zh");
  
    // Diffusive fluxes
    stats->addProf("udiff", "Diffusive flux of the " + u->longname, "m2 s-2", "zh");
    stats->addProf("vdiff", "Diffusive flux of the " + v->longname, "m2 s-2", "zh");
    for(FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->addProf(it->first+"diff", "Diffusive flux of the " + it->second->longname, it->second->unit + " m s-1", "zh");
  
    //Total fluxes
    stats->addProf("uflux", "Total flux of the " + u->longname, "m2 s-2", "zh");
    stats->addProf("vflux", "Total flux of the " + v->longname, "m2 s-2", "zh");
    for(FieldMap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->addProf(it->first+"flux", "Total flux of the " + it->second->longname, it->second->unit + " m s-1", "zh");
  }
 
  if(nerror)
    throw 1;
}

void Fields::save(int n)
{
  const double NoOffset = 0.;

  int nerror = 0;
  for(FieldMap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    char filename[256];
    std::sprintf(filename, "%s.%07d", it->second->name.c_str(), n);
    master->printMessage("Saving \"%s\" ... ", filename);

    // the offset is kept at zero, because otherwise bitwise identical restarts is not possible
    if(grid->saveField3d(it->second->data, atmp["tmp1"]->data, atmp["tmp2"]->data, filename, NoOffset))
    {
      master->printMessage("FAILED\n");
      ++nerror;
    }  
    else
    {
      master->printMessage("OK\n");
    }
  }

  if(nerror)
    throw 1;
}

#ifndef USECUDA
double Fields::checkMomentum()
{
  return calcMomentum_2nd(u->data, v->data, w->data, grid->dz);
}
#endif

#ifndef USECUDA
double Fields::checkTke()
{
  return calcTke_2nd(u->data, v->data, w->data, grid->dz);
}
#endif

#ifndef USECUDA
double Fields::checkMass()
{
  // CvH for now, do the mass check on the first scalar... Do we want to change this?
  FieldMap::const_iterator itProg=sp.begin();
  if(sp.begin() != sp.end())
    return calcMass(itProg->second->data, grid->dz);
  else
    return 0.;
}
#endif

double Fields::calcMass(double * restrict s, double * restrict dz)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->ijcells;

  double mass = 0;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        mass += s[ijk]*dz[k];
      }

  grid->getSum(&mass);

  mass /= (grid->itot*grid->jtot*grid->zsize);

  return mass;
}

double Fields::calcMomentum_2nd(double * restrict u, double * restrict v, double * restrict w, double * restrict dz)
{
  using fd::o2::interp2;

  int ijk,ii,jj,kk;

  ii = 1;
  jj = grid->icells;
  kk = grid->ijcells;

  double momentum;
  momentum = 0;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        momentum += (interp2(u[ijk], u[ijk+ii]) + interp2(v[ijk], v[ijk+jj]) + interp2(w[ijk], w[ijk+kk]))*dz[k];
      }

  grid->getSum(&momentum);

  momentum /= (grid->itot*grid->jtot*grid->zsize);

  return momentum;
}

double Fields::calcTke_2nd(double * restrict u, double * restrict v, double * restrict w, double * restrict dz)
{
  using fd::o2::interp2;

  int ijk,ii,jj,kk;

  ii = 1;
  jj = grid->icells;
  kk = grid->ijcells;

  double tke = 0;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        tke += ( interp2(u[ijk]*u[ijk], u[ijk+ii]*u[ijk+ii]) 
               + interp2(v[ijk]*v[ijk], v[ijk+jj]*v[ijk+jj]) 
               + interp2(w[ijk]*w[ijk], w[ijk+kk]*w[ijk+kk]))*dz[k];
      }

  grid->getSum(&tke);

  tke /= (grid->itot*grid->jtot*grid->zsize);
  tke *= 0.5;

  return tke;
}

void Fields::execCross()
{
  int nerror = 0;

  Cross *cross = model->cross;

  for(std::vector<std::string>::const_iterator it=crosssimple.begin(); it<crosssimple.end(); ++it)
    nerror += cross->crossSimple(a[*it]->data, atmp["tmp1"]->data, a[*it]->name);

  for(std::vector<std::string>::const_iterator it=crosslngrad.begin(); it<crosslngrad.end(); ++it)
    nerror += cross->crossLngrad(a[*it]->data, atmp["tmp1"]->data, atmp["tmp2"]->data, grid->dzi4, a[*it]->name + "lngrad");

  for(std::vector<std::string>::const_iterator it=crossfluxbot.begin(); it<crossfluxbot.end(); ++it)
    nerror += cross->crossPlane(a[*it]->datafluxbot, atmp["tmp1"]->data, a[*it]->name + "fluxbot");

  for(std::vector<std::string>::const_iterator it=crossfluxtop.begin(); it<crossfluxtop.end(); ++it)
    nerror += cross->crossPlane(a[*it]->datafluxtop, atmp["tmp1"]->data, a[*it]->name + "fluxtop");

  for(std::vector<std::string>::const_iterator it=crossbot.begin(); it<crossbot.end(); ++it)
    nerror += cross->crossPlane(a[*it]->databot, atmp["tmp1"]->data, a[*it]->name + "bot");

  for(std::vector<std::string>::const_iterator it=crosstop.begin(); it<crosstop.end(); ++it)
    nerror += cross->crossPlane(a[*it]->datatop, atmp["tmp1"]->data, a[*it]->name + "top");

  if(nerror)
    throw 1;
}
