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

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "model.h"
#include "stats.h"
#include "cross.h"
#include "diff_les2s.h"

#define NO_OFFSET 0.

cfields::cfields(cmodel *modelin)
{
  model  = modelin;
  grid   = model->grid;
  master = model->master;

  allocated = false;
  calcprofs = false;
}

cfields::~cfields()
{
  if(allocated)
  {
    // DEALLOCATE ALL THE FIELDS
    // deallocate the prognostic velocity fields
    for(fieldmap::iterator it=mp.begin(); it!=mp.end(); ++it)
      delete it->second;

    // deallocate the velocity tendency fields
    for(fieldmap::iterator it=mt.begin(); it!=mt.end(); ++it)
      delete it->second;

    // deallocate the prognostic scalar fields
    for(fieldmap::iterator it=sp.begin(); it!=sp.end(); ++it)
      delete it->second;

    // deallocate the scalar tendency fields
    for(fieldmap::iterator it=st.begin(); it!=st.end(); ++it)
      delete it->second;

    // deallocate the diagnostic scalars
    for(fieldmap::iterator it=sd.begin(); it!=sd.end(); ++it)
      delete it->second;

    // delete the arrays
    delete[] rhoref;
    delete[] umodel;
    delete[] vmodel;
  }
}

int cfields::readinifile(cinput *inputin)
{
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
    if(initpfld(*it,*it,"-"))
      return 1;
    nerror += inputin->getItem(&sp[*it]->visc, "fields", "svisc", *it);
  }

  // Read list of cross sections
  nerror += inputin->getList(&crosslist , "fields", "crosslist" , "");

  // initialize the basic set of fields
  nerror += initmomfld(u, ut, "u", "U velocity", "m s-1");
  nerror += initmomfld(v, vt, "v", "V velocity", "m s-1");
  nerror += initmomfld(w, wt, "w", "Vertical velocity", "m s-1");
  nerror += initdfld("p", "Pressure", "Pa");
  nerror += initdfld("tmp0", "", "");
  nerror += initdfld("tmp1", "", "");
  nerror += initdfld("tmp2", "", "");

  // \TODO check this later
  stats = model->stats;

  return nerror;
}

int cfields::init()
{
  if(master->mpiid == 0) std::printf("Initializing fields\n");

  int n = 0;

  // ALLOCATE ALL THE FIELDS
  // allocate the prognostic velocity fields
  for(fieldmap::iterator it=mp.begin(); it!=mp.end(); ++it)
    n += it->second->init();

  // allocate the velocity tendency fields
  for(fieldmap::iterator it=mt.begin(); it!=mt.end(); ++it)
    n += it->second->init();

  // allocate the prognostic scalar fields
  for(fieldmap::iterator it=sp.begin(); it!=sp.end(); ++it)
    n += it->second->init();

  // allocate the scalar tendency fields
  for(fieldmap::iterator it=st.begin(); it!=st.end(); ++it)
    n += it->second->init();

  // allocate the diagnostic scalars
  for(fieldmap::iterator it=sd.begin(); it!=sd.end(); ++it)
    n += it->second->init();

  if(n > 0)
    return 1;

  // \TODO Define a reference density. Needs to be replaced once anelastic is there
  rhoref = new double[grid->kcells];
  for (int k = grid->kstart; k<grid->kend; ++k)
    rhoref[k] = 1.;

  // allocate help arrays for statistics;
  umodel = new double[grid->kcells];
  vmodel = new double[grid->kcells];

  allocated = true;

  // Check different type of crosses and put them in their respective lists 
  for(fieldmap::iterator it=a.begin(); it!=a.end(); ++it)
  {
    checkaddcross(it->first, "",        &crosslist, &crosssimple);
    checkaddcross(it->first, "lngrad",  &crosslist, &crosslngrad);
    checkaddcross(it->first, "bot",     &crosslist, &crossbot);
    checkaddcross(it->first, "top",     &crosslist, &crosstop);
    checkaddcross(it->first, "fluxbot", &crosslist, &crossfluxbot);
    checkaddcross(it->first, "fluxtop", &crosslist, &crossfluxtop);
  }

  // If crosslist not empty, illegal variables or cross types were selected
  if(crosslist.size() > 0)
  {
    for(std::vector<std::string>::const_iterator it=crosslist.begin(); it!=crosslist.end(); ++it)
      if(master->mpiid == 0) std::printf("WARNING field %s in [fields][crosslist] is illegal\n", it->c_str());
  } 

  return 0;
}

int cfields::checkaddcross(std::string var, std::string type, std::vector<std::string> *crosslist, std::vector<std::string> *typelist)
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

  return 0;
}

int cfields::exec()
{
  // calculate the means for the prognostic scalars
  if(calcprofs)
  {
    for(fieldmap::iterator it=sp.begin(); it!=sp.end(); ++it)
      grid->calcmean(it->second->datamean, it->second->data, grid->kcells);
  }

  return 0;
}

int cfields::getfilter(cfield3d *ffield, filter *f)
{
  if(f->name == "wplus")
    calcfilterwplus(ffield->data, f->profs["area"].data, f->profs["areah"].data, stats->filtercount, w->data);
  else if(f->name == "wmin")                                                  
    calcfilterwmin (ffield->data, f->profs["area"].data, f->profs["areah"].data, stats->filtercount, w->data);
  return 0;
}

int cfields::calcfilterwplus(double * restrict fdata, double * restrict area, double * restrict areah,
                             int * restrict nfilter, double * restrict w)
{
  int ijk,jj,kk;
  int kstart,kend;

  jj = grid->icells;
  kk = grid->ijcells;
  kstart = grid->kstart;
  kend   = grid->kend;

  int ntmp;

  for(int k=grid->kstart; k<grid->kend; k++)
  {
    nfilter[k] = 0;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ntmp = (w[ijk] + w[ijk+kk]) > 0.;
        nfilter[k] += ntmp;
        fdata[ijk] = (double)ntmp;
      }
  }

  // set bc's for the filter (mirror)
  nfilter[kstart-1] = nfilter[kstart];
  nfilter[kend    ] = nfilter[kend-1];
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj + kstart*kk;
      fdata[ijk-kk] = fdata[ijk];
      ijk = i + j*jj + (kend-1)*kk;
      fdata[ijk+kk] = fdata[ijk];
    }

  grid->boundary_cyclic(fdata);

  int ijtot = grid->itot*grid->jtot;
  master->sum(nfilter, grid->kcells);

  for(int k=grid->kstart; k<grid->kend; k++)
    area[k] = (double)nfilter[k] / (double)ijtot;

  for(int k=grid->kstart; k<grid->kend+1; k++)
    areah[k] = 0.5*(area[k-1] + area[k]);

  return 0;
}

int cfields::calcfilterwmin(double * restrict fdata, double * restrict area, double * restrict areah,
                            int * restrict nfilter, double * restrict w)
{
  int ijk,jj,kk;
  int kstart,kend;

  jj = grid->icells;
  kk = grid->ijcells;
  kstart = grid->kstart;
  kend   = grid->kend;

  int ntmp;

  for(int k=grid->kstart; k<grid->kend; k++)
  {
    nfilter[k] = 0;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ntmp = (w[ijk] + w[ijk+kk]) <= 0.;
        nfilter[k] += ntmp;
        fdata[ijk] = (double)ntmp;
      }
  }

  // set bc's for the filter (mirror)
  nfilter[kstart-1] = nfilter[kstart];
  nfilter[kend    ] = nfilter[kend-1];
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj + kstart*kk;
      fdata[ijk-kk] = fdata[ijk];
      ijk = i + j*jj + (kend-1)*kk;
      fdata[ijk+kk] = fdata[ijk];
    }

  grid->boundary_cyclic(fdata);

  int ijtot = grid->itot*grid->jtot;
  master->sum(nfilter, grid->kcells);

  for(int k=grid->kstart; k<grid->kend; k++)
    area[k] = (double)nfilter[k] / (double)ijtot;

  for(int k=grid->kstart; k<grid->kend+1; k++)
    areah[k] = 0.5*(area[k-1] + area[k]);

  return 0;
}

int cfields::execstats(filter *f)
{
  // define locations
  const int uloc[] = {1,0,0};
  const int vloc[] = {0,1,0};
  const int wloc[] = {0,0,1};
  const int sloc[] = {0,0,0};

  // calculate the means
  stats->calcmean(u->data, f->profs["u"].data, grid->utrans, uloc, sd["tmp0"]->data, stats->filtercount);
  stats->calcmean(v->data, f->profs["v"].data, grid->vtrans, vloc, sd["tmp0"]->data, stats->filtercount);
  stats->calcmean(w->data, f->profs["w"].data, NO_OFFSET, wloc, sd["tmp0"]->data, stats->filtercount);
  for(fieldmap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
    stats->calcmean(it->second->data, f->profs[it->first].data, NO_OFFSET, sloc, sd["tmp0"]->data, stats->filtercount);

  stats->calcmean(s["p"]->data, f->profs["p"].data, NO_OFFSET, sloc, sd["tmp0"]->data, stats->filtercount);

  if(model->diff->getname() == "les2s")
    stats->calcmean(s["evisc"]->data, f->profs["evisc"].data, NO_OFFSET, sloc, sd["tmp0"]->data, stats->filtercount);

  // calculate model means without correction for transformation
  stats->calcmean(u->data, umodel, NO_OFFSET);
  stats->calcmean(v->data, vmodel, NO_OFFSET);

  // calculate the higher order moments
  for(int n=2; n<5; ++n)
  {
    std::stringstream ss;
    ss << n;
    std::string sn = ss.str();
    stats->calcmoment(u->data, umodel, f->profs["u"+sn].data, n, uloc, sd["tmp0"]->data, stats->filtercount);
    stats->calcmoment(v->data, vmodel, f->profs["v"+sn].data, n, vloc, sd["tmp0"]->data, stats->filtercount);
    stats->calcmoment(w->data, f->profs["w"].data, f->profs["w"+sn].data, n, wloc, sd["tmp0"]->data, stats->filtercount);
    for(fieldmap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->calcmoment(it->second->data, f->profs[it->first].data, f->profs[it->first+sn].data, n, sloc, sd["tmp0"]->data, stats->filtercount);
  }

  // calculate the gradients
  if(grid->swspatialorder == "2")
  {
    stats->calcgrad_2nd(u->data, f->profs["ugrad"].data, grid->dzhi, uloc,
                        sd["tmp0"]->data, stats->filtercount);
    stats->calcgrad_2nd(v->data, f->profs["vgrad"].data, grid->dzhi, vloc,
                        sd["tmp0"]->data, stats->filtercount);
    for(fieldmap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->calcgrad_2nd(it->second->data, f->profs[it->first+"grad"].data, grid->dzhi, sloc,
                          sd["tmp0"]->data, stats->filtercount);
  }
  else if(grid->swspatialorder == "4")
  {
    stats->calcgrad_4th(u->data, f->profs["ugrad"].data, grid->dzhi4, uloc,
                        sd["tmp0"]->data, stats->filtercount);
    stats->calcgrad_4th(v->data, f->profs["vgrad"].data, grid->dzhi4, vloc,
                        sd["tmp0"]->data, stats->filtercount);
    for(fieldmap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->calcgrad_4th(it->second->data, f->profs[it->first+"grad"].data, grid->dzhi4, sloc,
                          sd["tmp0"]->data, stats->filtercount);
  }

  // calculate the turbulent fluxes
  if(grid->swspatialorder == "2")
  {
    stats->calcflux_2nd(u->data, f->profs["u"].data, w->data, f->profs["w"].data,
                        f->profs["uw"].data, s["tmp1"]->data, uloc,
                        s["tmp0"]->data, stats->filtercount);
    stats->calcflux_2nd(v->data, f->profs["v"].data, w->data, f->profs["w"].data,
                        f->profs["vw"].data, s["tmp1"]->data, vloc,
                        s["tmp0"]->data, stats->filtercount);
    for(fieldmap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->calcflux_2nd(it->second->data, f->profs[it->first].data, w->data, f->profs["w"].data,
                          f->profs[it->first+"w"].data, s["tmp1"]->data, sloc,
                          s["tmp0"]->data, stats->filtercount);
  }
  else if(grid->swspatialorder == "4")
  {
    stats->calcflux_4th(u->data, w->data, f->profs["uw"].data, s["tmp1"]->data, 1, 0);
    stats->calcflux_4th(v->data, w->data, f->profs["vw"].data, s["tmp1"]->data, 0, 1);
    for(fieldmap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->calcflux_4th(it->second->data, w->data, f->profs[it->first+"w"].data, s["tmp1"]->data, 0, 0);
  }

  // calculate the diffusive fluxes
  if(grid->swspatialorder == "2")
  {
    if(model->diff->getname() == "les2s")
    {
      cdiff_les2s *diffptr = static_cast<cdiff_les2s *>(model->diff);
      stats->calcdiff_2nd(u->data, s["evisc"]->data, f->profs["udiff"].data, grid->dzhi, u->datafluxbot, u->datafluxtop, 1.);
      stats->calcdiff_2nd(v->data, s["evisc"]->data, f->profs["vdiff"].data, grid->dzhi, v->datafluxbot, v->datafluxtop, 1.);
      for(fieldmap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
        stats->calcdiff_2nd(it->second->data, s["evisc"]->data, f->profs[it->first+"diff"].data, grid->dzhi, it->second->datafluxbot, it->second->datafluxtop, diffptr->tPr);
    }
  }
  else if(grid->swspatialorder == "4")
  {
    stats->calcdiff_4th(u->data, f->profs["udiff"].data, grid->dzhi4, visc, uloc,
                        s["tmp0"]->data, stats->filtercount);
    stats->calcdiff_4th(v->data, f->profs["vdiff"].data, grid->dzhi4, visc, vloc,
                        s["tmp0"]->data, stats->filtercount);
    for(fieldmap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->calcdiff_4th(it->second->data, f->profs[it->first+"diff"].data, grid->dzhi4, it->second->visc, sloc,
                          s["tmp0"]->data, stats->filtercount);
  }

  // calculate the total fluxes
  stats->addfluxes(f->profs["uflux"].data, f->profs["uw"].data, f->profs["udiff"].data);
  stats->addfluxes(f->profs["vflux"].data, f->profs["vw"].data, f->profs["vdiff"].data);
  for(fieldmap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
    stats->addfluxes(f->profs[it->first+"flux"].data, f->profs[it->first+"w"].data, f->profs[it->first+"diff"].data);

  return 0;
}

int cfields::setcalcprofs(bool sw)
{
  calcprofs = sw;
  return 0;
}

int cfields::initmomfld(cfield3d *&fld, cfield3d *&fldt, std::string fldname, std::string longname, std::string unit)
{
  if (mp.find(fldname)!=mp.end())
  {
    std::printf("ERROR \"%s\" already exists\n", fldname.c_str());
    return 1;
  }

  // add a new prognostic momentum variable
  mp[fldname] = new cfield3d(grid, master, fldname, longname, unit);

  // add a new tendency for momentum variable
  std::string fldtname  = fldname + "t";
  std::string tunit     = unit + "s-1";
  std::string tlongname = "Tendency of " + longname;
  mt[fldname] = new cfield3d(grid, master, fldtname, tlongname, tunit);

  // TODO remove these from the model?
  fld  = mp[fldname];
  fldt = mt[fldname];

  // add the prognostic variable and its tendency to the collection
  // of all fields and tendencies
  a [fldname] = mp[fldname];
  ap[fldname] = mp[fldname];
  at[fldname] = mt[fldname];

  return 0;
}

int cfields::initpfld(std::string fldname, std::string longname, std::string unit)
{
  if (s.find(fldname)!=s.end())
  {
    std::printf("ERROR \"%s\" already exists\n", fldname.c_str());
    return 1;
  }
  
  // add a new scalar variable
  sp[fldname] = new cfield3d(grid, master, fldname,longname, unit);

  // add a new tendency for scalar variable
  std::string fldtname  = fldname + "t";
  std::string tlongname = "Tendency of " + longname;
  std::string tunit     = unit + "s-1";
  st[fldname] = new cfield3d(grid, master, fldtname,tlongname, tunit);

  // add the prognostic variable and its tendency to the collection
  // of all fields and tendencies
  a [fldname] = sp[fldname];
  s [fldname] = sp[fldname];
  ap[fldname] = sp[fldname];
  at[fldname] = st[fldname];

  return 0;
}

int cfields::initdfld(std::string fldname,std::string longname, std::string unit)
{
  if (s.find(fldname)!=s.end())
  {
    std::printf("ERROR \"%s\" already exists\n", fldname.c_str());
    return 1;
  }

  sd[fldname] = new cfield3d(grid, master, fldname, longname, unit);
  s [fldname] = sd[fldname];
  a [fldname] = sd[fldname];

  return 0;  
}

int cfields::create(cinput *inputin)
{
  if(master->mpiid == 0) std::printf("Creating fields\n");
  
  int n = 0;
  
  // Randomnize the momentum
  for(fieldmap::iterator it=mp.begin(); it!=mp.end(); ++it)
    n += randomnize(inputin, it->first, it->second->data);
  
  // Randomnize the scalars
  for(fieldmap::iterator it=sp.begin(); it!=sp.end(); ++it)
    n += randomnize(inputin, it->first, it->second->data);
  
  // Add Vortices
  n += addvortexpair(inputin);
  
  // Add the mean profiles to the fields
  n += addmeanprofile(inputin, "u", mp["u"]->data, grid->utrans);
  n += addmeanprofile(inputin, "v", mp["v"]->data, grid->vtrans);
 
  for(fieldmap::iterator it=sp.begin(); it!=sp.end(); ++it)
    n += addmeanprofile(inputin, it->first, it->second->data, 0.);
  
  // set w equal to zero at the boundaries, just to be sure
  int lbot = grid->kstart*grid->icells*grid->jcells;
  int ltop = grid->kend  *grid->icells*grid->jcells;
  for(int l=0; l<grid->icells*grid->jcells; ++l)
  {
    w->data[lbot+l] = 0.;
    w->data[ltop+l] = 0.;
  }
  
  return (n>0);
}

int cfields::randomnize(cinput *inputin, std::string fld, double * restrict data)
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
  kk = grid->icells*grid->jcells;
  
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
    printf("ERROR: randomnizer height rndz (%f) higher than domain top (%f)\n", grid->z[kendrnd],grid->zsize);
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

int cfields::addvortexpair(cinput *inputin)
{
  int nerror = 0;

  // add a double vortex to the initial conditions
  const double pi = std::acos((double)-1.);
  int ijk, jj, kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;
    
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

int cfields::addmeanprofile(cinput *inputin, std::string fld, double * restrict data, double offset)
{
  int ijk, jj, kk;
  double proftemp[grid->kmax];

  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  
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

int cfields::load(int n)
{
  int nerror = 0;

  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    // the offset is kept at zero, otherwise bitwise identical restarts is not possible
    char filename[256];
    std::sprintf(filename, "%s.%07d", it->second->name.c_str(), n);
    if(master->mpiid == 0) std::printf("Loading \"%s\" ... ", filename);
    if(grid->loadfield3d(it->second->data, sd["tmp1"]->data, sd["tmp2"]->data, filename, NO_OFFSET))
    {
      if(master->mpiid == 0) std::printf("FAILED\n");
      ++nerror;
    }
    else
    {
      if(master->mpiid == 0) std::printf("OK\n");
    }  
  }

  // add the profiles to te statistics
  if(stats->getsw() == "1")
  {
    // add variables to the statistics
    stats->addprof(u->name, u->longname, u->unit, "z" );
    stats->addprof(v->name, v->longname, v->unit, "z" );
    stats->addprof(w->name, w->longname, w->unit, "zh" );
  
    for(fieldmap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->addprof(it->first,it->second->longname, it->second->unit, "z");
    stats->addprof(sd["p"]->name, sd["p"]->longname, sd["p"]->unit, "z");
  
    if(model->diff->getname() == "les2s")
      stats->addprof(sd["evisc"]->name, sd["evisc"]->longname, sd["evisc"]->unit, "z");
  
    // moments
    for(int n=2; n<5; ++n)
    {
      std::stringstream ss;
      ss << n;
      std::string sn = ss.str();
      stats->addprof(u->name + sn,"Moment "+ sn + " of the " + u->longname,"(" + u->unit + ")"+sn, "z" );
      stats->addprof(v->name + sn,"Moment "+ sn + " of the " + v->longname,"(" + v->unit + ")"+sn, "z" );
      stats->addprof(w->name + sn,"Moment "+ sn + " of the " + w->longname,"(" + w->unit + ")"+sn, "zh" );
      for(fieldmap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
        stats->addprof(it->first + sn,"Moment "+ sn + " of the " + it->second->longname,"(" + it->second->unit + ")"+sn, "z" );
    }
  
    // gradients
    stats->addprof(u->name + "grad", "Gradient of the " + u->longname,"s-1","zh");
    stats->addprof(v->name + "grad", "Gradient of the " + v->longname,"s-1","zh");
    for(fieldmap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->addprof(it->first+"grad", "Gradient of the " + it->second->longname, it->second->unit + " m-1", "zh");
  
    // turbulent fluxes
    stats->addprof("uw", "Turbulent flux of the " + u->longname, "m2 s-2", "zh");
    stats->addprof("vw", "Turbulent flux of the " + v->longname, "m2 s-2", "zh");
    for(fieldmap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->addprof(it->first+"w", "Turbulent flux of the " + it->second->longname, it->second->unit + " m s-1", "zh");
  
    // Diffusive fluxes
    stats->addprof("udiff", "Diffusive flux of the " + u->longname, "m2 s-2", "zh");
    stats->addprof("vdiff", "Diffusive flux of the " + v->longname, "m2 s-2", "zh");
    for(fieldmap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->addprof(it->first+"diff", "Diffusive flux of the " + it->second->longname, it->second->unit + " m s-1", "zh");
  
    //Total fluxes
    stats->addprof("uflux", "Total flux of the " + u->longname, "m2 s-2", "zh");
    stats->addprof("vflux", "Total flux of the " + v->longname, "m2 s-2", "zh");
    for(fieldmap::const_iterator it=sp.begin(); it!=sp.end(); ++it)
      stats->addprof(it->first+"flux", "Total flux of the " + it->second->longname, it->second->unit + " m s-1", "zh");
  }
  
  return nerror;
}

int cfields::save(int n)
{
  int nerror = 0;
  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    char filename[256];
    std::sprintf(filename, "%s.%07d", it->second->name.c_str(), n);
    if(master->mpiid == 0) std::printf("Saving \"%s\" ... ", filename);

    // the offset is kept at zero, because otherwise bitwise identical restarts is not possible
    if(grid->savefield3d(it->second->data, sd["tmp1"]->data, sd["tmp2"]->data, filename, NO_OFFSET))
    {
      if(master->mpiid == 0) std::printf("FAILED\n");
      ++nerror;
    }  
    else
    {
      if(master->mpiid == 0) std::printf("OK\n");
    }
  }

  return nerror;
}

double cfields::checkmom()
{
  return calcmom_2nd(u->data, v->data, w->data, grid->dz);
}

double cfields::checktke()
{
  return calctke_2nd(u->data, v->data, w->data, grid->dz);
}

double cfields::checkmass()
{
  // CvH for now, do the mass check on the first scalar... Do we want to change this?
  fieldmap::iterator itProg=sp.begin();
  if(sp.begin() != sp.end())
    return calcmass(itProg->second->data, grid->dz);
  else
    return 0.;
}

double cfields::calcmass(double * restrict s, double * restrict dz)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double mass = 0;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        mass += s[ijk]*dz[k];
      }

  grid->getsum(&mass);

  mass /= (grid->itot*grid->jtot*grid->zsize);

  return mass;
}

double cfields::calcmom_2nd(double * restrict u, double * restrict v, double * restrict w, double * restrict dz)
{
  int ijk,ii,jj,kk;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

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

  grid->getsum(&momentum);

  momentum /= (grid->itot*grid->jtot*grid->zsize);

  return momentum;
}

double cfields::calctke_2nd(double * restrict u, double * restrict v, double * restrict w, double * restrict dz)
{
  int ijk,ii,jj,kk;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

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

  grid->getsum(&tke);

  tke /= (grid->itot*grid->jtot*grid->zsize);
  tke *= 0.5;

  return tke;
}

int cfields::execcross()
{
  int nerror = 0;

  for(std::vector<std::string>::iterator it=crosssimple.begin(); it<crosssimple.end(); ++it)
    nerror += model->cross->crosssimple(a[*it]->data, s["tmp1"]->data, a[*it]->name);

  for(std::vector<std::string>::iterator it=crosslngrad.begin(); it<crosslngrad.end(); ++it)
    nerror += model->cross->crosslngrad(a[*it]->data, s["tmp1"]->data, s["tmp2"]->data, grid->dzi4, a[*it]->name + "lngrad");

  for(std::vector<std::string>::iterator it=crossfluxbot.begin(); it<crossfluxbot.end(); ++it)
    nerror += model->cross->crossplane(a[*it]->datafluxbot, s["tmp1"]->data, a[*it]->name + "fluxbot");

  for(std::vector<std::string>::iterator it=crossfluxtop.begin(); it<crossfluxtop.end(); ++it)
    nerror += model->cross->crossplane(a[*it]->datafluxtop, s["tmp1"]->data, a[*it]->name + "fluxtop");

  for(std::vector<std::string>::iterator it=crossbot.begin(); it<crossbot.end(); ++it)
    nerror += model->cross->crossplane(a[*it]->databot, s["tmp1"]->data, a[*it]->name + "bot");

  for(std::vector<std::string>::iterator it=crosstop.begin(); it<crosstop.end(); ++it)
    nerror += model->cross->crossplane(a[*it]->datatop, s["tmp1"]->data, a[*it]->name + "top");

  return nerror; 
}

inline double cfields::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}

