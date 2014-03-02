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
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "thermo_moist.h"
#include "defines.h"
#include "model.h"
#include "diff_les2s.h"
#include "timeloop.h"
#include <netcdfcpp.h>

#define NO_OFFSET 0.

cstats::cstats(cmodel *modelin)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  master = model->master;

  allocated   = false;
  initialized = false;
}

cstats::~cstats()
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

int cstats::readinifile(cinput *inputin)
{
  int nerror = 0;
  nerror += inputin->getItem(&sampletime, "stats", "sampletime", "");
  return nerror;
}

int cstats::init(double ifactor)
{
  isampletime = (unsigned long)(ifactor * sampletime);

  umodel = new double[grid->kcells];
  vmodel = new double[grid->kcells];

  allocated = true;

  // set the number of stats to zero
  nstats = 0;

  return 0;
}

int cstats::create(int n)
{
  int nerror = 0;

  // create a NetCDF file for the statistics
  if(master->mpiid == 0)
  {
    char filename[256];
    std::sprintf(filename, "%s.%07d.nc", master->simname.c_str(), n);
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
  master->broadcast(&nerror, 1);
  if(nerror)
    return 1;

  // create dimensions
  if(master->mpiid == 0)
  {
    z_dim  = dataFile->add_dim("z" , grid->kmax);
    zh_dim = dataFile->add_dim("zh", grid->kmax+1);
    t_dim  = dataFile->add_dim("t");

    // create variables belonging to dimensions
    iter_var = dataFile->add_var("iter", ncInt   , t_dim );
    nerror+= iter_var->add_att("units", "-");
    nerror+= iter_var->add_att("longname", "Iteration oumber");
    t_var    = dataFile->add_var("t"   , ncDouble, t_dim );
    nerror+= t_var->add_att("units", "s");
    nerror+= t_var->add_att("longname", "Time");
    z_var    = dataFile->add_var("z"   , ncDouble, z_dim );
    nerror+= z_var->add_att("units", "m");
    nerror+= z_var->add_att("longname", "Full level height");
    zh_var   = dataFile->add_var("zh"  , ncDouble, zh_dim);
    nerror+= z_var->add_att("units", "m");
    nerror+= z_var->add_att("longname", "Half level height");
  }

  // means
  // addprof("u", "z" );
  // addprof("v", "z" );
  // addprof("w", "zh");
  // for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  //   addprof(it->first, "z");
  //   addprof("evisc", "z");
  // addprof("p", "z");

  // in case of moisture, add ql prof
  if(model->thermo->getname() == "moist")
  {
  //     addprof("ql", "z");
  //     addprof("cfrac", "z");
  }

  // // 2nd order
  // addprof("u2", "z" );
  // addprof("v2", "z" );
  // addprof("w2", "zh");
  // for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  //   addprof(it->first+"2", "z");

  // 3rd order
  // addprof("u3", "z" );
  // addprof("v3", "z" );
  // addprof("w3", "zh");
  // for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  //   addprof(it->first+"3", "z");

  // // 4th order
  // addprof("u4", "z" );
  // addprof("v4", "z" );
  // addprof("w4", "zh");
  // for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  //   addprof(it->first+"4", "z");

  // // gradients
  // addprof("ugrad", "zh");
  // addprof("vgrad", "zh");
  // for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  //   addprof(it->first+"grad", "zh");

  // // turbulent fluxes
  // addprof("uw", "zh");
  // addprof("vw", "zh");
  // for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  //   addprof(it->first+"w", "zh");

  // // diffusive fluxes
  // addprof("udiff", "zh");
  // addprof("vdiff", "zh");
  // for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  //   addprof(it->first+"diff", "zh");

  // // total fluxes
  // addprof("uflux", "zh");
  // addprof("vflux", "zh");
  // for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  //   addprof(it->first+"flux", "zh");

  if(master->mpiid == 0)
  {
    // save the grid variables
    z_var ->put(&grid->z [grid->kstart], grid->kmax  );
    zh_var->put(&grid->zh[grid->kstart], grid->kmax+1);

    dataFile->sync();
  }

  return 0;
}

unsigned long cstats::gettimelim(unsigned long itime)
{
  unsigned long idtlim = isampletime -  itime % isampletime;
  return idtlim;
}

int cstats::dostats()
{
  // check if time for execution
  if(model->timeloop->itime % isampletime != 0)
    return 0;

  // write message in case stats is triggered
  if(master->mpiid == 0) std::printf("Saving stats for time %f\n", model->timeloop->time);

  // return true such that stats are computed
  return 1;
}

int cstats::exec(int iteration, double time, unsigned long itime)
{
  // check if time for execution
  if(itime % isampletime != 0)
    return 0;

  // if(master->mpiid == 0) std::printf("Saving stats for time %f\n", time);

  if(master->mpiid == 0)
  {
    t_var   ->put_rec(&time     , nstats);
    iter_var->put_rec(&iteration, nstats);
  }

  // PROFILES
  // calculate means
  // calcmean(fields->u->data, profs["u"].data, grid->utrans);
  // calcmean(fields->v->data, profs["v"].data, grid->vtrans);
  // calcmean(fields->w->data, profs["w"].data, NO_OFFSET);
  // for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  //   calcmean(it->second->data, profs[it->first].data, NO_OFFSET);

  // calcmean(fields->s["p"]->data, profs["p"].data, NO_OFFSET);
  //   calcmean(fields->s["evisc"]->data, profs["evisc"].data, NO_OFFSET);

  // in case of moisture, calc ql mean
  if(model->thermo->getname() == "moist")
  {
    // use a static cast to get access to the thermo moist functions
    //static_cast<cthermo_moist *>(model->thermo)->getql(fields->s["tmp1"], fields->s["tmp2"]);
    static_cast<cthermo_moist *>(model->thermo)->getthermofield(fields->s["tmp1"], fields->s["tmp2"],"ql");
//     calcmean (fields->s["tmp1"]->data, profs["ql"].data, NO_OFFSET);
//     calccount(fields->s["tmp1"]->data, profs["cfrac"].data, 0.);
  }

  /*
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

  calcgrad(fields->u->data, profs["ugrad"].data, grid->dzhi);
  calcgrad(fields->v->data, profs["vgrad"].data, grid->dzhi);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    calcgrad(it->second->data, profs[it->first+"grad"].data, grid->dzhi);

  // calculate turbulent fluxes
  calcflux(fields->u->data, fields->w->data, profs["uw"].data, fields->s["tmp1"]->data, 1, 0);
  calcflux(fields->v->data, fields->w->data, profs["vw"].data, fields->s["tmp1"]->data, 0, 1);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    calcflux(it->second->data, fields->w->data, profs[it->first+"w"].data, fields->s["tmp1"]->data, 0, 0);

  // calculate diffusive fluxes
  // TODO find a prettier solution for this cast later
  cdiff_les2s *diffptr = static_cast<cdiff_les2s *>(model->diff);
  calcdiff(fields->u->data, fields->s["evisc"]->data, profs["udiff"].data, grid->dzhi, fields->u->datafluxbot, fields->u->datafluxtop, 1.);
  calcdiff(fields->v->data, fields->s["evisc"]->data, profs["vdiff"].data, grid->dzhi, fields->v->datafluxbot, fields->v->datafluxtop, 1.);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    calcdiff(it->second->data, fields->s["evisc"]->data, profs[it->first+"diff"].data, grid->dzhi, it->second->datafluxbot, it->second->datafluxtop, diffptr->tPr);

  addfluxes(profs["uflux"].data, profs["uw"].data, profs["udiff"].data);
  addfluxes(profs["vflux"].data, profs["vw"].data, profs["vdiff"].data);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    addfluxes(profs[it->first+"flux"].data, profs[it->first+"w"].data, profs[it->first+"diff"].data);
  */

  // put the data into the NetCDF file
  if(master->mpiid == 0)
  {
    for(profmap::const_iterator it=profs.begin(); it!=profs.end(); ++it)
      profs[it->first].ncvar->put_rec(&profs[it->first].data[grid->kstart], nstats);

    // sync the data
    dataFile->sync();
  }

  ++nstats;

  return 0;
}

int cstats::addprof(std::string name, std::string longname, std::string unit, std::string zloc)
{
  int nerror = 0;
  // create the NetCDF variable
  if(master->mpiid == 0)
  {
    if(zloc == "z")
      profs[name].ncvar = dataFile->add_var(name.c_str(), ncDouble, t_dim, z_dim );
    else if(zloc == "zh")
      profs[name].ncvar = dataFile->add_var(name.c_str(), ncDouble, t_dim, zh_dim);
    nerror+=profs[name].ncvar->add_att("units", unit.c_str());
    nerror+=profs[name].ncvar->add_att("longname", longname.c_str());
    nerror+=profs[name].ncvar->add_att("_FillValue", NC_FILL_DOUBLE);
  }

  // and allocate the memory and initialize at zero
  profs[name].data = new double[grid->kcells];
  for(int k=0; k<grid->kcells; ++k)
    profs[name].data[k] = 0.;
  std::printf("Prof %s\n", name.c_str());
  return (nerror>0);
}

// COMPUTATIONAL KERNELS BELOW
int cstats::calcmean(double * restrict data, double * restrict prof, double offset)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int k=0; k<grid->kcells; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += data[ijk] + offset;
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=0; k<grid->kcells; ++k)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

// COMPUTATIONAL KERNELS BELOW
int cstats::calccount(double * restrict data, double * restrict prof, double threshold)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int k=0; k<grid->kcells; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj + k*kk;
        if(data[ijk]>threshold)
        {
          prof[k] += 1.;
        }
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=0; k<grid->kcells; ++k)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats::calcmoment(double * restrict data, double * restrict datamean, double * restrict prof, double power, int a)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  
  for(int k=grid->kstart; k<grid->kend+a; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += std::pow(data[ijk]-datamean[k], power);
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+a; ++k)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats::calcflux(double * restrict data, double * restrict w, double * restrict prof, double * restrict tmp1, int locx, int locy)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  // set a pointer to the field that contains w, either interpolated or the original
  double * restrict calcw = w;
  if(locx == 1)
  {
    grid->interpolatex_2nd(tmp1, w, 0);
    calcw = tmp1;
  }
  else if(locy == 1)
  {
    grid->interpolatey_2nd(tmp1, w, 0);
    calcw = tmp1;
  }
  
  for(int k=grid->kstart; k<grid->kend+1; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += 0.5*(data[ijk-kk]+data[ijk])*calcw[ijk];
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+1; ++k)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats::calcgrad(double * restrict data, double * restrict prof, double * restrict dzhi)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  
  for(int k=grid->kstart; k<grid->kend+1; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += (data[ijk]-data[ijk-kk])*dzhi[k];
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+1; ++k)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats::calcdiff(double * restrict data, double * restrict evisc, double * restrict prof, double * restrict dzhi, double * restrict fluxbot, double * restrict fluxtop, double tPr)
{
  int ijk,ij,jj,kk,kstart,kend;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  kstart = grid->kstart;
  kend   = grid->kend;

  // CvH add horizontal interpolation for u and v and interpolate the eddy viscosity properly
  // bottom boundary
  prof[kstart] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij = i + j*jj;
      prof[kstart] += fluxbot[ij];
    }

  for(int k=grid->kstart+1; k<grid->kend; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += -0.5*(evisc[ijk-kk]+evisc[ijk])/tPr*(data[ijk]-data[ijk-kk])*dzhi[k];
      }
  }

  // top boundary
  prof[kend] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij = i + j*jj;
      prof[kend] += fluxtop[ij];
    }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+1; ++k)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats::addfluxes(double * restrict flux, double * restrict turb, double * restrict diff)
{
  for(int k=grid->kstart; k<grid->kend+1; ++k)
    flux[k] = turb[k] + diff[k];

  return 0;
}

