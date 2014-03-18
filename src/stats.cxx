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

  // set the pointers to NULL
  umodel = NULL;
  vmodel = NULL;

  filtercount = NULL;
}

cstats::~cstats()
{
  delete[] umodel;
  delete[] vmodel;
  delete[] filtercount;

  // delete the profiles
  for(filtermap::iterator it=filters.begin(); it!=filters.end(); ++it)
  {
    delete it->second.dataFile;
    for(profmap::const_iterator it2=it->second.profs.begin(); it2!=it->second.profs.end(); ++it2)
      delete[] it2->second.data;
  }
}

int cstats::readinifile(cinput *inputin)
{
  int nerror = 0;
  nerror += inputin->getItem(&swstats   , "stats", "swstats"   , "");
  nerror += inputin->getItem(&sampletime, "stats", "sampletime", "");

  if(!(swstats == "0" || swstats == "1" ))
  {
    ++nerror;
    if(master->mpiid == 0) std::printf("ERROR \"%s\" is an illegal value for swstats\n", swstats.c_str());
  }

  return nerror;
}

int cstats::init(double ifactor)
{
  // convenience pointers for short notation in class
  grid   = model->grid;
  fields = model->fields;
  master = model->master;

  isampletime = (unsigned long)(ifactor * sampletime);

  umodel = new double[grid->kcells];
  vmodel = new double[grid->kcells];

  filtercount = new int[grid->kcells];

  // add the default filter
  filters["default"].name = "default";
  filters["default"].dataFile = NULL;

  filters["wplus"].name = "wplus";
  filters["wplus"].dataFile = NULL;

  filters["wmin"].name = "wmin";
  filters["wmin"].dataFile = NULL;

  filters["ql"].name = "ql";
  filters["ql"].dataFile = NULL;

  filters["qlcore"].name = "qlcore";
  filters["qlcore"].dataFile = NULL;

  // set the number of stats to zero
  nstats = 0;

  return 0;
}

int cstats::create(int n)
{
  // do not create file if stats is disabled
  if(swstats == "0")
    return 0;

  int nerror = 0;

  for(filtermap::iterator it=filters.begin(); it!=filters.end(); ++it)
  {
    // shortcut
    filter *f = &it->second;

    // create a NetCDF file for the statistics
    if(master->mpiid == 0)
    {
      char filename[256];
      std::sprintf(filename, "%s.%s.%07d.nc", master->simname.c_str(), f->name.c_str(), n);
      f->dataFile = new NcFile(filename, NcFile::New);
      if(!f->dataFile->is_valid())
      {
        std::printf("ERROR cannot write statistics file\n");
        ++nerror;
      }
    }
    // crash on all processes in case the file could not be written
    master->broadcast(&nerror, 1);
    if(nerror)
      return 1;

    // create dimensions
    if(master->mpiid == 0)
    {
      f->z_dim  = f->dataFile->add_dim("z" , grid->kmax);
      f->zh_dim = f->dataFile->add_dim("zh", grid->kmax+1);
      f->t_dim  = f->dataFile->add_dim("t");

      NcVar *z_var, *zh_var;

      // create variables belonging to dimensions
      f->iter_var = f->dataFile->add_var("iter", ncInt, f->t_dim);
      f->iter_var->add_att("units", "-");
      f->iter_var->add_att("longname", "Iteration number");

      f->t_var = f->dataFile->add_var("t", ncDouble, f->t_dim);
      f->t_var->add_att("units", "s");
      f->t_var->add_att("longname", "Time");

      z_var = f->dataFile->add_var("z", ncDouble, f->z_dim);
      z_var->add_att("units", "m");
      z_var->add_att("longname", "Full level height");

      zh_var = f->dataFile->add_var("zh", ncDouble, f->zh_dim);
      zh_var->add_att("units", "m");
      zh_var->add_att("longname", "Half level height");

      // save the grid variables
      z_var ->put(&grid->z [grid->kstart], grid->kmax  );
      zh_var->put(&grid->zh[grid->kstart], grid->kmax+1);

      f->dataFile->sync();
    }

  }

  // for each filter add the area as a variable
  addprof("area" , "Fractional area contained in conditional statistics", "-", "z");
  addprof("areah", "Fractional area contained in conditional statistics", "-", "zh");

  return 0;
}

unsigned long cstats::gettimelim(unsigned long itime)
{
  unsigned long idtlim = isampletime -  itime % isampletime;
  return idtlim;
}

int cstats::dostats()
{
  // check if stats are enabled
  if(swstats == "0")
    return 0;

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
  // this function is only called when stats are enabled no need for swstats check

  // check if time for execution
  if(itime % isampletime != 0)
    return 0;

  for(filtermap::iterator it=filters.begin(); it!=filters.end(); ++it)
  {
    // shortcut
    filter *f = &it->second;

    // put the data into the NetCDF file
    if(master->mpiid == 0)
    {
      f->t_var   ->put_rec(&time     , nstats);
      f->iter_var->put_rec(&iteration, nstats);

      for(profmap::const_iterator it=f->profs.begin(); it!=f->profs.end(); ++it)
        f->profs[it->first].ncvar->put_rec(&f->profs[it->first].data[grid->kstart], nstats);

      for(tseriesmap::const_iterator it=f->tseries.begin(); it!=f->tseries.end(); ++it)
        f->tseries[it->first].ncvar->put_rec(&f->tseries[it->first].data, nstats);

      // sync the data
      f->dataFile->sync();
    }
  }

  ++nstats;

  return 0;
}

std::string cstats::getsw()
{
  return swstats;
}

int cstats::addprof(std::string name, std::string longname, std::string unit, std::string zloc)
{
  int nerror = 0;

  // add the profile to all files
  for(filtermap::iterator it=filters.begin(); it!=filters.end(); ++it)
  {
    // shortcut
    filter *f = &it->second;

    // create the NetCDF variable
    if(master->mpiid == 0)
    {
      if(zloc == "z")
      {
        f->profs[name].ncvar = f->dataFile->add_var(name.c_str(), ncDouble, f->t_dim, f->z_dim);
        f->profs[name].data = NULL;
      }
      else if(zloc == "zh")
      {
        f->profs[name].ncvar = f->dataFile->add_var(name.c_str(), ncDouble, f->t_dim, f->zh_dim);
        f->profs[name].data = NULL;
      }
      f->profs[name].ncvar->add_att("units", unit.c_str());
      f->profs[name].ncvar->add_att("long_name", longname.c_str());
      f->profs[name].ncvar->add_att("_FillValue", NC_FILL_DOUBLE);
    }

    // and allocate the memory and initialize at zero
    f->profs[name].data = new double[grid->kcells];
    for(int k=0; k<grid->kcells; ++k)
      f->profs[name].data[k] = 0.;
  }

  return nerror;
}

int cstats::addfixedprof(std::string name, std::string longname, std::string unit, std::string zloc, double * restrict prof)
{
  int nerror = 0;

  // add the profile to all files
  for(filtermap::iterator it=filters.begin(); it!=filters.end(); ++it)
  {
    // shortcut
    filter *f = &it->second;

    // create the NetCDF variable
    NcVar *var;
    if(master->mpiid == 0)
    {
      if(zloc == "z")
        var = f->dataFile->add_var(name.c_str(), ncDouble, f->z_dim);
      else if(zloc == "zh")
        var = f->dataFile->add_var(name.c_str(), ncDouble, f->zh_dim);
      var->add_att("units", unit.c_str());
      var->add_att("long_name", longname.c_str());
      var->add_att("_FillValue", NC_FILL_DOUBLE);

      if(zloc == "z")
        var->put(&prof[grid->kstart], grid->kmax);
      else if(zloc == "zh")
        var->put(&prof[grid->kstart], grid->kmax+1);
    }
  }

  return nerror;
}

int cstats::addtseries(std::string name, std::string longname, std::string unit)
{
  int nerror = 0;

  // add the series to all files
  for(filtermap::iterator it=filters.begin(); it!=filters.end(); ++it)
  {
    // shortcut
    filter *f = &it->second;

    // create the NetCDF variable
    if(master->mpiid == 0)
    {
      f->tseries[name].ncvar = f->dataFile->add_var(name.c_str(), ncDouble, f->t_dim);
      f->tseries[name].ncvar->add_att("units", unit.c_str());
      f->tseries[name].ncvar->add_att("long_name", longname.c_str());
      f->tseries[name].ncvar->add_att("_FillValue", NC_FILL_DOUBLE);
    }

    // and initialize at zero
    f->tseries[name].data = 0.;
  }

  return nerror;
}

int cstats::getfilter(cfield3d *ffield, filter *f)
{
  calcfilter(ffield->data, f->profs["area"].data, f->profs["areah"].data, filtercount);
  return 0;
}

// COMPUTATIONAL KERNELS BELOW
int cstats::calcfilter(double * restrict fdata, double * restrict area, double * restrict areah, int * restrict nfilter)
{
  int ijtot = grid->itot*grid->jtot;

  // set all the filter values to 1
  for(int n=0; n<grid->ncells; ++n)
    fdata[n] = 1.;

  for(int k=0; k<grid->kcells; ++k)
  {
    nfilter[k] = ijtot;
    area[k] = 1.;
    areah[k] = 1.;
  }

  return 0;
}


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

int cstats::calcmean(double * restrict data, double * restrict prof, double offset, const int loc[3],
                     double * restrict filter, int * restrict nfilter)
{
  int ijk,ii,jj,kk;
  double filterval;

  ii = 1;
  jj = grid->icells;
  kk = grid->ijcells;

  // interpolation offset, if locz = 1, which corresponds to half level, there is no interpolation
  int iif = loc[0]*ii;
  int jjf = loc[1]*jj;
  int kkf = loc[2]*kk;
  
  for(int k=1; k<grid->kcells; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        filterval = (1./6.)*(filter[ijk-iif] + filter[ijk-jjf] + filter[ijk-kkf] + 3.*filter[ijk]);
        prof[k] += filterval*(data[ijk] + offset);
      }
  }

  master->sum(prof, grid->kcells);

  // use the same interpolation trick as for the filter field, no interpolation on half levels
  int kf = loc[2];
  for(int k=1; k<grid->kcells; k++)
  {
    if(nfilter[k-kf]+nfilter[k] > 0)
      prof[k] /= (0.5*(double)(nfilter[k-kf] + nfilter[k]));
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}

/*
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
*/

// \TODO the count function assumes that the variable to count is at the filter location
int cstats::calccount(double * restrict data, double * restrict prof, double threshold,
                      double * restrict filter, int * restrict nfilter)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->ijcells;

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
          prof[k] += filter[ijk]*1.;
        }
      }
  }

  master->sum(prof, grid->kcells);

  for(int k=0; k<grid->kcells; k++)
  {
    if(nfilter[k] > 0)
      prof[k] /= (double)(nfilter[k]);
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}

/*
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
*/

int cstats::calcmoment(double * restrict data, double * restrict datamean, double * restrict prof, double power, const int loc[3],
                       double * restrict filter, int * restrict nfilter)
{
  int ijk,ii,jj,kk;
  double filterval;

  ii = 1;
  jj = grid->icells;
  kk = grid->ijcells;

  // interpolation offset, if locz = 1, which corresponds to half level, there is no interpolation
  int iif = loc[0]*ii;
  int jjf = loc[1]*jj;
  int kkf = loc[2]*kk;
 
  for(int k=grid->kstart; k<grid->kend+1; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj + k*kk;
        filterval = (1./6.)*(filter[ijk-iif] + filter[ijk-jjf] + filter[ijk-kkf] + 3.*filter[ijk]);
        prof[k] += filterval*std::pow(data[ijk]-datamean[k], power);
      }
  }

  master->sum(prof, grid->kcells);

  // use the same interpolation trick as for the filter field, no interpolation on half levels
  int kf = loc[2];

  double n = grid->itot*grid->jtot;
  for(int k=1; k<grid->kcells; k++)
  {
    if(nfilter[k-kf]+nfilter[k] > 0)
      prof[k] /= (0.5*(double)(nfilter[k-kf] + nfilter[k]));
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}

/*
int cstats::calcflux_2nd(double * restrict data, double * restrict w, double * restrict prof, double * restrict tmp1, int locx, int locy)
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
*/

int cstats::calcflux_2nd(double * restrict data, double * restrict datamean, double * restrict w, double * restrict wmean,
                         double * restrict prof, double * restrict tmp1, const int loc[3],
                         double * restrict filter, int * restrict nfilter)
{
  int ijk,ii,jj,kk;
  double filterval;

  ii = 1;
  jj = grid->icells;
  kk = grid->ijcells;

  // interpolation offset, if locz = 1, which corresponds to half level, there is no interpolation
  int iif = loc[0]*ii;
  int jjf = loc[1]*jj;
  int kkf = loc[2]*kk;
 
  // set a pointer to the field that contains w, either interpolated or the original
  double * restrict calcw = w;
  if(loc[0] == 1)
  {
    grid->interpolatex_2nd(tmp1, w, 0);
    calcw = tmp1;
  }
  else if(loc[1] == 1)
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
        filterval = (1./6.)*(filter[ijk-iif] + filter[ijk-jjf] + filter[ijk-kkf] + 3.*filter[ijk]);
        prof[k] += filterval*(0.5*(data[ijk-kk]+data[ijk])-datamean[k])*(calcw[ijk]-wmean[k]);
      }
  }

  master->sum(prof, grid->kcells);

  // use the same interpolation trick as for the filter field, no interpolation on half levels
  int kf = loc[2];

  double n = grid->itot*grid->jtot;
  for(int k=1; k<grid->kcells; k++)
  {
    if(nfilter[k-kf]+nfilter[k] > 0)
      prof[k] /= (0.5*(double)(nfilter[k-kf] + nfilter[k]));
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}

int cstats::calcflux_4th(double * restrict data, double * restrict w, double * restrict prof, double * restrict tmp1, int locx, int locy)
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
  
  for(int k=grid->kstart; k<grid->kend+1; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj + k*kk1;
        prof[k] += (ci0*data[ijk-kk2] + ci1*data[ijk-kk1] + ci2*data[ijk] + ci3*data[ijk+kk1])*calcw[ijk];
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+1; ++k)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

/*
int cstats::calcgrad_2nd(double * restrict data, double * restrict prof, double * restrict dzhi)
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
*/

int cstats::calcgrad_2nd(double * restrict data, double * restrict prof, double * restrict dzhi, const int loc[3],
                         double * restrict filter, int * restrict nfilter)
{
  int ijk,ii,jj,kk;
  double filterval;

  ii = 1;
  jj = grid->icells;
  kk = grid->ijcells;

  // interpolation offset, if locz = 1, which corresponds to half level, there is no interpolation
  int iif = loc[0]*ii;
  int jjf = loc[1]*jj;
  int kkf = loc[2]*kk;
 
  for(int k=grid->kstart; k<grid->kend+1; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj + k*kk;
        filterval = (1./6.)*(filter[ijk-iif] + filter[ijk-jjf] + filter[ijk-kkf] + 3.*filter[ijk]);
        prof[k] += filterval*(data[ijk]-data[ijk-kk])*dzhi[k];
      }
  }

  master->sum(prof, grid->kcells);

  // use the same interpolation trick as for the filter field, no interpolation on half levels
  int kf = loc[2];

  double n = grid->itot*grid->jtot;
  for(int k=1; k<grid->kcells; k++)
  {
    if(nfilter[k-kf]+nfilter[k] > 0)
      prof[k] /= (0.5*(double)(nfilter[k-kf] + nfilter[k]));
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}

int cstats::calcgrad_4th(double * restrict data, double * restrict prof, double * restrict dzhi4, const int loc[3],
                         double * restrict filter, int * restrict nfilter)
{
  int ijk,ii,jj,kk1,kk2;
  double filterval;

  ii  = 1;
  jj  = 1*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  // interpolation offset, if locz = 1, which corresponds to half level, there is no interpolation
  int iif = loc[0]*ii;
  int jjf = loc[1]*jj;
  int kkf = loc[2]*kk1;
 
  for(int k=grid->kstart; k<grid->kend+1; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj + k*kk1;
        filterval = (1./6.)*(filter[ijk-iif] + filter[ijk-jjf] + filter[ijk-kkf] + 3.*filter[ijk]);
        prof[k] += filterval*(cg0*data[ijk-kk2] + cg1*data[ijk-kk1] + cg2*data[ijk] + cg3*data[ijk+kk1])*dzhi4[k];
      }
  }

  master->sum(prof, grid->kcells);

  // use the same interpolation trick as for the filter field, no interpolation on half levels
  int kf = loc[2];

  double n = grid->itot*grid->jtot;
  for(int k=1; k<grid->kcells; k++)
  {
    if(nfilter[k-kf]+nfilter[k] > 0)
      prof[k] /= (0.5*(double)(nfilter[k-kf] + nfilter[k]));
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}

int cstats::calcdiff_4th(double * restrict data, double * restrict prof, double * restrict dzhi4, double visc, const int loc[3],
                         double * restrict filter, int * restrict nfilter)
{
  int ijk,ii,jj,kk1,kk2;
  double filterval;

  ii  = 1;
  jj  = 1*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  // interpolation offset, if locz = 1, which corresponds to half level, there is no interpolation
  int iif = loc[0]*ii;
  int jjf = loc[1]*jj;
  int kkf = loc[2]*kk1;
 
  for(int k=grid->kstart; k<grid->kend+1; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj + k*kk1;
        filterval = (1./6.)*(filter[ijk-iif] + filter[ijk-jjf] + filter[ijk-kkf] + 3.*filter[ijk]);
        prof[k] -= filterval*visc*(cg0*data[ijk-kk2] + cg1*data[ijk-kk1] + cg2*data[ijk] + cg3*data[ijk+kk1])*dzhi4[k];
      }
  }

  // use the same interpolation trick as for the filter field, no interpolation on half levels
  int kf = loc[2];

  double n = grid->itot*grid->jtot;
  for(int k=1; k<grid->kcells; k++)
  {
    if(nfilter[k-kf]+nfilter[k] > 0)
      prof[k] /= (0.5*(double)(nfilter[k-kf] + nfilter[k]));
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}

int cstats::calcdiff_2nd(double * restrict data, double * restrict evisc, double * restrict prof, double * restrict dzhi, double * restrict fluxbot, double * restrict fluxtop, double tPr)
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

int cstats::calcpath(double * restrict data, double * restrict path)
{
  int ijk,jj,kk;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  int kstart = grid->kstart;
  int nerror = 0;


  *path = 0.;
  // Integrate with height
  for(int k=kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        *path += fields->rhoref[k] * data[ijk] * grid->dz[k];
      }

  *path /= 1.0*grid->imax*grid->jmax;

  grid->getprof(path,1);

  return nerror;
}


int cstats::calccover(double * restrict data, double * restrict cover, double threshold)
{
  int ijk,jj,kk;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  int kstart = grid->kstart;
  int nerror = 0;

  *cover = 0.;
  // Integrate with height
  for(int j=grid->jstart; j<grid->jend; j++)
    for(int i=grid->istart; i<grid->iend; i++)
      for(int k=kstart; k<grid->kend; k++)
      {
        ijk  = i + j*jj + k*kk;
        if (data[ijk]>threshold)
        {
          *cover += 1.;
          break;
        }
      }

  *cover /= grid->imax*grid->jmax;

  grid->getprof(cover,1);

  return nerror;
}

