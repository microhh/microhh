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
#include <cmath>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "thermo_moist.h"
#include "defines.h"
#include "constants.h"
#include "fd.h"
#include "model.h"
#include "diff_smag2.h"
#include "timeloop.h"
#include <netcdfcpp.h>

#define NTHRES 0

Stats::Stats(Model *modelin, Input *inputin)
{
  model = modelin;

  // set the pointers to zero
  nmask  = 0;
  nmaskh = 0;

  int nerror = 0;
  nerror += inputin->getItem(&swstats, "stats", "swstats", "", "0");

  if(swstats == "1")
    nerror += inputin->getItem(&sampletime, "stats", "sampletime", "");

  if(!(swstats == "0" || swstats == "1" ))
  {
    ++nerror;
    master->printError("\"%s\" is an illegal value for swstats\n", swstats.c_str());
  }

  if(nerror)
    throw 1;
}

Stats::~Stats()
{
  delete[] nmask;
  delete[] nmaskh;

  // delete the profiles
  for(MaskMap::iterator it=masks.begin(); it!=masks.end(); ++it)
  {
    delete it->second.dataFile;
    for(profmap::const_iterator it2=it->second.profs.begin(); it2!=it->second.profs.end(); ++it2)
      delete[] it2->second.data;
  }
}

void Stats::init(double ifactor)
{
  // convenience pointers for short notation in class
  grid   = model->grid;
  fields = model->fields;
  master = model->master;

  // add the default mask
  addmask("default");

  isampletime = (unsigned long)(ifactor * sampletime);

  nmask  = new int[grid->kcells];
  nmaskh = new int[grid->kcells];

  // set the number of stats to zero
  nstats = 0;
}

void Stats::create(int n)
{
  // do not create file if stats is disabled
  if(swstats == "0")
    return;

  int nerror = 0;

  for(MaskMap::iterator it=masks.begin(); it!=masks.end(); ++it)
  {
    // shortcut
    Mask *m = &it->second;

    // create a NetCDF file for the statistics
    if(master->mpiid == 0)
    {
      char filename[256];
      std::sprintf(filename, "%s.%s.%07d.nc", master->simname.c_str(), m->name.c_str(), n);
      m->dataFile = new NcFile(filename, NcFile::New);
      if(!m->dataFile->is_valid())
      {
        master->printError("cannot write statistics file\n");
        ++nerror;
      }
    }
    // crash on all processes in case the file could not be written
    master->broadcast(&nerror, 1);
    if(nerror)
      throw 1;

    // create dimensions
    if(master->mpiid == 0)
    {
      m->z_dim  = m->dataFile->add_dim("z" , grid->kmax);
      m->zh_dim = m->dataFile->add_dim("zh", grid->kmax+1);
      m->t_dim  = m->dataFile->add_dim("t");

      NcVar *z_var, *zh_var;

      // create variables belonging to dimensions
      m->iter_var = m->dataFile->add_var("iter", ncInt, m->t_dim);
      m->iter_var->add_att("units", "-");
      m->iter_var->add_att("long_name", "Iteration number");

      m->t_var = m->dataFile->add_var("t", ncDouble, m->t_dim);
      m->t_var->add_att("units", "s");
      m->t_var->add_att("long_name", "Time");

      z_var = m->dataFile->add_var("z", ncDouble, m->z_dim);
      z_var->add_att("units", "m");
      z_var->add_att("long_name", "Full level height");

      zh_var = m->dataFile->add_var("zh", ncDouble, m->zh_dim);
      zh_var->add_att("units", "m");
      zh_var->add_att("long_name", "Half level height");

      // save the grid variables
      z_var ->put(&grid->z [grid->kstart], grid->kmax  );
      zh_var->put(&grid->zh[grid->kstart], grid->kmax+1);

      m->dataFile->sync();
    }

  }

  // for each mask add the area as a variable
  addprof("area" , "Fractional area contained in mask", "-", "z");
  addprof("areah", "Fractional area contained in mask", "-", "zh");
}

unsigned long Stats::getTimeLimit(unsigned long itime)
{
  if(swstats == "0")
    return constants::ulhuge;

  unsigned long idtlim = isampletime - itime % isampletime;
  return idtlim;
}

int Stats::dostats()
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

int Stats::exec(int iteration, double time, unsigned long itime)
{
  // this function is only called when stats are enabled no need for swstats check

  // check if time for execution
  if(itime % isampletime != 0)
    return 0;

  for(MaskMap::iterator it=masks.begin(); it!=masks.end(); ++it)
  {
    // shortcut
    Mask *m = &it->second;

    // put the data into the NetCDF file
    if(master->mpiid == 0)
    {
      m->t_var   ->put_rec(&time     , nstats);
      m->iter_var->put_rec(&iteration, nstats);

      for(profmap::const_iterator it=m->profs.begin(); it!=m->profs.end(); ++it)
        m->profs[it->first].ncvar->put_rec(&m->profs[it->first].data[grid->kstart], nstats);

      for(tseriesmap::const_iterator it=m->tseries.begin(); it!=m->tseries.end(); ++it)
        m->tseries[it->first].ncvar->put_rec(&m->tseries[it->first].data, nstats);

      // sync the data
      m->dataFile->sync();
    }
  }

  ++nstats;

  return 0;
}

std::string Stats::getSwitch()
{
  return swstats;
}

void Stats::addmask(const std::string maskname)
{
  masks[maskname].name = maskname;
  masks[maskname].dataFile = 0;
}

int Stats::addprof(std::string name, std::string longname, std::string unit, std::string zloc)
{
  int nerror = 0;

  // add the profile to all files
  for(MaskMap::iterator it=masks.begin(); it!=masks.end(); ++it)
  {
    // shortcut
    Mask *m = &it->second;

    // create the NetCDF variable
    if(master->mpiid == 0)
    {
      if(zloc == "z")
      {
        m->profs[name].ncvar = m->dataFile->add_var(name.c_str(), ncDouble, m->t_dim, m->z_dim);
        m->profs[name].data = NULL;
      }
      else if(zloc == "zh")
      {
        m->profs[name].ncvar = m->dataFile->add_var(name.c_str(), ncDouble, m->t_dim, m->zh_dim);
        m->profs[name].data = NULL;
      }
      m->profs[name].ncvar->add_att("units", unit.c_str());
      m->profs[name].ncvar->add_att("long_name", longname.c_str());
      m->profs[name].ncvar->add_att("_FillValue", NC_FILL_DOUBLE);
    }

    // and allocate the memory and initialize at zero
    m->profs[name].data = new double[grid->kcells];
    for(int k=0; k<grid->kcells; ++k)
      m->profs[name].data[k] = 0.;
  }

  return nerror;
}

int Stats::addfixedprof(std::string name, std::string longname, std::string unit, std::string zloc, double * restrict prof)
{
  int nerror = 0;

  // add the profile to all files
  for(MaskMap::iterator it=masks.begin(); it!=masks.end(); ++it)
  {
    // shortcut
    Mask *m = &it->second;

    // create the NetCDF variable
    NcVar *var;
    if(master->mpiid == 0)
    {
      if(zloc == "z")
        var = m->dataFile->add_var(name.c_str(), ncDouble, m->z_dim);
      else if(zloc == "zh")
        var = m->dataFile->add_var(name.c_str(), ncDouble, m->zh_dim);
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

int Stats::addtseries(std::string name, std::string longname, std::string unit)
{
  int nerror = 0;

  // add the series to all files
  for(MaskMap::iterator it=masks.begin(); it!=masks.end(); ++it)
  {
    // shortcut
    Mask *m = &it->second;

    // create the NetCDF variable
    if(master->mpiid == 0)
    {
      m->tseries[name].ncvar = m->dataFile->add_var(name.c_str(), ncDouble, m->t_dim);
      m->tseries[name].ncvar->add_att("units", unit.c_str());
      m->tseries[name].ncvar->add_att("long_name", longname.c_str());
      m->tseries[name].ncvar->add_att("_FillValue", NC_FILL_DOUBLE);
    }

    // and initialize at zero
    m->tseries[name].data = 0.;
  }

  return nerror;
}

int Stats::getMask(Field3d *mfield, Field3d *mfieldh, Mask *m)
{
  calcmask(mfield->data, mfieldh->data, mfieldh->databot,
           nmask, nmaskh, &nmaskbot);
  return 0;
}

// COMPUTATIONAL KERNELS BELOW
int Stats::calcmask(double * restrict mask, double * restrict maskh, double * restrict maskbot,
                     int * restrict nmask, int * restrict nmaskh, int * restrict nmaskbot)
{
  int ijtot = grid->itot*grid->jtot;

  // set all the mask values to 1
  for(int n=0; n<grid->ncells; ++n)
    mask[n] = 1.;

  for(int n=0; n<grid->ncells; ++n)
    maskh[n] = 1.;

  for(int n=0; n<grid->ijcells; ++n)
    maskbot[n] = 1.;

  for(int k=0; k<grid->kcells; ++k)
  {
    nmask [k] = ijtot;
    nmaskh[k] = ijtot;
  }
  *nmaskbot = ijtot;

  return 0;
}

/*
void Stats::calcmean(double * const restrict prof, const double * const restrict data,
                      const double offset)
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
        prof[k] += data[ijk] + offset;
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=0; k<grid->kcells; ++k)
    prof[k] /= n;

  grid->getProf(prof, grid->kcells);
}
*/

int Stats::calcarea(double * restrict area, const int loc[3], int * restrict nmask)
{
  int ijtot = grid->itot*grid->jtot;

  for(int k=grid->kstart; k<grid->kend+loc[2]; k++)
  {
    if(nmask[k] > NTHRES)
      area[k] = (double)(nmask[k]) / (double)ijtot;
    else
      area[k] = 0.;
  }

  return 0;
}

void Stats::calcmean(double * const restrict prof, const double * const restrict data,
                      const double offset, const int loc[3],
                      const double * const restrict mask, const int * const restrict nmask)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->ijcells;

  for(int k=1; k<grid->kcells; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += mask[ijk]*(data[ijk] + offset);
      }
  }

  master->sum(prof, grid->kcells);

  for(int k=1; k<grid->kcells; k++)
  {
    if(nmask[k] > NTHRES)
      prof[k] /= (double)(nmask[k]);
    else
      prof[k] = NC_FILL_DOUBLE;
  }
}

void Stats::calcmean2d(double * const restrict mean, const double * const restrict data,
                        const double offset,
                        const double * const restrict mask, const int * const restrict nmask)
{
  int ij,jj;
  jj = grid->icells;

  if(*nmask > NTHRES)
  {
    *mean = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ij  = i + j*jj;
        *mean += mask[ij]*(data[ij] + offset);
      }
    master->sum(mean,1);
    *mean /= (double)*nmask;
  }
  else
    *mean = NC_FILL_DOUBLE; 
}

int Stats::calcsortprof(double * restrict data, double * restrict bin, double * restrict prof)
{
  int ijk,jj,kk,index,kstart,kend;
  double minval,maxval,range;

  jj = grid->icells;
  kk = grid->ijcells;
  kstart = grid->kstart;
  kend = grid->kend;

  minval =  constants::dhuge;
  maxval = -constants::dhuge;

  // first, get min and max
  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        if(data[ijk] < minval)
          minval = data[ijk];
        if(data[ijk] > maxval)
          maxval = data[ijk];
      }

  master->min(&minval, 1);
  master->max(&maxval, 1);

  // make sure that the max ends up in the last bin (introduce 1E-9 error)
  maxval *= (1.+constants::dsmall);

  range = maxval-minval;

  // create bins, equal to the number of grid cells per proc
  // make sure that bins is not larger than the memory of one 3d field
  int bins = grid->nmax;

  // calculate bin width, subtract one to make the minimum and maximum 
  // are in the middle of the bin range and add half a bin size on both sides
  // |----x----|----x----|----x----|
  double dbin = range / (double)(bins-1);
  minval -= 0.5*dbin;
  maxval += 0.5*dbin;

  // set the bin array to zero
  for(int n=0; n<bins; ++n)
    bin[n] = 0;

  // calculate the division factor of one equivalent height unit
  // (the total volume saved is itot*jtot*zsize)
  double nslice = (double)(grid->itot*grid->jtot);

  // check in which bin each value falls and increment the bin count
  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
      // do not add a ivdep pragma here, because multiple instances could write the same bin[index]
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        index = (int)((data[ijk] - minval) / dbin);
        bin[index] += grid->dz[k] / nslice;
      }

  // get the bin count
  master->sum(bin, bins);

  // set the starting values of the loop
  index = 0;
  double zbin = 0.5*bin[index];
  double profval = minval + 0.5*dbin;
  double dzfrac;
  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    // Integrate the profile up to the bin count.
    // Escape the while loop when the integrated profile 
    // exceeds the next grid point.
    while(zbin < grid->z[k])
    {
      zbin += 0.5*(bin[index]+bin[index+1]);
      profval += dbin;
      ++index;
    }

    dzfrac = (zbin-grid->z[k]) / (0.5*(bin[index-1]+bin[index]));

    prof[k] = profval - dzfrac*dbin;
  }

  // now calculate the ghost cells
  // \TODO this is not accurate enough, extrapolate properly
  double profbot = minval;
  double proftop = maxval;

  if(grid->swspatialorder == "2")
  {
    prof[kstart-1] = 2.*profbot - prof[kstart];
    prof[kend]     = 2.*proftop - prof[kend-1];
  }
  else if(grid->swspatialorder == "4")
  {
    prof[kstart-1] = (8./3.)*profbot - 2.*prof[kstart] + (1./3.)*prof[kstart+1];
    prof[kstart-2] = 8.*profbot      - 9.*prof[kstart] + 2.*prof[kstart+1];
    prof[kend]     = (8./3.)*proftop - 2.*prof[kend-1] + (1./3.)*prof[kend-2];
    prof[kend+1]   = 8.*proftop      - 9.*prof[kend-1] + 2.*prof[kend-2];
  }

  return 0;
}

/*
int Stats::calccount(double * restrict data, double * restrict prof, double threshold)
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
          prof[k] += 1.;
        }
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=0; k<grid->kcells; ++k)
    prof[k] /= n;

  grid->getProf(prof, grid->kcells);

  return 0;
}
*/

// \TODO the count function assumes that the variable to count is at the mask location
int Stats::calccount(double * restrict data, double * restrict prof, double threshold,
                      double * restrict mask, int * restrict nmask)
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
        if(data[ijk] > threshold)
          prof[k] += mask[ijk]*1.;
      }
  }

  master->sum(prof, grid->kcells);

  for(int k=0; k<grid->kcells; k++)
  {
    if(nmask[k] > NTHRES)
      prof[k] /= (double)(nmask[k]);
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}

/*
int Stats::calcmoment(double * restrict data, double * restrict datamean, double * restrict prof, double power, int a)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->ijcells;
  
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

  grid->getProf(prof, grid->kcells);

  return 0;
}
*/

int Stats::calcmoment(double * restrict data, double * restrict datamean, double * restrict prof, double power, const int loc[3],
                       double * restrict mask, int * restrict nmask)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->ijcells;
 
  for(int k=grid->kstart; k<grid->kend+1; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        prof[k] += mask[ijk]*std::pow(data[ijk]-datamean[k], power);
      }
  }

  master->sum(prof, grid->kcells);

  for(int k=1; k<grid->kcells; k++)
  {
    if(nmask[k] > NTHRES)
      prof[k] /= (double)(nmask[k]);
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}

/*
int Stats::calcflux_2nd(double * restrict data, double * restrict w, double * restrict prof, double * restrict tmp1, int locx, int locy)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->ijcells;

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

  grid->getProf(prof, grid->kcells);

  return 0;
}
*/

int Stats::calcflux_2nd(double * restrict data, double * restrict datamean, double * restrict w, double * restrict wmean,
                        double * restrict prof, double * restrict tmp1, const int loc[3],
                        double * restrict mask, int * restrict nmask)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->ijcells;
 
  // set a pointer to the field that contains w, either interpolated or the original
  double * restrict calcw = w;

  // define the locations
  const int wloc [3] = {0,0,1};
  const int uwloc[3] = {1,0,1};
  const int vwloc[3] = {0,1,1};

  if(loc[0] == 1)
  {
    grid->interpolate_2nd(tmp1, w, wloc, uwloc);
    calcw = tmp1;
  }
  else if(loc[1] == 1)
  {
    grid->interpolate_2nd(tmp1, w, wloc, vwloc);
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
        prof[k] += mask[ijk]*(0.5*(data[ijk-kk]+data[ijk])-0.5*(datamean[k-1]+datamean[k]))*(calcw[ijk]-wmean[k]);
        // prof[k] += mask[ijk]*0.5*(data[ijk-kk]+data[ijk])*calcw[ijk];
      }
  }

  master->sum(prof, grid->kcells);

  for(int k=1; k<grid->kcells; ++k)
  {
    if(nmask[k] > NTHRES && datamean[k-1] != NC_FILL_DOUBLE && datamean[k] != NC_FILL_DOUBLE)
      prof[k] /= (double)(nmask[k]);
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}

int Stats::calcflux_4th(double * restrict data, double * restrict w, double * restrict prof, double * restrict tmp1, const int loc[3],
                         double * restrict mask, int * restrict nmask)
{
  using namespace fd::o4;

  int ijk,jj,kk1,kk2;

  jj  = 1*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  // set a pointer to the field that contains w, either interpolated or the original
  double * restrict calcw = w;

  // define the locations
  const int wloc [3] = {0,0,1};
  const int uwloc[3] = {1,0,1};
  const int vwloc[3] = {0,1,1};

  if(loc[0] == 1)
  {
    grid->interpolate_4th(tmp1, w, wloc, uwloc);
    calcw = tmp1;
  }
  else if(loc[1] == 1)
  {
    grid->interpolate_4th(tmp1, w, wloc, vwloc);
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
        prof[k] += mask[ijk]*(ci0*data[ijk-kk2] + ci1*data[ijk-kk1] + ci2*data[ijk] + ci3*data[ijk+kk1])*calcw[ijk];
      }
  }

  master->sum(prof, grid->kcells);

  for(int k=1; k<grid->kcells; k++)
  {
    if(nmask[k] > NTHRES)
      prof[k] /= (double)(nmask[k]);
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}

/*
int Stats::calcgrad_2nd(double * restrict data, double * restrict prof, double * restrict dzhi)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->ijcells;
  
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

  grid->getProf(prof, grid->kcells);

  return 0;
}
*/

int Stats::calcgrad_2nd(double * restrict data, double * restrict prof, double * restrict dzhi, const int loc[3],
                         double * restrict mask, int * restrict nmask)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->ijcells;

  for(int k=grid->kstart; k<grid->kend+1; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += mask[ijk]*(data[ijk]-data[ijk-kk])*dzhi[k];
      }
  }

  master->sum(prof, grid->kcells);

  for(int k=1; k<grid->kcells; k++)
  {
    if(nmask[k] > NTHRES)
      prof[k] /= (double)(nmask[k]);
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}

int Stats::calcgrad_4th(double * restrict data, double * restrict prof, double * restrict dzhi4, const int loc[3],
                         double * restrict mask, int * restrict nmask)
{
  using namespace fd::o4;

  int ijk,jj,kk1,kk2;

  jj  = 1*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  for(int k=grid->kstart; k<grid->kend+1; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj + k*kk1;
        prof[k] += mask[ijk]*(cg0*data[ijk-kk2] + cg1*data[ijk-kk1] + cg2*data[ijk] + cg3*data[ijk+kk1])*dzhi4[k];
      }
  }

  master->sum(prof, grid->kcells);

  for(int k=1; k<grid->kcells; k++)
  {
    if(nmask[k] > NTHRES)
      prof[k] /= (double)(nmask[k]);
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}

int Stats::calcdiff_4th(double * restrict data, double * restrict prof, double * restrict dzhi4, double visc, const int loc[3],
                         double * restrict mask, int * restrict nmask)
{
  using namespace fd::o4;

  int ijk,jj,kk1,kk2;

  jj  = 1*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;
 
  for(int k=grid->kstart; k<grid->kend+1; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj + k*kk1;
        prof[k] -= mask[ijk]*visc*(cg0*data[ijk-kk2] + cg1*data[ijk-kk1] + cg2*data[ijk] + cg3*data[ijk+kk1])*dzhi4[k];
      }
  }

  master->sum(prof, grid->kcells);

  for(int k=1; k<grid->kcells; k++)
  {
    if(nmask[k] > NTHRES)
      prof[k] /= (double)(nmask[k]);
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}

int Stats::calcdiff_2nd(double * restrict data, double * restrict prof, double * restrict dzhi, double visc, const int loc[3],
                         double * restrict mask, int * restrict nmask)
{
  int ijk,jj,kk;

  jj = 1*grid->icells;
  kk = 1*grid->ijcells;
 
  for(int k=grid->kstart; k<grid->kend+1; ++k)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] -= mask[ijk]*visc*(data[ijk] - data[ijk-kk])*dzhi[k];
      }
  }

  master->sum(prof, grid->kcells);

  for(int k=1; k<grid->kcells; k++)
  {
    if(nmask[k] > NTHRES)
      prof[k] /= (double)(nmask[k]);
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}


int Stats::calcdiff_2nd(double * restrict data, double * restrict w, double * restrict evisc,
                         double * restrict prof, double * restrict dzhi,
                         double * restrict fluxbot, double * restrict fluxtop, double tPr, const int loc[3],
                         double * restrict mask, int * restrict nmask)
{
  int ijk,ij,ii,jj,kk,kstart,kend;
  double eviscu,eviscv,eviscs;
  double dxi,dyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->ijcells;
  kstart = grid->kstart;
  kend = grid->kend;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  // bottom boundary
  prof[kstart] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      prof[kstart] += mask[ijk]*fluxbot[ij];
    }

  // calculate the interior
  if(loc[0] == 1)
  {
    for(int k=grid->kstart+1; k<grid->kend; ++k)
    {
      prof[k] = 0.;
      for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for(int i=grid->istart; i<grid->iend; ++i)
        {
          ijk  = i + j*jj + k*kk;
          // evisc * (du/dz + dw/dx)
          eviscu = 0.25*(evisc[ijk-ii-kk]+evisc[ijk-ii]+evisc[ijk-kk]+evisc[ijk]);
          prof[k] += -mask[ijk]*eviscu*( (data[ijk]-data[ijk-kk])*dzhi[k] + (w[ijk]-w[ijk-ii])*dxi );
        }
    }
  }
  else if(loc[1] == 1)
  {
    for(int k=grid->kstart+1; k<grid->kend; ++k)
    {
      prof[k] = 0.;
      for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for(int i=grid->istart; i<grid->iend; ++i)
        {
          ijk  = i + j*jj + k*kk;
          // evisc * (dv/dz + dw/dy)
          eviscv = 0.25*(evisc[ijk-jj-kk]+evisc[ijk-jj]+evisc[ijk-kk]+evisc[ijk]);
          prof[k] += -mask[ijk]*eviscv*( (data[ijk]-data[ijk-kk])*dzhi[k] + (w[ijk]-w[ijk-jj])*dyi );
        }
    }
  }
  else
  {
    for(int k=grid->kstart+1; k<grid->kend; ++k)
    {
      prof[k] = 0.;
      for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for(int i=grid->istart; i<grid->iend; ++i)
        {
          ijk  = i + j*jj + k*kk;
          eviscs = 0.5*(evisc[ijk-kk]+evisc[ijk])/tPr;
          prof[k] += -mask[ijk]*eviscs*(data[ijk]-data[ijk-kk])*dzhi[k];
        }
    }
  }

  // top boundary
  prof[kend] = 0.;
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kend*kk;
      prof[kend] += mask[ijk]*fluxtop[ij];
    }

  master->sum(prof, grid->kcells);

  for(int k=1; k<grid->kcells; k++)
  {
    if(nmask[k] > NTHRES)
      prof[k] /= (double)(nmask[k]);
    else
      prof[k] = NC_FILL_DOUBLE;
  }

  return 0;
}

int Stats::addfluxes(double * restrict flux, double * restrict turb, double * restrict diff)
{
  for(int k=grid->kstart; k<grid->kend+1; ++k)
  {
    if(turb[k] == NC_FILL_DOUBLE || diff[k] == NC_FILL_DOUBLE)
      flux[k] = NC_FILL_DOUBLE;
    else
      flux[k] = turb[k] + diff[k];
  }

  return 0;
}

/**
 * This function calculates the total domain integrated path of variable data over maskbot
 */
int Stats::calcpath(double * restrict data, double * restrict maskbot, int * restrict nmaskbot, double * restrict path)
{
  int ijk,ij,jj,kk;
  jj = grid->icells;
  kk = grid->ijcells;
  int kstart = grid->kstart;

  *path = 0.;

  if(*nmaskbot > NTHRES)
  {
    // Integrate liquid water
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ij  = i + j*jj;
        if(maskbot[ij] == 1)
          for(int k=kstart; k<grid->kend; k++)
          {
            ijk = i + j*jj + k*kk;
            *path += fields->rhoref[k] * data[ijk] * grid->dz[k];
          }
      }
    *path /= (double)*nmaskbot;
    grid->getProf(path,1);
  }
  else
    *path = NC_FILL_DOUBLE;

  return 0;
}

/**
 * This function calculates the vertical projected cover of variable data over maskbot
 */
int Stats::calccover(double * restrict data, double * restrict maskbot, int * restrict nmaskbot, double * restrict cover, double threshold)
{
  int ijk,ij,jj,kk;
  jj = grid->icells;
  kk = grid->ijcells;
  int kstart = grid->kstart;

  *cover = 0.;

  if(*nmaskbot > NTHRES)
  {
    // Per column, check if cloud present
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ij  = i + j*jj;
        if(maskbot[ij] == 1)
          for(int k=kstart; k<grid->kend; k++)
          {
            ijk = i + j*jj + k*kk;
            if(data[ijk]>threshold)
            {
              *cover += 1.;
              break;
            }
          }
      }
    *cover /= (double)*nmaskbot;
    grid->getProf(cover,1);
  }
  else
    *cover = NC_FILL_DOUBLE;

  return 0;
}

