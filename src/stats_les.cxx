#include <cstdio>
#include <cmath>
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "defines.h"
#include <netcdfcpp.h>

#define NO_OFFSET 0.

cstats_les::cstats_les(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;

  allocated   = false;
  initialized = false;
}

cstats_les::~cstats_les()
{
  if(initialized)
    delete dataFile;

  if(allocated)
  {
    delete[] uabs;
    delete[] vabs;

    // delete the profiles
    for(profmap::const_iterator it=profs.begin(); it!=profs.end(); ++it)
      delete[] it->second.data;
  }
}

int cstats_les::init()
{
  uabs = new double[grid->kcells];
  vabs = new double[grid->kcells];

  allocated = true;

  // set the number of stats to zero
  nstats = 0;

  return 0;
}

int cstats_les::create(int n)
{
  int nerror = 0;

  // create a NetCDF file for the statistics
  if(mpi->mpiid == 0)
  {
    char filename[256];
    std::sprintf(filename, "%s.%07d.nc", mpi->simname.c_str(), n);
    dataFile = new NcFile(filename, NcFile::New);
    if(!dataFile->is_valid())
    {
      std::printf("ERROR cannot write statistics file\n");
      ++nerror;
    }
    else
    {
      // create dimensions
      z_dim  = dataFile->add_dim("z" , grid->kmax);
      zh_dim = dataFile->add_dim("zh", grid->kmax+1);
      t_dim  = dataFile->add_dim("t");

      // create variables
      iter_var = dataFile->add_var("iter", ncInt   , t_dim );
      t_var    = dataFile->add_var("t"   , ncDouble, t_dim );
      z_var    = dataFile->add_var("z"   , ncDouble, z_dim );
      zh_var   = dataFile->add_var("zh"  , ncDouble, zh_dim);

      // means
      addprof("u", "z" );
      addprof("v", "z" );
      addprof("w", "zh");
      for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
        addprof(it->first, "z");
      addprof("evisc", "z");
      addprof("p", "z");

      // 2nd order
      addprof("u2", "z" );
      addprof("v2", "z" );
      addprof("w2", "zh");
      for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
        addprof(it->first+"2", "z");

      // 3rd order
      addprof("u3", "z" );
      addprof("v3", "z" );
      addprof("w3", "zh");
      for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
        addprof(it->first+"3", "z");

      // 4th order
      addprof("u4", "z" );
      addprof("v4", "z" );
      addprof("w4", "zh");
      for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
        addprof(it->first+"4", "z");

      // gradients
      addprof("ugrad", "zh");
      addprof("vgrad", "zh");
      for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
        addprof(it->first+"grad", "zh");

      // turbulent fluxes
      addprof("uw", "zh");
      addprof("vw", "zh");
      for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
        addprof(it->first+"w", "zh");

      // diffusive fluxes
      addprof("udiff", "zh");
      addprof("vdiff", "zh");
      for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
        addprof(it->first+"diff", "zh");

      // total fluxes
      addprof("uflux", "zh");
      addprof("vflux", "zh");
      for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
        addprof(it->first+"flux", "zh");

      // save the grid variables
      z_var ->put(&grid->z [grid->kstart], grid->kmax  );
      zh_var->put(&grid->zh[grid->kstart], grid->kmax+1);

      dataFile->sync();
    }

    initialized = true;
  }

  // crash on all processes in case the file could not be written
  mpi->broadcast(&nerror, 1);

  return (nerror > 0);
}

int cstats_les::exec(int iteration, double time)
{
  if(mpi->mpiid == 0) std::printf("Saving stats for time %f\n", time);

  if(mpi->mpiid == 0)
  {
    t_var   ->put_rec(&time     , nstats);
    iter_var->put_rec(&iteration, nstats);
  }

  // PROFILES
  // calculate means
  calcmean(fields->u->data, profs["u"].data, NO_OFFSET);
  calcmean(fields->v->data, profs["v"].data, NO_OFFSET);
  calcmean(fields->w->data, profs["w"].data, NO_OFFSET);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    calcmean(it->second->data, profs[it->first].data, NO_OFFSET);

  calcmean(fields->s["p"]->data, profs["p"].data, NO_OFFSET);
  calcmean(fields->s["evisc"]->data, profs["evisc"].data, NO_OFFSET);

  // calculate absolute means
  calcmean(fields->u->data, uabs, grid->u);
  calcmean(fields->v->data, vabs, grid->v);

  // 2nd order
  calcmoment(fields->u->data, profs["u"].data, profs["u2"].data, 2., 0);
  calcmoment(fields->v->data, profs["v"].data, profs["v2"].data, 2., 0);
  calcmoment(fields->w->data, profs["w"].data, profs["w2"].data, 2., 1);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    calcmoment(it->second->data, profs[it->first].data, profs[it->first+"2"].data, 2., 0);

  // 3rd order
  calcmoment(fields->u->data, profs["u"].data, profs["u3"].data, 3., 0);
  calcmoment(fields->v->data, profs["v"].data, profs["v3"].data, 3., 0);
  calcmoment(fields->w->data, profs["w"].data, profs["w3"].data, 3., 1);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    calcmoment(it->second->data, profs[it->first].data, profs[it->first+"3"].data, 3., 0);

  // 4th order
  calcmoment(fields->u->data, profs["u"].data, profs["u4"].data, 4., 0);
  calcmoment(fields->v->data, profs["v"].data, profs["v4"].data, 4., 0);
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
  calcdiff(fields->u->data, fields->s["evisc"]->data, profs["udiff"].data, grid->dzhi, fields->u->datafluxbot, fields->u->datafluxtop, 1.);
  calcdiff(fields->v->data, fields->s["evisc"]->data, profs["vdiff"].data, grid->dzhi, fields->v->datafluxbot, fields->v->datafluxtop, 1.);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    calcdiff(it->second->data, fields->s["evisc"]->data, profs[it->first+"diff"].data, grid->dzhi, it->second->datafluxbot, it->second->datafluxtop, fields->tPr);

  addfluxes(profs["uflux"].data, profs["uw"].data, profs["udiff"].data);
  addfluxes(profs["vflux"].data, profs["vw"].data, profs["vdiff"].data);
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    addfluxes(profs[it->first+"flux"].data, profs[it->first+"w"].data, profs[it->first+"diff"].data);

  // put the data into the NetCDF file
  if(mpi->mpiid == 0)
  {
    for(profmap::const_iterator it=profs.begin(); it!=profs.end(); ++it)
      profs[it->first].ncvar->put_rec(&profs[it->first].data[grid->kstart], nstats);
  }

  // sync the data
  if(mpi->mpiid == 0) 
    dataFile->sync();

  ++nstats;

  return 0;
}

int cstats_les::addprof(std::string name, std::string zloc)
{
  // create the NetCDF variable
  if(zloc == "z")
    profs[name].ncvar = dataFile->add_var(name.c_str(), ncDouble, t_dim, z_dim );
  else if(zloc == "zh")
    profs[name].ncvar = dataFile->add_var(name.c_str(), ncDouble, t_dim, zh_dim);

  // and allocate the memory
  profs[name].data = new double[grid->kcells];

  return 0;
}

// COMPUTATIONAL KERNELS BELOW
int cstats_les::calcmean(double * restrict data, double * restrict prof, double offset)
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

int cstats_les::calcmoment(double * restrict data, double * restrict datamean, double * restrict prof, double power, int a)
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

int cstats_les::calcflux(double * restrict data, double * restrict w, double * restrict prof, double * restrict tmp1, int locx, int locy)
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

int cstats_les::calcgrad(double * restrict data, double * restrict prof, double * restrict dzhi)
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

int cstats_les::calcdiff(double * restrict data, double * restrict evisc, double * restrict prof, double * restrict dzhi, double * restrict fluxbot, double * restrict fluxtop, double tPr)
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

int cstats_les::addfluxes(double * restrict flux, double * restrict turb, double * restrict diff)
{
  for(int k=grid->kstart; k<grid->kend+1; ++k)
    flux[k] = turb[k] + diff[k];

  return 0;
}

