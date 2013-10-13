#include <cstdio>
#include <cmath>
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "defines.h"
#include <netcdfcpp.h>

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
    delete[] u.data;
    delete[] v.data;
    delete[] w.data;
    delete[] s.data;

    delete[] u2.data;
    delete[] v2.data;
    delete[] w2.data;
    delete[] s2.data;

    delete[] u3.data;
    delete[] v3.data;
    delete[] w3.data;
    delete[] s3.data;

    delete[] evisc.data;

    delete[] ugrad.data;
    delete[] vgrad.data;
    delete[] sgrad.data;

    delete[] wu.data;
    delete[] wv.data;
    delete[] ws.data;

    delete[] udiff.data;
    delete[] vdiff.data;
    delete[] sdiff.data;

    delete[] uflux.data;
    delete[] vflux.data;
    delete[] sflux.data;
  }
}

int cstats_les::init()
{
  u.data = new double[grid->kcells];
  v.data = new double[grid->kcells];
  w.data = new double[grid->kcells];
  s.data = new double[grid->kcells];

  u2.data = new double[grid->kcells];
  v2.data = new double[grid->kcells];
  w2.data = new double[grid->kcells];
  s2.data = new double[grid->kcells];

  u3.data = new double[grid->kcells];
  v3.data = new double[grid->kcells];
  w3.data = new double[grid->kcells];
  s3.data = new double[grid->kcells];

  evisc.data = new double[grid->kcells];

  ugrad.data = new double[grid->kcells];
  vgrad.data = new double[grid->kcells];
  sgrad.data = new double[grid->kcells];

  wu.data = new double[grid->kcells];
  wv.data = new double[grid->kcells];
  ws.data = new double[grid->kcells];

  udiff.data = new double[grid->kcells];
  vdiff.data = new double[grid->kcells];
  sdiff.data = new double[grid->kcells];

  uflux.data = new double[grid->kcells];
  vflux.data = new double[grid->kcells];
  sflux.data = new double[grid->kcells];

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

      u.ncvar = dataFile->add_var("u", ncDouble, t_dim, z_dim );
      v.ncvar = dataFile->add_var("v", ncDouble, t_dim, z_dim );
      w.ncvar = dataFile->add_var("w", ncDouble, t_dim, zh_dim);
      s.ncvar = dataFile->add_var("s", ncDouble, t_dim, z_dim );

      u2.ncvar = dataFile->add_var("u2", ncDouble, t_dim, z_dim );
      v2.ncvar = dataFile->add_var("v2", ncDouble, t_dim, z_dim );
      w2.ncvar = dataFile->add_var("w2", ncDouble, t_dim, zh_dim);
      s2.ncvar = dataFile->add_var("s2", ncDouble, t_dim, z_dim );

      u3.ncvar = dataFile->add_var("u3", ncDouble, t_dim, z_dim );
      v3.ncvar = dataFile->add_var("v3", ncDouble, t_dim, z_dim );
      w3.ncvar = dataFile->add_var("w3", ncDouble, t_dim, zh_dim);
      s3.ncvar = dataFile->add_var("s3", ncDouble, t_dim, z_dim );

      ugrad.ncvar = dataFile->add_var("ugrad", ncDouble, t_dim, zh_dim );
      vgrad.ncvar = dataFile->add_var("vgrad", ncDouble, t_dim, zh_dim );
      sgrad.ncvar = dataFile->add_var("sgrad", ncDouble, t_dim, zh_dim );

      wu.ncvar = dataFile->add_var("uw", ncDouble, t_dim, zh_dim );
      wv.ncvar = dataFile->add_var("vw", ncDouble, t_dim, zh_dim );
      ws.ncvar = dataFile->add_var("sw", ncDouble, t_dim, zh_dim );

      udiff.ncvar = dataFile->add_var("udiff", ncDouble, t_dim, zh_dim );
      vdiff.ncvar = dataFile->add_var("vdiff", ncDouble, t_dim, zh_dim );
      sdiff.ncvar = dataFile->add_var("sdiff", ncDouble, t_dim, zh_dim );

      uflux.ncvar = dataFile->add_var("uflux", ncDouble, t_dim, zh_dim );
      vflux.ncvar = dataFile->add_var("vflux", ncDouble, t_dim, zh_dim );
      sflux.ncvar = dataFile->add_var("sflux", ncDouble, t_dim, zh_dim );

      evisc.ncvar = dataFile->add_var("evisc", ncDouble, t_dim, z_dim );

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
  calcmean(fields->u->data, u.data);
  calcmean(fields->v->data, v.data);
  calcmean(fields->w->data, w.data);
  calcmean(fields->s["s"]->data, s.data);
  calcmean(fields->s["evisc"]->data, evisc.data);

  // calc variances
  calcmoment(fields->u->data, u.data, u2.data, 2., 0);
  calcmoment(fields->v->data, v.data, v2.data, 2., 0);
  calcmoment(fields->w->data, w.data, w2.data, 2., 1);
  calcmoment(fields->s["s"]->data, s.data, s2.data, 2., 0);

  // calc skewnesses
  calcmoment(fields->u->data, u.data, u3.data, 3., 0);
  calcmoment(fields->v->data, v.data, v3.data, 3., 0);
  calcmoment(fields->w->data, w.data, w3.data, 3., 1);
  calcmoment(fields->s["s"]->data, s.data, s3.data, 3., 0);

  calcgrad(fields->u->data, ugrad.data, grid->dzhi);
  calcgrad(fields->v->data, vgrad.data, grid->dzhi);
  calcgrad(fields->s["s"]->data, sgrad.data, grid->dzhi);

  // calculate turbulent fluxes
  calcflux(fields->u->data, fields->w->data, wu.data, fields->s["tmp1"]->data, 1, 0);
  calcflux(fields->v->data, fields->w->data, wv.data, fields->s["tmp1"]->data, 0, 1);
  calcflux(fields->s["s"]->data, fields->w->data, ws.data, fields->s["tmp1"]->data, 0, 0);

  // calculate diffusive fluxes
  calcdiff(fields->u->data, fields->s["evisc"]->data, udiff.data, grid->dzhi, fields->u->datafluxbot, fields->u->datafluxtop, 1.);
  calcdiff(fields->v->data, fields->s["evisc"]->data, vdiff.data, grid->dzhi, fields->v->datafluxbot, fields->v->datafluxtop, 1.);
  calcdiff(fields->s["s"]->data, fields->s["evisc"]->data, sdiff.data, grid->dzhi, fields->s["s"]->datafluxbot, fields->s["s"]->datafluxtop, fields->tPr);

  // add the turbulent and diffusive fluxes
  for(int k=grid->kstart; k<grid->kend+1; k++)
  {
    uflux.data[k] = wu.data[k] + udiff.data[k];
    vflux.data[k] = wv.data[k] + vdiff.data[k];
    sflux.data[k] = ws.data[k] + sdiff.data[k];
  }

  // put the data into the NetCDF file
  if(mpi->mpiid == 0)
  {
    u.ncvar->put_rec(&u.data[grid->kstart], nstats);
    v.ncvar->put_rec(&v.data[grid->kstart], nstats);
    w.ncvar->put_rec(&w.data[grid->kstart], nstats);
    s.ncvar->put_rec(&s.data[grid->kstart], nstats);

    evisc.ncvar->put_rec(&evisc.data[grid->kstart], nstats);

    u2.ncvar->put_rec(&u2.data[grid->kstart], nstats);
    v2.ncvar->put_rec(&v2.data[grid->kstart], nstats);
    w2.ncvar->put_rec(&w2.data[grid->kstart], nstats);
    s2.ncvar->put_rec(&s2.data[grid->kstart], nstats);

    u3.ncvar->put_rec(&u3.data[grid->kstart], nstats);
    v3.ncvar->put_rec(&v3.data[grid->kstart], nstats);
    w3.ncvar->put_rec(&w3.data[grid->kstart], nstats);
    s3.ncvar->put_rec(&s3.data[grid->kstart], nstats);

    ugrad.ncvar->put_rec(&ugrad.data[grid->kstart], nstats);
    vgrad.ncvar->put_rec(&vgrad.data[grid->kstart], nstats);
    sgrad.ncvar->put_rec(&sgrad.data[grid->kstart], nstats);

    wu.ncvar->put_rec(&wu.data[grid->kstart], nstats);
    wv.ncvar->put_rec(&wv.data[grid->kstart], nstats);
    ws.ncvar->put_rec(&ws.data[grid->kstart], nstats);

    udiff.ncvar->put_rec(&udiff.data[grid->kstart], nstats);
    vdiff.ncvar->put_rec(&vdiff.data[grid->kstart], nstats);
    sdiff.ncvar->put_rec(&sdiff.data[grid->kstart], nstats);

    uflux.ncvar->put_rec(&uflux.data[grid->kstart], nstats);
    vflux.ncvar->put_rec(&vflux.data[grid->kstart], nstats);
    sflux.ncvar->put_rec(&sflux.data[grid->kstart], nstats);
  }

  // sync the data
  if(mpi->mpiid == 0) 
    dataFile->sync();

  nstats++;

  return 0;
}

// BELOW ARE ALL COMPUTATIONAL KERNELS
int cstats_les::calcmean(double * restrict data, double * restrict prof)
{
  int ijk,ii,jj,kk;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  
  for(int k=0; k<grid->kcells; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += data[ijk];
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=0; k<grid->kcells; k++)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats_les::calcmoment(double * restrict data, double * restrict datamean, double * restrict prof, double power, int a)
{
  int ijk,ii,jj,kk;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  
  for(int k=grid->kstart; k<grid->kend+a; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += std::pow(data[ijk]-datamean[k], power);
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+a; k++)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats_les::calcflux(double * restrict data, double * restrict w, double * restrict prof, double * restrict tmp1, int locx, int locy)
{
  int ijk,ii,jj,kk;

  ii = 1;
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
  
  for(int k=grid->kstart; k<grid->kend+1; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += 0.5*(data[ijk-kk]+data[ijk])*calcw[ijk];
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+1; k++)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats_les::calcgrad(double * restrict data, double * restrict prof, double * restrict dzhi)
{
  int ijk,ii,jj,kk;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  
  for(int k=grid->kstart; k<grid->kend+1; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += (data[ijk]-data[ijk-kk])*dzhi[k];
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+1; k++)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats_les::calcdiff(double * restrict data, double * restrict evisc, double * restrict prof, double * restrict dzhi, double * restrict fluxbot, double * restrict fluxtop, double tPr)
{
  int ijk,ij,ii,jj,kk,kstart,kend;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  kstart = grid->kstart;
  kend   = grid->kend;

  // CvH add horizontal interpolation for u and v and interpolate the eddy viscosity properly
  // bottom boundary
  prof[kstart] = 0.;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ij = i + j*jj;
      prof[kstart] += fluxbot[ij];
    }

  for(int k=grid->kstart+1; k<grid->kend; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += -0.5*(evisc[ijk-kk]+evisc[ijk])/tPr*(data[ijk]-data[ijk-kk])*dzhi[k];
      }
  }

  // top boundary
  prof[kend] = 0.;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ij = i + j*jj;
      prof[kend] += fluxtop[ij];
    }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+1; k++)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

