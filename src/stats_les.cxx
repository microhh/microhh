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
    delete[] u;
    delete[] v;
    delete[] w;
    delete[] s;

    delete[] u2;
    delete[] v2;
    delete[] w2;
    delete[] s2;

    delete[] u3;
    delete[] v3;
    delete[] w3;
    delete[] s3;

    delete[] evisc;

    delete[] ugrad;
    delete[] vgrad;
    delete[] sgrad;

    delete[] wu;
    delete[] wv;
    delete[] ws;

    delete[] udiff;
    delete[] vdiff;
    delete[] sdiff;

    delete[] uflux;
    delete[] vflux;
    delete[] sflux;
  }
}

int cstats_les::init()
{
  u = new double[grid->kcells];
  v = new double[grid->kcells];
  w = new double[grid->kcells];
  s = new double[grid->kcells];

  u2 = new double[grid->kcells];
  v2 = new double[grid->kcells];
  w2 = new double[grid->kcells];
  s2 = new double[grid->kcells];

  u3 = new double[grid->kcells];
  v3 = new double[grid->kcells];
  w3 = new double[grid->kcells];
  s3 = new double[grid->kcells];

  evisc = new double[grid->kcells];

  ugrad = new double[grid->kcells];
  vgrad = new double[grid->kcells];
  sgrad = new double[grid->kcells];

  wu = new double[grid->kcells];
  wv = new double[grid->kcells];
  ws = new double[grid->kcells];

  udiff = new double[grid->kcells];
  vdiff = new double[grid->kcells];
  sdiff = new double[grid->kcells];

  uflux = new double[grid->kcells];
  vflux = new double[grid->kcells];
  sflux = new double[grid->kcells];

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

      u_var = dataFile->add_var("u", ncDouble, t_dim, z_dim );
      v_var = dataFile->add_var("v", ncDouble, t_dim, z_dim );
      w_var = dataFile->add_var("w", ncDouble, t_dim, zh_dim);
      s_var = dataFile->add_var("s", ncDouble, t_dim, z_dim );

      u2_var = dataFile->add_var("u2", ncDouble, t_dim, z_dim );
      v2_var = dataFile->add_var("v2", ncDouble, t_dim, z_dim );
      w2_var = dataFile->add_var("w2", ncDouble, t_dim, zh_dim);
      s2_var = dataFile->add_var("s2", ncDouble, t_dim, z_dim );

      u3_var = dataFile->add_var("u3", ncDouble, t_dim, z_dim );
      v3_var = dataFile->add_var("v3", ncDouble, t_dim, z_dim );
      w3_var = dataFile->add_var("w3", ncDouble, t_dim, zh_dim);
      s3_var = dataFile->add_var("s3", ncDouble, t_dim, z_dim );

      ugrad_var = dataFile->add_var("ugrad", ncDouble, t_dim, zh_dim );
      vgrad_var = dataFile->add_var("vgrad", ncDouble, t_dim, zh_dim );
      sgrad_var = dataFile->add_var("sgrad", ncDouble, t_dim, zh_dim );

      wu_var = dataFile->add_var("uw", ncDouble, t_dim, zh_dim );
      wv_var = dataFile->add_var("vw", ncDouble, t_dim, zh_dim );
      ws_var = dataFile->add_var("sw", ncDouble, t_dim, zh_dim );

      udiff_var = dataFile->add_var("udiff", ncDouble, t_dim, zh_dim );
      vdiff_var = dataFile->add_var("vdiff", ncDouble, t_dim, zh_dim );
      sdiff_var = dataFile->add_var("sdiff", ncDouble, t_dim, zh_dim );

      uflux_var = dataFile->add_var("uflux", ncDouble, t_dim, zh_dim );
      vflux_var = dataFile->add_var("vflux", ncDouble, t_dim, zh_dim );
      sflux_var = dataFile->add_var("sflux", ncDouble, t_dim, zh_dim );

      evisc_var = dataFile->add_var("evisc", ncDouble, t_dim, z_dim );

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
  calcmean(fields->u->data, u);
  calcmean(fields->v->data, v);
  calcmean(fields->w->data, w);
  calcmean(fields->s["s"]->data, s);
  calcmean(fields->s["evisc"]->data, evisc);

  // calc variances
  calcmoment(fields->u->data, u, u2, 2., 0);
  calcmoment(fields->v->data, v, v2, 2., 0);
  calcmoment(fields->w->data, w, w2, 2., 1);
  calcmoment(fields->s["s"]->data, s, s2, 2., 0);

  // calc skewnesses
  calcmoment(fields->u->data, u, u3, 3., 0);
  calcmoment(fields->v->data, v, v3, 3., 0);
  calcmoment(fields->w->data, w, w3, 3., 1);
  calcmoment(fields->s["s"]->data, s, s3, 3., 0);

  calcgrad(fields->u->data, ugrad, grid->dzhi);
  calcgrad(fields->v->data, vgrad, grid->dzhi);
  calcgrad(fields->s["s"]->data, sgrad, grid->dzhi);

  // calculate turbulent fluxes
  calcflux(fields->u->data, fields->w->data, wu, fields->s["tmp1"]->data, 1, 0);
  calcflux(fields->v->data, fields->w->data, wv, fields->s["tmp1"]->data, 0, 1);
  calcflux(fields->s["s"]->data, fields->w->data, ws, fields->s["tmp1"]->data, 0, 0);

  // calculate diffusive fluxes
  calcdiff(fields->u->data, fields->s["evisc"]->data, udiff, grid->dzhi, fields->u->datafluxbot, fields->u->datafluxtop, 1.);
  calcdiff(fields->v->data, fields->s["evisc"]->data, vdiff, grid->dzhi, fields->v->datafluxbot, fields->v->datafluxtop, 1.);
  calcdiff(fields->s["s"]->data, fields->s["evisc"]->data, sdiff, grid->dzhi, fields->s["s"]->datafluxbot, fields->s["s"]->datafluxtop, fields->tPr);

  // add the turbulent and diffusive fluxes
  for(int k=grid->kstart; k<grid->kend+1; k++)
  {
    uflux[k] = wu[k] + udiff[k];
    vflux[k] = wv[k] + vdiff[k];
    sflux[k] = ws[k] + sdiff[k];
  }

  // put the data into the NetCDF file
  if(mpi->mpiid == 0)
  {
    u_var->put_rec(&u[grid->kstart], nstats);
    v_var->put_rec(&v[grid->kstart], nstats);
    w_var->put_rec(&w[grid->kstart], nstats);
    s_var->put_rec(&s[grid->kstart], nstats);

    evisc_var->put_rec(&evisc[grid->kstart], nstats);

    u2_var->put_rec(&u2[grid->kstart], nstats);
    v2_var->put_rec(&v2[grid->kstart], nstats);
    w2_var->put_rec(&w2[grid->kstart], nstats);
    s2_var->put_rec(&s2[grid->kstart], nstats);

    u3_var->put_rec(&u3[grid->kstart], nstats);
    v3_var->put_rec(&v3[grid->kstart], nstats);
    w3_var->put_rec(&w3[grid->kstart], nstats);
    s3_var->put_rec(&s3[grid->kstart], nstats);

    ugrad_var->put_rec(&ugrad[grid->kstart], nstats);
    vgrad_var->put_rec(&vgrad[grid->kstart], nstats);
    sgrad_var->put_rec(&sgrad[grid->kstart], nstats);

    wu_var->put_rec(&wu[grid->kstart], nstats);
    wv_var->put_rec(&wv[grid->kstart], nstats);
    ws_var->put_rec(&ws[grid->kstart], nstats);

    udiff_var->put_rec(&udiff[grid->kstart], nstats);
    vdiff_var->put_rec(&vdiff[grid->kstart], nstats);
    sdiff_var->put_rec(&sdiff[grid->kstart], nstats);

    uflux_var->put_rec(&uflux[grid->kstart], nstats);
    vflux_var->put_rec(&vflux[grid->kstart], nstats);
    sflux_var->put_rec(&sflux[grid->kstart], nstats);
  }

  // sync the data
  if(mpi->mpiid == 0) 
    dataFile->sync();

  nstats++;

  return 0;
}

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

