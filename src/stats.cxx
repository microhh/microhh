#include <cstdio>
#include <netcdfcpp.h>
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "defines.h"

cstats::cstats(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cstats::~cstats()
{
  if(mpi->mpiid == 0)
    delete dataFile;
}

int cstats::readinifile(cinput *inputin)
{
  int n = 0;

  n += inputin->getItem(&istats, "postproc", "stats");

  if(n > 0)
    return 1;

  return 0;
}

int cstats::init()
{
  // create a NetCDF file for the statistics
  if(mpi->mpiid == 0)
  {
    dataFile = new NcFile("test.nc", NcFile::Replace);
    if(!dataFile->is_valid())
    {
      std::printf("ERROR: can't write statistics file");
      return 1;
    }

    // create dimensions
    z_dim  = dataFile->add_dim("z" , grid->kmax);
    zh_dim = dataFile->add_dim("zh", grid->kmax+1);
    t_dim  = dataFile->add_dim("t");

    // create variables
    t_var  = dataFile->add_var("t" , ncDouble, t_dim );
    z_var  = dataFile->add_var("z" , ncDouble, z_dim );
    zh_var = dataFile->add_var("zh", ncDouble, zh_dim);

    u_var = dataFile->add_var("u", ncDouble, t_dim, z_dim );
    v_var = dataFile->add_var("v", ncDouble, t_dim, z_dim );
    w_var = dataFile->add_var("w", ncDouble, t_dim, zh_dim);
    s_var = dataFile->add_var("s", ncDouble, t_dim, z_dim );

    u2_var = dataFile->add_var("u2", ncDouble, t_dim, z_dim );
    v2_var = dataFile->add_var("v2", ncDouble, t_dim, z_dim );
    w2_var = dataFile->add_var("w2", ncDouble, t_dim, zh_dim);
    s2_var = dataFile->add_var("s2", ncDouble, t_dim, z_dim );

    wu_var = dataFile->add_var("wu", ncDouble, t_dim, zh_dim );
    wv_var = dataFile->add_var("wv", ncDouble, t_dim, zh_dim );
    ws_var = dataFile->add_var("ws", ncDouble, t_dim, zh_dim );

    wud_var = dataFile->add_var("wud", ncDouble, t_dim, zh_dim );
    wvd_var = dataFile->add_var("wvd", ncDouble, t_dim, zh_dim );
    wsd_var = dataFile->add_var("wsd", ncDouble, t_dim, zh_dim );

    uflux_var = dataFile->add_var("uflux", ncDouble, t_dim, zh_dim );
    vflux_var = dataFile->add_var("vflux", ncDouble, t_dim, zh_dim );
    sflux_var = dataFile->add_var("sflux", ncDouble, t_dim, zh_dim );

    // save the grid variables
    z_var ->put(&grid->z [grid->kstart], grid->kmax  );
    zh_var->put(&grid->zh[grid->kstart], grid->kmax+1);

    dataFile->sync();
  }

  u = new double[grid->kcells];
  v = new double[grid->kcells];
  w = new double[grid->kcells];
  s = new double[grid->kcells];

  u2 = new double[grid->kcells];
  v2 = new double[grid->kcells];
  w2 = new double[grid->kcells];
  s2 = new double[grid->kcells];

  wu = new double[grid->kcells];
  wv = new double[grid->kcells];
  ws = new double[grid->kcells];

  wud = new double[grid->kcells];
  wvd = new double[grid->kcells];
  wsd = new double[grid->kcells];

  uflux = new double[grid->kcells];
  vflux = new double[grid->kcells];
  sflux = new double[grid->kcells];

  // set the number of stats to zero
  nstats = 0;

  return 0;
}

int cstats::exec(int iteration, double time)
{
  if(mpi->mpiid == 0) std::printf("Saving stats for iteration %d\n", iteration);

  if(mpi->mpiid == 0)
    t_var->put_rec(&time, nstats);

  // PROFILES
  // calculate means
  calcmean((*fields->u).data, u, 0);
  calcmean((*fields->v).data, v, 0);
  calcmean((*fields->w).data, w, 1);
  calcmean((*fields->s).data, s, 0);

  calcvar((*fields->u).data, u, u2, 0);
  calcvar((*fields->v).data, v, v2, 0);
  calcvar((*fields->w).data, w, w2, 1);
  calcvar((*fields->s).data, s, s2, 0);

  calcflux((*fields->u).data, (*fields->w).data, wu);
  calcflux((*fields->v).data, (*fields->w).data, wv);
  calcflux((*fields->s).data, (*fields->w).data, ws);

  calcdiff((*fields->u).data, wud, grid->dzhi4, fields->visc );
  calcdiff((*fields->v).data, wvd, grid->dzhi4, fields->visc );
  calcdiff((*fields->s).data, wsd, grid->dzhi4, fields->viscs);

  for(int k=grid->kstart; k<grid->kend+1; k++)
  {
    uflux[k] = wu[k] + wud[k];
    vflux[k] = wv[k] + wvd[k];
    sflux[k] = ws[k] + wsd[k];
  }

  if(mpi->mpiid == 0)
  {
    u_var->put_rec(&u[grid->kstart], nstats);
    v_var->put_rec(&v[grid->kstart], nstats);
    w_var->put_rec(&w[grid->kstart], nstats);
    s_var->put_rec(&s[grid->kstart], nstats);

    u2_var->put_rec(&u2[grid->kstart], nstats);
    v2_var->put_rec(&v2[grid->kstart], nstats);
    w2_var->put_rec(&w2[grid->kstart], nstats);
    s2_var->put_rec(&s2[grid->kstart], nstats);

    wu_var->put_rec(&wu[grid->kstart], nstats);
    wv_var->put_rec(&wv[grid->kstart], nstats);
    ws_var->put_rec(&ws[grid->kstart], nstats);

    wud_var->put_rec(&wud[grid->kstart], nstats);
    wvd_var->put_rec(&wvd[grid->kstart], nstats);
    wsd_var->put_rec(&wsd[grid->kstart], nstats);

    uflux_var->put_rec(&uflux[grid->kstart], nstats);
    vflux_var->put_rec(&vflux[grid->kstart], nstats);
    sflux_var->put_rec(&sflux[grid->kstart], nstats);
  }

  // calculate variance and higher order moments
  // save variances

  // calculate flux
  // save flux
  
  // SCALARS

  // sync the data
  if(mpi->mpiid == 0) 
    dataFile->sync();

  nstats++;

  return 0;
}

int cstats::calcmean(double * restrict data, double * restrict prof, int a)
{
  int    ijk,ii,jj,kk,kstart;
  double dxi,dyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  
  kstart = grid->kstart;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend+a; k++)
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

  for(int k=grid->kstart; k<grid->kend+a; k++)
    prof[k] /= n;

  grid->getprof(prof);

  return 0;
}

int cstats::calcvar(double * restrict data, double * restrict datamean, double * restrict prof, int a)
{
  int    ijk,ii,jj,kk,kstart;
  double dxi,dyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  
  kstart = grid->kstart;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend+a; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += (data[ijk]-datamean[k])*(data[ijk]-datamean[k]);
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+a; k++)
    prof[k] /= n;

  grid->getprof(prof);

  return 0;
}

int cstats::calcflux(double * restrict data, double * restrict w, double * restrict prof)
{
  int    ijk,ii,jj,kk1,kk2,kstart;
  double dxi,dyi;

  ii  = 1;
  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  
  kstart = grid->kstart;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend+1; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk1;
        prof[k] += (ci0*data[ijk-kk2] + ci1*data[ijk-kk1] + ci2*data[ijk] + ci3*data[ijk+kk1])*w[ijk];
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+1; k++)
    prof[k] /= n;

  grid->getprof(prof);

  return 0;
}

int cstats::calcdiff(double * restrict data, double * restrict prof, double * restrict dzhi4, double visc)
{
  int    ijk,ii,jj,kk1,kk2,kstart;
  double dxi,dyi;

  ii  = 1;
  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  
  kstart = grid->kstart;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend+1; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk1;
        prof[k] += -visc*(cg0*data[ijk-kk2] + cg1*data[ijk-kk1] + cg2*data[ijk] + cg3*data[ijk+kk1])*dzhi4[k];
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+1; k++)
    prof[k] /= n;

  grid->getprof(prof);

  return 0;
}
