#include <cstdio>
#include <cmath>
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

  // optional, by default switch stats on
  n += inputin->getItem(&istats, "postproc", "stats", 1);

  if(n > 0)
    return 1;

  return 0;
}

int cstats::init()
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

  u2_shear  = new double[grid->kcells];
  v2_shear  = new double[grid->kcells];
  tke_shear = new double[grid->kcells];

  u2_turb  = new double[grid->kcells];
  v2_turb  = new double[grid->kcells];
  tke_turb = new double[grid->kcells];

  // set the number of stats to zero
  nstats = 0;

  return 0;
}

int cstats::create(std::string simname, int n)
{
  // create a NetCDF file for the statistics
  if(mpi->mpiid == 0)
  {
    char filename[256];
    std::sprintf(filename, "%s.%07d.nc", simname.c_str(), n);
    dataFile = new NcFile(filename, NcFile::New);
    if(!dataFile->is_valid())
    {
      std::printf("ERROR cannot write statistics file\n");
      return 1;
    }

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

    wu_var = dataFile->add_var("wu", ncDouble, t_dim, zh_dim );
    wv_var = dataFile->add_var("wv", ncDouble, t_dim, zh_dim );
    ws_var = dataFile->add_var("ws", ncDouble, t_dim, zh_dim );

    udiff_var = dataFile->add_var("udiff", ncDouble, t_dim, zh_dim );
    vdiff_var = dataFile->add_var("vdiff", ncDouble, t_dim, zh_dim );
    sdiff_var = dataFile->add_var("sdiff", ncDouble, t_dim, zh_dim );

    uflux_var = dataFile->add_var("uflux", ncDouble, t_dim, zh_dim );
    vflux_var = dataFile->add_var("vflux", ncDouble, t_dim, zh_dim );
    sflux_var = dataFile->add_var("sflux", ncDouble, t_dim, zh_dim );

    u2_shear_var  = dataFile->add_var("u2_shear" , ncDouble, t_dim, z_dim );
    v2_shear_var  = dataFile->add_var("v2_shear" , ncDouble, t_dim, z_dim );
    tke_shear_var = dataFile->add_var("tke_shear", ncDouble, t_dim, z_dim );

    u2_turb_var  = dataFile->add_var("u2_turb" , ncDouble, t_dim, z_dim );
    v2_turb_var  = dataFile->add_var("v2_turb" , ncDouble, t_dim, z_dim );
    tke_turb_var = dataFile->add_var("tke_turb", ncDouble, t_dim, z_dim );

    // save the grid variables
    z_var ->put(&grid->z [grid->kstart], grid->kmax  );
    zh_var->put(&grid->zh[grid->kstart], grid->kmax+1);

    dataFile->sync();
  }

  return 0;
}

int cstats::exec(int iteration, double time)
{
  if(mpi->mpiid == 0) std::printf("Saving stats for iteration %d\n", iteration);

  if(mpi->mpiid == 0)
  {
    t_var   ->put_rec(&time     , nstats);
    iter_var->put_rec(&iteration, nstats);
  }

  // PROFILES
  // calculate means
  calcmean((*fields->u).data, u);
  calcmean((*fields->v).data, v);
  calcmean((*fields->w).data, w);
  calcmean((*fields->s).data, s);

  // calc variances
  calcmoment((*fields->u).data, u, u2, 2., 0);
  calcmoment((*fields->v).data, v, v2, 2., 0);
  calcmoment((*fields->w).data, w, w2, 2., 1);
  calcmoment((*fields->s).data, s, s2, 2., 0);

  // calc skewnesses
  calcmoment((*fields->u).data, u, u3, 3., 0);
  calcmoment((*fields->v).data, v, v3, 3., 0);
  calcmoment((*fields->w).data, w, w3, 3., 1);
  calcmoment((*fields->s).data, s, s3, 3., 0);

  // calculate gradients
  calcgrad((*fields->u).data, ugrad, grid->dzhi4);
  calcgrad((*fields->v).data, vgrad, grid->dzhi4);
  calcgrad((*fields->s).data, sgrad, grid->dzhi4);

  // calculate turbulent fluxes
  calcflux((*fields->u).data, (*fields->w).data, wu, (*fields->tmp1).data, 1, 0);
  calcflux((*fields->v).data, (*fields->w).data, wv, (*fields->tmp1).data, 0, 1);
  calcflux((*fields->s).data, (*fields->w).data, ws, (*fields->tmp1).data, 0, 0);

  // calculate diffusive fluxes
  calcdiff((*fields->u).data, udiff, grid->dzhi4, fields->visc );
  calcdiff((*fields->v).data, vdiff, grid->dzhi4, fields->visc );
  calcdiff((*fields->s).data, sdiff, grid->dzhi4, fields->viscs);

  // add the turbulent and diffusive fluxes
  for(int k=grid->kstart; k<grid->kend+1; k++)
  {
    uflux[k] = wu[k] + udiff[k];
    vflux[k] = wv[k] + vdiff[k];
    sflux[k] = ws[k] + sdiff[k];
  }

  // calculate the TKE budget
  calctkebudget((*fields->u).data, (*fields->v).data, (*fields->w).data, (*fields->tmp1).data, (*fields->tmp2).data,
                u, v,
                u2_shear, v2_shear, tke_shear,
                u2_turb , v2_turb , tke_turb,
                grid->dzi4);

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

    u2_shear_var ->put_rec(&u2_shear [grid->kstart], nstats);
    v2_shear_var ->put_rec(&v2_shear [grid->kstart], nstats);
    tke_shear_var->put_rec(&tke_shear[grid->kstart], nstats);

    u2_turb_var ->put_rec(&u2_turb [grid->kstart], nstats);
    v2_turb_var ->put_rec(&v2_turb [grid->kstart], nstats);
    tke_turb_var->put_rec(&tke_turb[grid->kstart], nstats);
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

int cstats::calcmean(double * restrict data, double * restrict prof)
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

int cstats::calcmoment(double * restrict data, double * restrict datamean, double * restrict prof, double power, int a)
{
  int ijk,ii,jj,kk;
  double dxi,dyi;

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

int cstats::calcflux(double * restrict data, double * restrict w, double * restrict prof, double * restrict tmp1, int locx, int locy)
{
  int ijk,ii,jj,kk1,kk2;

  ii  = 1;
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
  
  for(int k=grid->kstart; k<grid->kend+1; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk1;
        prof[k] += (ci0*data[ijk-kk2] + ci1*data[ijk-kk1] + ci2*data[ijk] + ci3*data[ijk+kk1])*calcw[ijk];
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+1; k++)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats::calcgrad(double * restrict data, double * restrict prof, double * restrict dzhi4)
{
  int ijk,ii,jj,kk1,kk2;

  ii  = 1;
  jj  = 1*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  
  for(int k=grid->kstart; k<grid->kend+1; k++)
  {
    prof[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk1;
        prof[k] += (cg0*data[ijk-kk2] + cg1*data[ijk-kk1] + cg2*data[ijk] + cg3*data[ijk+kk1])*dzhi4[k];
      }
  }

  double n = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend+1; k++)
    prof[k] /= n;

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats::calcdiff(double * restrict data, double * restrict prof, double * restrict dzhi4, double visc)
{
  int ijk,ii,jj,kk1,kk2;

  ii  = 1;
  jj  = 1*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  
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

  grid->getprof(prof, grid->kcells);

  return 0;
}

int cstats::calctkebudget(double * restrict u, double * restrict v, double * restrict w, double * restrict wx, double * restrict wy,
                          double * restrict umean, double * restrict vmean,
                          double * restrict u2_shear, double * restrict v2_shear, double * restrict tke_shear,
                          double * restrict u2_turb , double * restrict v2_turb , double * restrict tke_turb,
                          double * restrict dzi4)
{
  // get w on the x and y location
  grid->interpolatex_4th(wx, w, 0);
  grid->interpolatey_4th(wy, w, 0);

  int ijk,ii,jj,kk1,kk2,kk3;

  ii  = 1;
  jj  = 1*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  double n = grid->imax*grid->jmax;
 
  // calculate the shear term u'w*dumean/dz
  for(int k=grid->kstart; k<grid->kend; k++)
  {
    u2_shear [k] = 0.;
    v2_shear [k] = 0.;
    tke_shear[k] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk1;
        u2_shear[k] += -2.*(u[ijk]-umean[k])*(ci0*wx[ijk-kk1] + ci1*wx[ijk] + ci2*wx[ijk+kk1] + ci3*wx[ijk+kk2])
                     * ( cg0*(ci0*umean[k-3] + ci1*umean[k-2] + ci2*umean[k-1] + ci3*umean[k  ])
                       + cg1*(ci0*umean[k-2] + ci1*umean[k-1] + ci2*umean[k  ] + ci3*umean[k+1])
                       + cg2*(ci0*umean[k-1] + ci1*umean[k  ] + ci2*umean[k+1] + ci3*umean[k+2])
                       + cg3*(ci0*umean[k  ] + ci1*umean[k+1] + ci2*umean[k+2] + ci3*umean[k+3])) * dzi4[k];

        v2_shear[k] += -2.*(v[ijk]-vmean[k])*(ci0*wy[ijk-kk1] + ci1*wy[ijk] + ci2*wy[ijk+kk1] + ci3*wy[ijk+kk2])
                     * ( cg0*(ci0*vmean[k-3] + ci1*vmean[k-2] + ci2*vmean[k-1] + ci3*vmean[k  ])
                       + cg1*(ci0*vmean[k-2] + ci1*vmean[k-1] + ci2*vmean[k  ] + ci3*vmean[k+1])
                       + cg2*(ci0*vmean[k-1] + ci1*vmean[k  ] + ci2*vmean[k+1] + ci3*vmean[k+2])
                       + cg3*(ci0*vmean[k  ] + ci1*vmean[k+1] + ci2*vmean[k+2] + ci3*vmean[k+3])) * dzi4[k];

        tke_shear[k] = 0.5*(u2_shear[k] + v2_shear[k]);
      }
  }

  for(int k=grid->kstart; k<grid->kend; k++)
  {
    u2_shear [k] /= n;
    v2_shear [k] /= n;
    tke_shear[k] /= n;
  }

  grid->getprof(u2_shear , grid->kcells);
  grid->getprof(v2_shear , grid->kcells);
  grid->getprof(tke_shear, grid->kcells);


  // calculate the turbulent transport term
  for(int k=grid->kstart; k<grid->kend; k++)
  {
    u2_turb [k] = 0.;
    v2_turb [k] = 0.;
    tke_turb[k] = 0.;

    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk1;
        u2_turb[k] += ( cg0*(std::pow(ci0*(u[ijk-kk3]-umean[k-3]) + ci1*(u[ijk-kk2]-umean[k-2]) + ci2*(u[ijk-kk1]-umean[k-1]) + ci3*(u[ijk    ]-umean[k  ]), 2.)*wx[ijk-kk1])
                      + cg1*(std::pow(ci0*(u[ijk-kk2]-umean[k-2]) + ci1*(u[ijk-kk1]-umean[k-1]) + ci2*(u[ijk    ]-umean[k  ]) + ci3*(u[ijk+kk1]-umean[k+1]), 2.)*wx[ijk    ])
                      + cg2*(std::pow(ci0*(u[ijk-kk1]-umean[k-1]) + ci1*(u[ijk    ]-umean[k  ]) + ci2*(u[ijk+kk1]-umean[k+1]) + ci3*(u[ijk+kk2]-umean[k+2]), 2.)*wx[ijk+kk1])
                      + cg3*(std::pow(ci0*(u[ijk    ]-umean[k  ]) + ci1*(u[ijk+kk1]-umean[k+1]) + ci2*(u[ijk+kk2]-umean[k+2]) + ci3*(u[ijk+kk3]-umean[k+3]), 2.)*wx[ijk+kk2])) * dzi4[k];

        v2_turb[k] += ( cg0*(std::pow(ci0*(v[ijk-kk3]-vmean[k-3]) + ci1*(v[ijk-kk2]-vmean[k-2]) + ci2*(v[ijk-kk1]-vmean[k-1]) + ci3*(v[ijk    ]-vmean[k  ]), 2.)*wy[ijk-kk1])
                      + cg1*(std::pow(ci0*(v[ijk-kk2]-vmean[k-2]) + ci1*(v[ijk-kk1]-vmean[k-1]) + ci2*(v[ijk    ]-vmean[k  ]) + ci3*(v[ijk+kk1]-vmean[k+1]), 2.)*wy[ijk    ])
                      + cg2*(std::pow(ci0*(v[ijk-kk1]-vmean[k-1]) + ci1*(v[ijk    ]-vmean[k  ]) + ci2*(v[ijk+kk1]-vmean[k+1]) + ci3*(v[ijk+kk2]-vmean[k+2]), 2.)*wy[ijk+kk1])
                      + cg3*(std::pow(ci0*(v[ijk    ]-vmean[k  ]) + ci1*(v[ijk+kk1]-vmean[k+1]) + ci2*(v[ijk+kk2]-vmean[k+2]) + ci3*(v[ijk+kk3]-vmean[k+3]), 2.)*wy[ijk+kk2])) * dzi4[k];

        tke_turb[k] = 0.5*(u2_turb[k] + v2_turb[k]);
      }
  }

  for(int k=grid->kstart; k<grid->kend; k++)
  {
    u2_turb [k] /= n;
    v2_turb [k] /= n;
    tke_turb[k] /= n;
  }

  grid->getprof(u2_turb , grid->kcells);
  grid->getprof(v2_turb , grid->kcells);
  grid->getprof(tke_turb, grid->kcells);

  return 0;
}

