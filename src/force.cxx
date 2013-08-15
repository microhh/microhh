#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "force.h"
#include "defines.h"

cforce::cforce(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;

  allocated = false;
}

cforce::~cforce()
{
  if(allocated)
  {
    if(swforce == "2")
    {
      delete[] ug;
      delete[] vg;
    }
  }
}

int cforce::readinifile(cinput *inputin)
{
  int n = 0;

  n += inputin->getItem(&swforce, "force", "swforce", "");
  
  if(swforce == "1")
    n += inputin->getItem(&uflow, "force", "uflow", "");

  if(swforce == "2")
    n += inputin->getItem(&fc, "force", "fc", "");

  if(n > 0)
    return 1;

  return 0;
}

int cforce::init()
{
  if(swforce == "2")
  {
    ug = new double[grid->kcells];
    vg = new double[grid->kcells];
  }

  allocated = true;

  return 0;
}

int cforce::create(cinput *inputin)
{
  int n = 0;

  if(swforce == "2")
  {
    n += inputin->getProf(&ug[grid->kstart], "ug", grid->kmax);
    n += inputin->getProf(&vg[grid->kstart], "vg", grid->kmax);
  }

  if(n > 0)
    return 1;

  return 0;
}

int cforce::exec(double dt)
{
  if(swforce == "0")
    return 0;

  else if(swforce == "1")
    flux((*fields->ut).data, (*fields->u).data, grid->dz, dt);

  else if(swforce == "2")
  {
    if(grid->swspatialorder == "2")
      coriolis_2nd(fields->ut->data, fields->vt->data, fields->u->data, fields->v->data, ug, vg);
    else if(grid->swspatialorder == "4")
      coriolis_4th(fields->ut->data, fields->vt->data, fields->u->data, fields->v->data, ug, vg);
  }

  return 0;
}

int cforce::save()
{
  // TODO add subsidence to the same save file
  if(swforce == "2")
  {
    char filename[256];
    std::sprintf(filename, "%s.%07d", "force", 0);

    if(mpi->mpiid == 0)
    {
      std::printf("Saving \"%s\"\n", filename);
      FILE *pFile;
      pFile = fopen(filename, "wb");

      if(pFile == NULL)
      {
        std::printf("ERROR \"%s\" cannot be written", filename);
        return 1;
      }

      fwrite(&ug[grid->kstart], sizeof(double), grid->kmax, pFile);
      fwrite(&vg[grid->kstart], sizeof(double), grid->kmax, pFile);
      
      fclose(pFile);
    }
  }

  return 0;
}

int cforce::load()
{
  int nerror = 0;

  if(swforce == "2")
  {
    char filename[256];
    std::sprintf(filename, "%s.%07d", "force", 0);

    if(mpi->mpiid == 0)
    {
      std::printf("Loading \"%s\"\n", filename);

      FILE *pFile;
      pFile = fopen(filename, "rb");

      if(pFile == NULL)
      {
        std::printf("ERROR \"%s\" does not exist\n", filename);
        ++nerror;
      }
      else
      {
        fread(&ug[grid->kstart], sizeof(double), grid->kmax, pFile);
        fread(&vg[grid->kstart], sizeof(double), grid->kmax, pFile);
      
        fclose(pFile);
      }
    }

    mpi->broadcast(&nerror, 1);
    if(nerror)
      return 1;

    // send the buffers to all processes
    mpi->broadcast(ug, grid->kmax);
    mpi->broadcast(vg, grid->kmax);
  }

  return 0;
}

int cforce::flux(double * restrict ut, double * restrict u, double * restrict dz, double dt)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  
  double uavg, utavg;

  uavg  = 0.;
  utavg = 0.;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        uavg  = uavg  + u [ijk]*dz[k];
        utavg = utavg + ut[ijk]*dz[k];
      }

  grid->getsum(&uavg);
  grid->getsum(&utavg);

  uavg  = uavg  / (grid->itot*grid->jtot*grid->zsize);
  utavg = utavg / (grid->itot*grid->jtot*grid->zsize);

  double fbody; 
  fbody = (uflow - uavg) / dt - utavg;

  for(int n=0; n<grid->ncells; n++)
    ut[n] += fbody;

  return 0;
}

int cforce::coriolis_2nd(double * restrict ut, double * restrict vt,
                         double * restrict u , double * restrict v ,
                         double * restrict ug, double * restrict vg)
{
  int ijk,ii,jj,kk;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ut[ijk] += fc * (0.25*(v[ijk-ii] + v[ijk] + v[ijk-ii+jj] + v[ijk+jj]) - vg[k]);
      }

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        vt[ijk] -= fc * (0.25*(u[ijk-jj] + u[ijk] + u[ijk+ii-jj] + u[ijk+ii]) - ug[k]);
      }

  return 0;
}

int cforce::coriolis_4th(double * restrict ut, double * restrict vt,
                         double * restrict u , double * restrict v ,
                         double * restrict ug, double * restrict vg)
{
  int ijk,ii1,ii2,jj1,jj2,kk1;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        ut[ijk] += fc * ( ( ci0*(ci0*v[ijk-ii2-jj1] + ci1*v[ijk-ii1-jj1] + ci2*v[ijk-jj1] + ci3*v[ijk+ii1-jj1])
                          + ci1*(ci0*v[ijk-ii2    ] + ci1*v[ijk-ii1    ] + ci2*v[ijk    ] + ci3*v[ijk+ii1    ])
                          + ci2*(ci0*v[ijk-ii2+jj1] + ci1*v[ijk-ii1+jj1] + ci2*v[ijk+jj1] + ci3*v[ijk+ii1+jj1])
                          + ci3*(ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii1+jj2] + ci2*v[ijk+jj2] + ci3*v[ijk+ii1+jj2]) ) - vg[k]);
      }

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        vt[ijk] -= fc * ( ( ci0*(ci0*u[ijk-ii1-jj2] + ci1*u[ijk-jj2] + ci2*u[ijk+ii1-jj2] + ci3*u[ijk+ii2-jj2])
                          + ci1*(ci0*u[ijk-ii1-jj1] + ci1*u[ijk-jj1] + ci2*u[ijk+ii1-jj1] + ci3*u[ijk+ii2-jj1])
                          + ci2*(ci0*u[ijk-ii1    ] + ci1*u[ijk    ] + ci2*u[ijk+ii1    ] + ci3*u[ijk+ii2    ])
                          + ci3*(ci0*u[ijk-ii1+jj1] + ci1*u[ijk+jj1] + ci2*u[ijk+ii1+jj1] + ci3*u[ijk+ii2+jj1]) ) - ug[k]);
      }

  return 0;
}

