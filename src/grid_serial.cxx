#ifndef PARALLEL
#include <fftw3.h>
#include <cstdio>
#include "grid.h"
#include "defines.h"

// MPI functions
int cgrid::initmpi()
{
  mpitypes = true;
  return 0;
} 

int cgrid::exitmpi()
{
  return 0;
}

int cgrid::boundary_cyclic(double * restrict data)
{
  int ncount = 1;
  int ijk0,ijk1,jj,kk;

  jj = icells;
  kk = icells*jcells;

  // first, east west boundaries
  for(int k=0; k<kcells; k++)
    for(int j=0; j<jcells; j++)
#pragma ivdep
      for(int i=0; i<igc; i++)
      {
        ijk0 = i          + j*jj + k*kk;
        ijk1 = iend-igc+i + j*jj + k*kk;
        data[ijk0] = data[ijk1];
      }

  for(int k=0; k<kcells; k++)
    for(int j=0; j<jcells; j++)
#pragma ivdep
      for(int i=0; i<igc; i++)
      {
        ijk0 = i+iend   + j*jj + k*kk;
        ijk1 = i+istart + j*jj + k*kk;
        data[ijk0] = data[ijk1];
      }

  // if the run is 3D, apply the BCs
  if(jtot > 1)
  {
    // second, send and receive the ghost cells in the north-south direction
    for(int k=0; k<kcells; k++)
      for(int j=0; j<jgc; j++)
#pragma ivdep
        for(int i=0; i<icells; i++)
        {
          ijk0 = i + j           *jj + k*kk;
          ijk1 = i + (jend-jgc+j)*jj + k*kk;
          data[ijk0] = data[ijk1];
        }

    for(int k=0; k<kcells; k++)
      for(int j=0; j<jgc; j++)
#pragma ivdep
        for(int i=0; i<icells; i++)
        {
          ijk0 = i + (j+jend  )*jj + k*kk;
          ijk1 = i + (j+jstart)*jj + k*kk;
          data[ijk0] = data[ijk1];
        }
  }
  // in case of 2D, fill all the ghost cells with the current value
  else
  {
    // 2d essential variables
    int ijkref,ijknorth,ijksouth,jj,kk;

    jj = icells;
    kk = icells*jcells;

    for(int k=kstart; k<kend; k++)
      for(int j=0; j<jgc; j++)
#pragma ivdep
        for(int i=istart; i<iend; i++)
        {
          ijkref   = i + jstart*jj   + k*kk;
          ijknorth = i + j*jj        + k*kk;
          ijksouth = i + (jend+j)*jj + k*kk;
          data[ijknorth] = data[ijkref];
          data[ijksouth] = data[ijkref];
        }
  }

  return 0;
}

int cgrid::boundary_cyclic2d(double * restrict data)
{
  int ij0,ij1,jj;

  jj = icells;

  // first, east west boundaries
  for(int j=0; j<jcells; j++)
#pragma ivdep
    for(int i=0; i<igc; i++)
    {
      ij0 = i          + j*jj;
      ij1 = iend-igc+i + j*jj;
      data[ij0] = data[ij1];
    }

  for(int j=0; j<jcells; j++)
#pragma ivdep
    for(int i=0; i<igc; i++)
    {
      ij0 = i+iend   + j*jj;
      ij1 = i+istart + j*jj;
      data[ij0] = data[ij1];
    }

  // if the run is 3D, apply the BCs
  if(jtot > 1)
  {
    // second, send and receive the ghost cells in the north-south direction
    for(int j=0; j<jgc; j++)
#pragma ivdep
      for(int i=0; i<icells; i++)
      {
        ij0 = i + j           *jj;
        ij1 = i + (jend-jgc+j)*jj;
        data[ij0] = data[ij1];
      }

    for(int j=0; j<jgc; j++)
#pragma ivdep
      for(int i=0; i<icells; i++)
      {
        ij0 = i + (j+jend  )*jj;
        ij1 = i + (j+jstart)*jj;
        data[ij0] = data[ij1];
      }
  }
  // in case of 2D, fill all the ghost cells with the current value
  else
  {
    // 2d essential variables
    int ijref,ijnorth,ijsouth,jj;

    jj = icells;

    for(int j=0; j<jgc; j++)
#pragma ivdep
      for(int i=istart; i<iend; i++)
      {
        ijref   = i + jstart*jj;
        ijnorth = i + j*jj;
        ijsouth = i + (jend+j)*jj;
        data[ijnorth] = data[ijref];
        data[ijsouth] = data[ijref];
      }
  }

  return 0;
}

int cgrid::transposezx(double * restrict ar, double * restrict as)
{
  int ijk,jj,kk;
  jj = imax;
  kk = imax*jmax;

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk = i + j*jj + k*kk;
        ar[ijk] = as[ijk];
      }

  return 0;
}

int cgrid::transposexz(double * restrict ar, double * restrict as)
{
  int ijk,jj,kk;
  jj = imax;
  kk = imax*jmax;

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk = i + j*jj + k*kk;
        ar[ijk] = as[ijk];
      }

  return 0;
}

int cgrid::transposexy(double * restrict ar, double * restrict as)
{
  int ijk,jj,kk;
  jj = imax;
  kk = imax*jmax;

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk = i + j*jj + k*kk;
        ar[ijk] = as[ijk];
      }

  return 0;
}

int cgrid::transposeyx(double * restrict ar, double * restrict as)
{
  int ijk,jj,kk;
  jj = imax;
  kk = imax*jmax;

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk = i + j*jj + k*kk;
        ar[ijk] = as[ijk];
      }

  return 0;
}

int cgrid::transposeyz(double * restrict ar, double * restrict as)
{
  int ijk,jj,kk;
  jj = imax;
  kk = imax*jmax;

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk = i + j*jj + k*kk;
        ar[ijk] = as[ijk];
      }

  return 0;
}

int cgrid::transposezy(double * restrict ar, double * restrict as)
{
  int ijk,jj,kk;
  jj = imax;
  kk = imax*jmax;

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk = i + j*jj + k*kk;
        ar[ijk] = as[ijk];
      }

  return 0;
}

int cgrid::getmax(double *var)
{
  return 0;
}

int cgrid::getsum(double *var)
{
  return 0;
}

int cgrid::getprof(double *prof, int kcellsin)
{
  return 0;
}

// IO functions
int cgrid::save()
{
  // SAVE THE GRID
  FILE *pFile;
  char filename[256];
  std::sprintf(filename, "%s.%07d", "grid", 0);
  pFile = fopen(filename, "wb");

  if(pFile == NULL)
  {
    std::printf("ERROR \"%s\" cannot be written\n", filename);
    return 1;
  }
  else
    std::printf("Saving \"%s\"\n", filename);

  fwrite(&x [istart], sizeof(double), itot, pFile);
  fwrite(&xh[istart], sizeof(double), itot, pFile);
  fwrite(&y [jstart], sizeof(double), jtot, pFile);
  fwrite(&yh[jstart], sizeof(double), jtot, pFile);
  fwrite(&z [kstart], sizeof(double), ktot, pFile);
  fwrite(&zh[kstart], sizeof(double), ktot, pFile);
  fclose(pFile);

  // SAVE THE FFTW PLAN
  iplanf = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_R2HC, FFTW_EXHAUSTIVE);
  iplanb = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_HC2R, FFTW_EXHAUSTIVE);
  jplanf = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_R2HC, FFTW_EXHAUSTIVE);
  jplanb = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_HC2R, FFTW_EXHAUSTIVE);

  fftwplan = true;

  if(mpi->mpiid == 0)
  {
    char filename[256];
    std::sprintf(filename, "%s.%07d", "fftwplan", 0);

    std::printf("Saving \"%s\"\n", filename);

    int n = fftw_export_wisdom_to_filename(filename);
    if(n == 0)
    {
      std::printf("ERROR \"%s\" cannot be saved\n", filename);
      return 1;
    }
  }

  return 0;
}

int cgrid::load()
{
  // LOAD THE GRID
  FILE *pFile;
  char filename[256];
  std::sprintf(filename, "%s.%07d", "grid", 0);
  pFile = fopen(filename, "rb");

  if(pFile == NULL)
  {
    std::printf("ERROR \"%s\" does not exist\n", filename);
    return 1;
  }
  else
    std::printf("Loading \"%s\"\n", filename);

  fread(&x [istart], sizeof(double), itot, pFile);
  fread(&xh[istart], sizeof(double), itot, pFile);
  fread(&y [jstart], sizeof(double), jtot, pFile);
  fread(&yh[jstart], sizeof(double), jtot, pFile);
  fread(&z [kstart], sizeof(double), ktot, pFile);
  fread(&zh[kstart], sizeof(double), ktot, pFile);
  fclose(pFile);

  // calculate the missing coordinates
  calculate();

  // LOAD THE FFTW PLAN
  std::sprintf(filename, "%s.%07d", "fftwplan", 0);

  if(mpi->mpiid == 0)
    std::printf("Loading \"%s\"\n", filename);

  int n = fftw_import_wisdom_from_filename(filename);
  if(n == 0)
  {
    if(mpi->mpiid == 0)
      std::printf("ERROR \"%s\" does not exist\n", filename);
    return 1;
  }

  iplanf = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_R2HC, FFTW_EXHAUSTIVE);
  iplanb = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_HC2R, FFTW_EXHAUSTIVE);
  jplanf = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_R2HC, FFTW_EXHAUSTIVE);
  jplanb = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_HC2R, FFTW_EXHAUSTIVE);

  fftwplan = true;

  fftw_forget_wisdom();

  return 0;
}

int cgrid::savefield3d(double * restrict data, double * restrict tmp1, double * restrict tmp2, char *filename)
{
  FILE *pFile;
  pFile = fopen(filename, "wb");

  if(pFile == NULL)
    return 1;

  int ijk,jj,kk;

  jj = icells;
  kk = icells*jcells;

  for(int k=kstart; k<kend; k++)
    for(int j=jstart; j<jend; j++)
      {
        ijk = istart + j*jj + k*kk;
        fwrite(&data[ijk], sizeof(double), imax, pFile);
      }

  fclose(pFile);

  return 0;
}

int cgrid::loadfield3d(double * restrict data, double * restrict tmp1, double * restrict tmp2, char *filename)
{
  FILE *pFile;
  pFile = fopen(filename, "rb");

  if(pFile == NULL)
    return 1;

  int ijk,jj,kk;

  jj = icells;
  kk = icells*jcells;

  for(int k=kstart; k<kend; k++)
    for(int j=jstart; j<jend; j++)
      {
        ijk = istart + j*jj + k*kk;
        fread(&data[ijk], sizeof(double), imax, pFile);
      }

  fclose(pFile);

  return 0;
}

int cgrid::fftforward(double * restrict data,   double * restrict tmp1,
                      double * restrict fftini, double * restrict fftouti,
                      double * restrict fftinj, double * restrict fftoutj)
{
  int ijk,jj,kk;

  jj = itot;
  kk = itot*jmax;

  // do the first fourier transform
  for(int k=0; k<kblock; k++)
    for(int j=0; j<jmax; j++)
    {
#pragma ivdep
      for(int i=0; i<itot; i++)
      {
        ijk = i + j*jj + k*kk;
        fftini[i] = data[ijk];
      }

      fftw_execute(iplanf);

#pragma ivdep
      for(int i=0; i<itot; i++)
      {
        ijk = i + j*jj + k*kk;
        data[ijk] = fftouti[i];
      }
    }

  jj = iblock;
  kk = iblock*jtot;

  // do the second fourier transform
  for(int k=0; k<kblock; k++)
    for(int i=0; i<iblock; i++)
    {
      for(int j=0; j<jtot; j++)
      {
        ijk = i + j*jj + k*kk;
        fftinj[j] = data[ijk];
      }

      fftw_execute(jplanf);

      for(int j=0; j<jtot; j++)
      {
        ijk = i + j*jj + k*kk;
        // shift to use p in pressure solver
        data[ijk] = fftoutj[j];
      }
    }

  return 0;
}

int cgrid::fftbackward(double * restrict data,   double * restrict tmp1,
                       double * restrict fftini, double * restrict fftouti,
                       double * restrict fftinj, double * restrict fftoutj)
{
  int ijk,jj,kk;

  jj = iblock;
  kk = iblock*jtot;

  // transform the second transform back
  for(int k=0; k<kblock; k++)
    for(int i=0; i<iblock; i++)
    {
      for(int j=0; j<jtot; j++)
      {
        ijk = i + j*jj + k*kk;
        fftinj[j] = data[ijk];
      }

      fftw_execute(jplanb);

      for(int j=0; j<jtot; j++)
      {
        ijk = i + j*jj + k*kk;
        data[ijk] = fftoutj[j] / jtot;
      }
    }

  jj = itot;
  kk = itot*jmax;

  // transform the first transform back
  for(int k=0; k<kblock; k++)
    for(int j=0; j<jmax; j++)
    {
#pragma ivdep
      for(int i=0; i<itot; i++)
      {
        ijk = i + j*jj + k*kk;
        fftini[i] = data[ijk];
      }

      fftw_execute(iplanb);

#pragma ivdep
      for(int i=0; i<itot; i++)
      {
        ijk = i + j*jj + k*kk;
        // swap array here to avoid unncessary 3d loop
        tmp1[ijk] = fftouti[i] / itot;
      }
    }

  return 0;
}

int cgrid::savexzslice(double * restrict data, double * restrict tmp, int jslice, char *filename)
{
  // extract the data from the 3d field without the ghost cells
  int ijk,jj,kk;
  int ijkb,jjb,kkb;

  jj  = icells;
  kk  = icells*jcells;
  jjb = imax;
  kkb = imax;

  int count = imax*kmax;

  for(int k=0; k<kmax; k++)
#pragma ivdep
    for(int i=0; i<imax; i++)
    {
      // take the modulus of jslice and jmax to have the right offset within proc
      ijk  = i+igc + (jslice+jgc)*jj + (k+kgc)*kk;
      ijkb = i + k*kkb;
      tmp[ijkb] = data[ijk];
    }

  FILE *pFile;
  pFile = fopen(filename, "wb");
  if(pFile == NULL)
    return 1;

  fwrite(tmp, sizeof(double), count, pFile);
  fclose(pFile);

  return 0;
}

int cgrid::savexyslice(double * restrict data, double * restrict tmp, int kslice, char *filename)
{
  // extract the data from the 3d field without the ghost cells
  int ijk,jj,kk;
  int ijkb,jjb,kkb;

  jj  = icells;
  kk  = icells*jcells;
  jjb = imax;

  int count = imax*jmax;

  for(int j=0; j<jmax; j++)
#pragma ivdep
    for(int i=0; i<imax; i++)
    {
      // take the modulus of jslice and jmax to have the right offset within proc
      ijk  = i+igc + (j+jgc)*jj + (kslice+kgc)*kk;
      ijkb = i + j*jjb;
      tmp[ijkb] = data[ijk];
    }

  FILE *pFile;
  pFile = fopen(filename, "wb");
  if(pFile == NULL)
    return 1;

  fwrite(tmp, sizeof(double), count, pFile);
  fclose(pFile);

  return 0;
}
#endif
